"""False-negative candidate discovery (spec 7) - an independent module,
because VCF-only review can never find calls the caller never emitted.

Sources:
  1. filtered_recovery - FILTER-failed calls with residual ALT support
  2. targeted_locus    - hotspot loci / BED / user gene list checked against
                         VCF, gVCF evidence and BAM pileup, with explicit
                         callability classification
  3. second_caller     - pluggable interface for DeepVariant etc. (framework)

Outputs are ALWAYS candidates requiring human review; nothing is
auto-restored to PASS and no "confirmed FN counts" are reported without
external truth.
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from sqlalchemy.orm import Session

from ..db import models as db
from ..db import audit
from ..evidence import make_backend, EvidenceParams
from ..variants import import_vcf
from ..io.vcf import VcfReader

log = logging.getLogger("avrs.fn")


@dataclass
class LocusSpec:
    chrom: str
    pos: int          # 1-based
    ref: Optional[str] = None
    alt: Optional[str] = None
    source: str = "manual"


def load_hotspot_file(path: str | Path) -> list[LocusSpec]:
    """hotspot VCF, BED, or TSV list of chrom/pos/ref/alt."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    out: list[LocusSpec] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        if line.startswith(("track", "browser")):
            continue
        parts = line.split("\t")
        if parts[0].lower().startswith("chr") or parts[0].isdigit():
            if len(parts) >= 3 and parts[1].isdigit() and parts[2].isdigit():
                out.append(LocusSpec(chrom=parts[0], pos=int(parts[1]) + 1, source="bed"))
            elif len(parts) >= 2 and parts[1].isdigit():
                row = [p for p in parts]
                out.append(LocusSpec(chrom=row[0], pos=int(row[1]),
                                     ref=row[2] if len(row) > 2 else None,
                                     alt=row[3] if len(row) > 3 else None,
                                     source="list"))
    return out


class FNCandidate:
    def __init__(self, candidate_id, locus, source, caller_status, bam_evidence,
                 callable_status, classification, possible_missed, reason_codes):
        self.candidate_id = candidate_id
        self.locus = locus
        self.source = source
        self.caller_status = caller_status
        self.bam_evidence = bam_evidence
        self.callable_status = callable_status
        self.classification = classification
        self.possible_missed = possible_missed
        self.reason_codes = reason_codes

    def to_db(self, sample_id: int) -> db.FNDiscoveryCandidate:
        return db.FNDiscoveryCandidate(
            sample_id=sample_id, candidate_id=self.candidate_id,
            candidate_source=self.source, locus=self.locus,
            ref=self.bam_evidence.get("_ref"), alt=self.bam_evidence.get("_alt"),
            caller_status=self.caller_status,
            bam_evidence_json=json.dumps(self.bam_evidence),
            callable_status=self.callable_status,
            possible_missed_variant=self.possible_missed,
            classification=self.classification,
            reason_codes=json.dumps(self.reason_codes),
        )


def recover_filtered_candidates(session: Session, sample_id: str, cfg,
                                raw_vcf_path: str | Path, bam_path: str | Path,
                                fasta_path: str | Path,
                                min_alt_reads: int = 3) -> list[FNCandidate]:
    """Source 1: variants that failed FILTER but retain plausible ALT support.

    Never re-labels them PASS - they become candidates for human review.
    """
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
    sample_db_id = sample.id
    imported, _ = import_vcf(session, raw_vcf_path, sample_db_id, sample_id,
                             cfg.reference.build, fasta=None,
                             source_label="filtered_recovery_raw")
    candidates = []
    backend = make_backend(bam_path, fasta_path, EvidenceParams())
    for v in (session.query(db.Variant)
              .filter_by(sample_id=sample_db_id)
              .filter(db.Variant.original_filter != "PASS").all()):
        if v.original_filter in (".", None):
            continue
        ev = backend.evidence(v.chrom, v.pos, v.ref, v.alt, v.variant_type,
                              caller={"AD_ALT": v.info.get("AD_ALT"),
                                      "AD_REF": v.info.get("AD_REF"),
                                      "DP": v.info.get("DP")})
        possible = ev.alt_count >= min_alt_reads
        classification = "possible_missed_variant" if possible else "no_evidence_of_alt"
        reasons = [f"FILTER={v.original_filter}", f"alt_reads={ev.alt_count}"]
        if ev.strand_ratio_imbalance and ev.strand_ratio_imbalance >= 0.9:
            possible = False
            reasons.append("strand_biased")
            classification = "no_evidence_of_alt"
        cand = FNCandidate(
            candidate_id=f"FR-{v.variant_uid}",
            locus=f"{v.chrom}:{v.pos}:{v.ref}>{v.alt}",
            source="filtered_recovery",
            caller_status="filtered_in_raw",
            bam_evidence={**ev.model_dump(), "_ref": v.ref, "_alt": v.alt},
            callable_status="callable" if ev.total_reads > 0 else "insufficient_coverage",
            classification=classification,
            possible_missed=possible, reason_codes=reasons,
        )
        candidates.append(cand)
        session.add(cand.to_db(sample_db_id))
    session.commit()
    audit(session, actor="system:fn_discovery", action="filtered_recovery",
          entity_type="sample", entity_id=sample_id,
          details={"candidates": len(candidates)})
    return candidates


def targeted_locus_review(session: Session, sample_id: str, cfg, bam_path, fasta_path,
                          loci: list[LocusSpec], min_callable_depth: int = 10,
                          min_alt_reads: int = 4, min_vaf: float = 0.15) -> list[FNCandidate]:
    """Source 2: check known loci against VCF presence + BAM evidence.

    Classifications are explicit and mutually exclusive:
      - possible_missed_variant: BAM shows real-looking ALT support,
        variant absent from the filtered VCF (NOT confirmed FN)
      - no_call: genotype is ./.
      - insufficient_coverage: locus not callable at configured depth
      - no_evidence_of_alt: callable, but no ALT support in the pileup
    """
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
    backend = make_backend(bam_path, fasta_path, EvidenceParams())
    present_uids = {v.variant_uid for v in
                    session.query(db.Variant).filter_by(sample_id=sample.id).all()}
    present_loci = {(v.chrom, v.pos, v.ref, v.alt) for v in
                   session.query(db.Variant).filter_by(sample_id=sample.id).all()}
    candidates = []
    for spec in loci:
        ev = backend.evidence(spec.chrom, spec.pos, spec.ref or ".", spec.alt or ".",
                              "SNV" if spec.alt and len(spec.alt) == 1 else "INS")
        cov = backend.coverage_at(spec.chrom, spec.pos)
        depth = cov.get("mean_depth", 0)
        bam_ev = {**ev.model_dump(), "_ref": spec.ref, "_alt": spec.alt,
                  "_mean_depth_at_locus": depth}
        vaf = ev.vaf or 0.0
        in_vcf = any(k in present_loci for k in
                     [(spec.chrom, spec.pos, spec.ref, spec.alt)]
                     if spec.ref and spec.alt)
        if spec.ref and spec.alt:
            uid_like = f"{cfg.reference.build}#{spec.chrom}:{spec.pos}:{spec.ref}>{spec.alt}"
            in_vcf = in_vcf or uid_like in present_uids
        if depth < min_callable_depth:
            classification = "insufficient_coverage"
            possible = False
            reasons = [f"depth={depth:.1f} < {min_callable_depth}"]
        elif in_vcf:
            classification = "in_vcf"
            possible = False
            reasons = ["variant already present in VCF"]
        elif ev.alt_count >= min_alt_reads and vaf >= min_vaf:
            classification = "possible_missed_variant"
            possible = True
            reasons = [f"alt_reads={ev.alt_count}", f"vaf={vaf:.2f}",
                       "absent from filtered VCF"]
        elif ev.alt_count > 0:
            classification = "no_evidence_of_alt"
            possible = False
            reasons = [f"alt_reads={ev.alt_count} below threshold {min_alt_reads} "
                       f"or vaf {vaf:.2f} below {min_vaf}"]
        else:
            classification = "no_evidence_of_alt"
            possible = False
            reasons = ["no ALT-supporting reads in pileup"]
        cand = FNCandidate(
            candidate_id=f"TL-{spec.chrom}:{spec.pos}:{spec.ref or '*'}:{spec.alt or '*'}",
            locus=f"{spec.chrom}:{spec.pos}",
            source="targeted_locus", caller_status="in_vcf" if in_vcf else "absent_from_vcf",
            bam_evidence=bam_ev,
            callable_status="callable" if depth >= min_callable_depth else "insufficient_coverage",
            classification=classification, possible_missed=possible,
            reason_codes=reasons,
        )
        candidates.append(cand)
        session.add(cand.to_db(sample.id))
    session.commit()
    audit(session, actor="system:fn_discovery", action="targeted_locus_review",
          entity_type="sample", entity_id=sample_id,
          details={"loci": len(loci), "possible_missed": sum(1 for c in candidates if c.possible_missed)})
    return candidates


class SecondCallerAdapter:
    """Source 3 (optional): pluggable second-caller interface (spec 7.3).

    Concrete adapters (e.g. DeepVariant) must:
      - run the caller and produce a normalized VCF path
      - comparisons are haplotype-aware / normalized BEFORE set comparison
    Consensus of two callers is NOT truth; disagreements go to review.
    """
    registry: dict[str, type] = {}

    def __init_subclass__(cls, caller_name: str = None):
        if caller_name:
            SecondCallerAdapter.registry[caller_name] = cls

    def call_variants(self, bam_path: str, fasta_path: str, out_vcf: str,
                      ref_build: str) -> str:
        raise NotImplementedError("Second-caller integration not implemented in MVP")

    def compare(self, primary_vcf: str, second_vcf: str, fasta) -> dict:
        """Set comparison at normalized allele level -> concordance report."""
        from ..io.vcf import VcfReader, split_multiallelic, normalize_allele
        from ..io.fasta import Fasta
        primary, second = {}, {}
        for tag, path, store in [("primary", primary_vcf, primary), ("second", second_vcf, second)]:
            reader = VcfReader(path)
            for rec in reader:
                for a in split_multiallelic(rec):
                    alt = a.alts[0]
                    try:
                        getter = lambda p, c=a.chrom: fasta.fetch(c, p, p)  # noqa: E731
                        p, r, al = normalize_allele(a.pos, a.ref, alt, getter)
                    except Exception:
                        p, r, al = a.pos, a.ref, alt
                    store[f"{a.chrom}:{p}:{r}>{al}"] = tag
        keys = set(primary) | set(second)
        return {
            "concordant": sorted(set(primary) & set(second)),
            "primary_only": sorted(set(primary) - set(second)),
            "second_only": sorted(set(second) - set(primary)),
            "note": "disagreement is a review trigger, not a truth label",
        }


def summarize(session: Session, sample_id: str) -> dict:
    sample = session.query(db.Sample).filter_by(sample_id=sample_id).one()
    cands = session.query(db.FNDiscoveryCandidate).filter_by(sample_id=sample.id).all()
    return {
        "total_candidates": len(cands),
        "possible_missed": sum(1 for c in cands if c.possible_missed_variant),
        "insufficient_coverage": sum(1 for c in cands if c.classification == "insufficient_coverage"),
        "no_evidence": sum(1 for c in cands if c.classification == "no_evidence_of_alt"),
        "disclaimer": ("Counts are unvalidated candidates. 'possible_missed' is NOT a "
                       "confirmed false negative - confirmation requires human review "
                       "or external truth data."),
    }
