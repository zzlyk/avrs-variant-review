"""VCF import: parse, split multiallelics, normalize (allele level),
assign stable build-scoped variant IDs, and persist with the original
representation mapping (spec 3.5)."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from sqlalchemy.orm import Session

from ..db import models as db
from ..db import audit
from ..io.vcf import VcfReader, split_multiallelic, normalize_allele, variant_uid, classify_variant


class _RefGetter:
    def __init__(self, fasta):
        self.fasta = fasta

    def base_at(self, chrom: str):
        def get(pos1: int) -> str:
            return self.fasta.fetch(chrom, pos1, pos1)
        return get


def parse_ad(rec) -> tuple[Optional[int], Optional[int]]:
    """Caller-reported AD for the (first) ALT of this record."""
    for fmt in rec.formats.values():
        ad = fmt.get("AD")
        if ad and "," in ad:
            parts = ad.split(",")
            try:
                r = int(parts[0])
                a = int(parts[1]) if len(parts) > 1 else None
                return r, a
            except ValueError:
                return None, None
        break
    return None, None


def import_vcf(session: Session, vcf_path: str | Path, sample_db_id: int,
               sample_id: str, reference_build: str, fasta=None,
               chr_prefix_policy: str = "strict",
               source_label: str = "raw_vcf") -> tuple[int, int]:
    """Import a VCF into the variants table.

    Returns (imported_count, skipped_count).
    chr_prefix_policy:
      - strict: VCF contig names must match the reference contigs exactly
        (no blind chr1<->1 rewriting; mismatches raise).
      - suffix_match: accept 'chr'-less names if the reference uses the
        'chr' prefix (or vice versa) and record the mapping explicitly.
    """
    vcf_path = Path(vcf_path)
    reader = VcfReader(vcf_path)
    ref_names = set()
    fasta_map = {}
    if fasta is not None:
        ref_names = set(fasta.contig_names)
    else:
        ref_names = set(reader.header.contigs)

    def _map_chrom(c: str) -> str:
        if c in ref_names:
            return c
        if chr_prefix_policy == "suffix_match":
            if f"chr{c}" in ref_names:
                return f"chr{c}"
            if c.startswith("chr") and c[3:] in ref_names:
                return c[3:]
        raise ValueError(
            f"Chromosome {c!r} from {vcf_path.name} not found in reference "
            f"contigs. Reference/VCF chromosome naming mismatch must be fixed "
            f"upstream (e.g. via bcftools annotate/reheader) - AVRS never "
            f"blind-rewrites contig names.")

    imported = 0
    skipped = 0
    # idempotent import: re-running against the same sample must not violate
    # the (sample_id, variant_uid) uniqueness constraint
    existing = {uid for (uid,) in session.query(db.Variant.variant_uid)
                .filter_by(sample_id=sample_db_id).all()}
    for rec in reader:
        # keep even filtered records: FN discovery rescues them later
        for allele_rec in split_multiallelic(rec):
            alt = allele_rec.alts[0] if allele_rec.alts else None
            if not alt:
                skipped += 1
                continue
            if alt == "*":
                skipped += 1
                continue
            chrom = _map_chrom(allele_rec.chrom)
            npos, nref, nalt = allele_rec.pos, allele_rec.ref, alt
            if fasta is not None:
                getter = _RefGetter(fasta).base_at(chrom)
                try:
                    if fasta.fetch(chrom, allele_rec.pos, allele_rec.pos + len(allele_rec.ref) - 1).upper() != allele_rec.ref.upper():
                        raise ValueError(
                            f"REF mismatch at {chrom}:{allele_rec.pos}: VCF REF "
                            f"{allele_rec.ref} != reference "
                            f"{fasta.fetch(chrom, allele_rec.pos, allele_rec.pos + len(allele_rec.ref) - 1)}")
                    npos, nref, nalt = normalize_allele(allele_rec.pos, allele_rec.ref, alt, getter)
                except KeyError:
                    pass  # contig missing from FASTA: keep original repr
            uid = variant_uid(reference_build, chrom, npos, nref, nalt)
            if uid in existing:
                skipped += 1
                continue
            ad_r, ad_a = parse_ad(allele_rec)
            gt = None
            for fmt in allele_rec.formats.values():
                gt = fmt.get("GT")
                break
            info = {
                "DP": allele_rec.info_int("DP"),
                "AD_REF": ad_r,
                "AD_ALT": ad_a,
                "GQ": _first_int(allele_rec, "GQ"),
                "FILTER": allele_rec.filt,
            }
            v = db.Variant(
                variant_uid=uid,
                sample_id=sample_db_id,
                original_chrom=allele_rec.chrom,
                original_pos=allele_rec.pos,
                original_ref=allele_rec.ref,
                original_alt=alt,
                original_filter=allele_rec.filt,
                chrom=chrom,
                pos=npos,
                ref=nref,
                alt=nalt,
                variant_type=classify_variant(nref, nalt),
                genotype=gt,
                qual=allele_rec.qual,
                info_json=__import__("json").dumps(info),
            )
            session.add(v)
            existing.add(uid)
            imported += 1
    session.flush()
    audit(session, actor="system:variants.importer", action="import_vcf",
          entity_type="sample", entity_id=sample_id,
          details={"path": str(vcf_path), "imported": imported, "skipped": skipped,
                   "source": source_label})
    return imported, skipped


def _first_int(rec, key: str):
    for fmt in rec.formats.values():
        v = fmt.get(key)
        if v is not None:
            try:
                return int(float(v))
            except ValueError:
                return None
    return None
