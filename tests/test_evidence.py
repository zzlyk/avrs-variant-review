"""Unit tests: BAM io round-trip, evidence extraction (SNV + indel), rules."""
import pytest

from avrs.io.bam import BamHeader, BamRecord, iter_records, fetch_region, write_bam
from avrs.evidence import PurePythonEvidenceBackend
from avrs.evidence.core import ReadView, compute_evidence, EvidenceParams
from avrs.io.fasta import write_fasta, write_fai, Fasta
from avrs.triage import evaluate_rules, prioritize, should_review, suspiciousness_score
from avrs.config import TriageProfile


def _ref_seq(contig_len=2000):
    """Single source of truth for the synthetic reference used by the tests."""
    return "".join("ACGT"[(i * 7) % 4] for i in range(contig_len))


def _rec(name, pos0, seq, flag=0, mapq=60, cigar=None, qual=None):
    cigar = cigar or [(len(seq), "M")]
    return BamRecord(qname=name, flag=flag, ref_id=0, pos=pos0, mapq=mapq,
                     cigar=cigar, seq=seq, qual=qual or [37] * len(seq))


def _write_demo_bam(tmp_path, recs, contig_len=2000):
    fasta = tmp_path / "ref.fa"
    seq = _ref_seq(contig_len)
    write_fasta(fasta, {"chr1": seq})
    write_fai(fasta, {"chr1": seq})
    header = BamHeader(text="@HD\tVN:1.6\tSO:coordinate\n",
                       refs=[("chr1", contig_len)])
    bam = tmp_path / "t.bam"
    write_bam(bam, header, recs)
    return bam, fasta, seq


def test_bam_roundtrip_and_index(tmp_path):
    recs = [_rec(f"r{i}", 100 + i * 10, "ACGTACGTAC") for i in range(5)]
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    assert (tmp_path / "t.bam.bai").exists()
    got = list(iter_records(bam))
    assert len(got) == 5
    assert got[0][0].refs[0][0] == "chr1"
    assert got[0][1].qname == "r0"
    # indexed region fetch must agree with a sequential scan
    region = fetch_region(bam, "chr1", 100, 140)
    scanned = [r for _, r in iter_records(bam)
               if r.pos < 140 and r.pos + r.reference_length() > 100]
    assert len(region) == len(scanned) >= 2
    assert {r.qname for r in region} == {r.qname for r in scanned}


def test_bam_index_matches_scan_for_many_records(tmp_path):
    recs = [_rec(f"r{i}", 100 + i * 7, "ACGTACGTACGTACGTACGT") for i in range(200)]
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    for start in (100, 500, 900, 1200):
        end = start + 100
        region = fetch_region(bam, "chr1", start, end)
        scanned = [r for _, r in iter_records(bam)
                   if r.pos < end and r.pos + r.reference_length() > start]
        assert {r.qname for r in region} == {r.qname for r in scanned}, f"mismatch at {start}"


def test_bam_qualities_and_cigar_parse(tmp_path):
    recs = [_rec("r1", 50, "ACGT", cigar=[(2, "M"), (2, "I"), (2, "M")], qual=[30, 30, 30, 30, 30, 30])]
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    rec = list(iter_records(bam))[0][1]
    assert rec.cigar == [(2, "M"), (2, "I"), (2, "M")]
    assert rec.reference_length() == 4
    assert rec.qual[0] == 30


def test_evidence_snv_counts(tmp_path):
    ref_seq = _ref_seq()
    pos = 500
    ref_base = ref_seq[pos - 1]
    alt_base = "A" if ref_base != "A" else "C"
    recs = []
    for i in range(10):
        s = list(ref_seq[pos - 21: pos + 21])
        if i >= 5:
            s[20] = alt_base
        recs.append(_rec(f"r{i}", pos - 21, "".join(s)))
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    backend = PurePythonEvidenceBackend(bam, fasta, EvidenceParams())
    ev = backend.evidence("chr1", pos, ref_base, alt_base, "SNV",
                           caller={"DP": 10, "AD_REF": 5, "AD_ALT": 5, "GT": "0/1"})
    assert ev.total_reads == 10
    assert ev.alt_count == 5
    assert ev.ref_count == 5
    assert abs(ev.vaf - 0.5) < 1e-6
    assert ev.strand_forward == 5 and ev.strand_reverse == 0
    assert ev.ad_consistency is True


def test_evidence_indel_uses_cigar_not_base_pileup(tmp_path):
    """A deletion must be detected through the D CIGAR op, not by base counts.

    VCF convention: POS is the anchor base, REF = anchor + deleted bases.
    Read: 10M (1-based 590-599), 2D (600-601 deleted), 30M (602-631).
    """
    ref_seq = _ref_seq()
    pos = 599           # 1-based anchor
    ref_allele = ref_seq[pos - 1: pos + 2]     # anchor + 2 deleted bases
    assert len(ref_allele) == 3
    start = pos - 10 + 1                        # 0-based 590 -> 1-based 591
    start0 = pos - 10                            # 0-based 589
    read_seq = ref_seq[start0: start0 + 10] + ref_seq[start0 + 12: start0 + 42]
    assert len(read_seq) == 40
    recs = [_rec(f"d{i}", start0, read_seq, cigar=[(10, "M"), (2, "D"), (30, "M")])
            for i in range(6)]
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    backend = PurePythonEvidenceBackend(bam, fasta, EvidenceParams())
    ev = backend.evidence("chr1", pos, ref_allele, ref_allele[0], "DEL")
    assert ev.total_reads == 6
    assert ev.alt_count == 6, "deletion support must come from CIGAR D ops"
    assert ev.ref_count == 0


def test_evidence_empty_coverage(tmp_path):
    recs = [_rec("r1", 10, "ACGTACGTAC")]
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    backend = PurePythonEvidenceBackend(bam, fasta, EvidenceParams())
    ev = backend.evidence("chr1", 1500, "A", "C", "SNV")
    assert ev.total_reads == 0
    assert ev.alt_count == 0 and ev.ref_count == 0
    assert ev.vaf is None
    assert "no_reads_overlap_variant" in ev.warnings


def test_evidence_strand_bias_detection(tmp_path):
    ref_seq = _ref_seq()
    pos = 700
    ref_base = ref_seq[pos - 1]
    alt_base = "A" if ref_base != "A" else "C"
    recs = []
    for i in range(12):
        s = list(ref_seq[pos - 21: pos + 21])
        if i >= 6:
            s[20] = alt_base
        recs.append(_rec(f"r{i}", pos - 21, "".join(s), flag=0x10 if i >= 6 else 0))
    bam, fasta, seq = _write_demo_bam(tmp_path, recs)
    backend = PurePythonEvidenceBackend(bam, fasta, EvidenceParams())
    ev = backend.evidence("chr1", pos, ref_base, alt_base, "SNV")
    assert ev.strand_forward == 0 and ev.strand_reverse == 6
    assert ev.strand_ratio_imbalance == 1.0


def test_rule_engine_flags_and_priority():
    prof = TriageProfile()
    from avrs.evidence.metrics import VariantEvidence
    ev = VariantEvidence(chrom="chr1", pos=10, ref="A", alt="C", variant_type="SNV",
                         total_reads=3, usable_reads=3, ref_count=1, alt_count=2,
                         other_count=0, vaf=0.67, strand_forward=0, strand_reverse=2,
                         strand_ratio_imbalance=1.0, mean_mapq_alt=8.0,
                         mean_baseq_alt=7.0, alt_start_cluster_frac=0.9,
                         dup_frac=0.5, caller_gq=5, ad_consistency=False,
                         genotype="0/1", vaf_expected_for_gt=0.5, vaf_deviation=0.17)
    hits = evaluate_rules(ev, prof)
    ids = {h.rule_id for h in hits}
    assert "LOW_DP" in ids
    assert "LOW_GQ" in ids
    assert "LOW_MQ" in ids
    assert "STRAND_BIAS" in ids
    assert "START_POS_CLUSTER" in ids
    assert "AD_PILEUP_MISMATCH" in ids
    assert "HIGH_DUP_RATE" in ids
    for h in hits:
        assert h.rule_id and h.reason and h.severity
    score = suspiciousness_score(hits)
    assert 0 <= score <= 100
    assert prioritize(hits) in ("HIGH", "MEDIUM", "LOW")


def test_review_modes():
    hits = []
    assert should_review("all", hits, "LOW", None, [], False, False, "u", []) is True
    assert should_review("suspicious_only", hits, "LOW", None, [], False, False, "u", []) is False
    assert should_review("suspicious_only", hits, "HIGH", None, [], False, False, "u", []) is True
    assert should_review("target_genes", hits, "LOW", "BRCA1", ["BRCA1"], False, False, "u", []) is True
    assert should_review("filtered_calls", hits, "LOW", None, [], True, False, "u", []) is True
    assert should_review("bed_regions", hits, "LOW", None, [], False, True, "u", []) is True
    assert should_review("selected_variants", hits, "LOW", None, [], False, False, "uid", ["uid"]) is True
