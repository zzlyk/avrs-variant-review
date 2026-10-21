"""Unit tests: VCF parsing, normalization, variant ID stability."""
from pathlib import Path

import pytest

from avrs.io.vcf import (VcfReader, classify_variant, normalize_allele,
                         split_multiallelic, variant_uid)

VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tG\t50\tPASS\tDP=30\tGT:AD:DP\t0/1:15,15:30
chr1\t200\t.\tC\tT,G\t60\tPASS\tDP=20\tGT:AD:DP\t1/2:6,7,7:20
chr1\t300\t.\tT\t.\t.\t.\t.\tGT\t./.
chr1\t400\t.\tGGA\tG\t40\tLowQual\tDP=12\tGT:AD:DP\t0/1:8,4:12
"""


def _write(tmp_path, text=VCF):
    p = tmp_path / "t.vcf"
    p.write_text(text, encoding="utf-8")
    return p


def test_vcf_header_and_records(tmp_path):
    reader = VcfReader(_write(tmp_path))
    assert reader.header.contigs == ["chr1"]
    assert reader.header.sample_names == ["S1"]
    recs = list(reader)
    assert len(recs) == 4
    assert recs[0].chrom == "chr1" and recs[0].pos == 100
    assert recs[0].info["DP"] == "30"
    assert recs[0].genotype() == "0/1"
    assert recs[0].format_value("AD") == "15,15"


def test_multiallelic_split(tmp_path):
    reader = VcfReader(_write(tmp_path))
    recs = list(reader)
    split = split_multiallelic(recs[1])
    assert len(split) == 2
    assert split[0].alts == ["T"] and split[1].alts == ["G"]
    # AD is subset per ALT
    assert split[0].formats["S1"]["AD"] == "6,7"
    assert split[1].formats["S1"]["AD"] == "6,7"


def test_missing_alt_is_skipped_not_crashed(tmp_path):
    reader = VcfReader(_write(tmp_path))
    recs = list(reader)
    assert recs[2].alts == []
    assert recs[2].genotype() == "./."


def test_normalize_snv_unchanged():
    ref = lambda p: "A"  # noqa: E731
    assert normalize_allele(100, "A", "G", ref) == (100, "A", "G")


def test_normalize_indel_left_aligns_in_homopolymer():
    # reference: pos1='T', pos2..7='G' ; VCF at pos2 REF=GGG ALT=G (delete 2 G's)
    seq = "T" + "G" * 6
    def ref(p):
        return seq[p - 1]
    pos, r, a = normalize_allele(2, "GGG", "G", ref)
    # left alignment moves the anchor to the leftmost possible position (pos 1)
    assert pos == 1
    assert r == "TGG" and a == "T"


def test_normalize_no_shift_needed_in_unique_context():
    seq = "ACGTACGTAC"
    pos, r, a = normalize_allele(3, "G", "GTT", lambda p: seq[p - 1])
    assert r == "G" and a.startswith("G")
    assert len(a) > len(r)


def test_variant_uid_includes_build_and_is_stable():
    a = variant_uid("GRCh38", "chr1", 100, "A", "G")
    b = variant_uid("GRCh38", "chr1", 100, "A", "G")
    c = variant_uid("GRCh37", "chr1", 100, "A", "G")
    assert a == b
    assert a != c
    assert "GRCh38" in a and "chr1:100" in a and "A>G" in a


def test_classify_variant():
    assert classify_variant("A", "G") == "SNV"
    assert classify_variant("AG", "A") == "DEL"
    assert classify_variant("A", "AGG") == "INS"
    assert classify_variant("AG", "TT") == "MNV"
    assert classify_variant("AGG", "AT") == "COMPLEX"
