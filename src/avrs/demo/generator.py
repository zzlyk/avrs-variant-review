"""Reproducible demo dataset generator (spec 11.3).

Builds a small synthetic genome (no patient data) with 7 review scenarios:
  1 clean ALT-supported SNV
  2 low-quality / low-depth SNV
  3 strand-biased SNV
  4 indel inside a homopolymer
  5 insufficient coverage target locus
  6 variant present in BAM but absent from the VCF (caller miss)
  7 equivocal variant (cannot be resolved)

Simulated labels are scenario labels for pipeline testing - they are NOT
clinical truth and must never be used as a truth set.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

SEED = 20261020
BASES = "ACGT"


@dataclass
class Scenario:
    name: str
    chrom: str
    pos: int          # 1-based
    ref: str
    alt: str
    in_vcf: bool
    description: str
    n_reads_ref: int = 12
    n_reads_alt: int = 12
    alt_baseq: int = 37
    ref_baseq: int = 37
    alt_strand: str = "both"     # both | reverse | forward
    same_start: bool = False
    dup_rate: float = 0.0
    mapq: int = 60
    insert_seq: str = ""
    read_len: int = 100
    genotype: str = "0/1"
    vcf_filter: str = "PASS"
    vcf_gq: int = 99
    hotspot: bool = False
    # override the caller-reported AD to simulate a caller/pileup discrepancy
    vcf_ad_override: Optional[tuple[int, int]] = None
    # left-aligned representation (computed from the reference); the VCF keeps
    # the original representation, the reads are built at the aligned one
    norm_pos: int = 0
    norm_ref: str = ""
    norm_alt: str = ""


def make_reference(lengths: dict[str, int], seed: int = SEED) -> dict[str, str]:
    rng = random.Random(seed)
    contigs = {}
    for name, length in lengths.items():
        seq = [rng.choice(BASES) for _ in range(length)]
        contigs[name] = "".join(seq)
    return contigs


def plant_homopolymer(contigs: dict[str, str], chrom: str, pos: int, base: str, length: int):
    seq = list(contigs[chrom])
    seq[pos - 1: pos - 1 + length] = list(base * length)
    contigs[chrom] = "".join(seq)


def plant_repeat(contigs: dict[str, str], chrom: str, pos: int, unit: str, copies: int):
    seq = list(contigs[chrom])
    seq[pos - 1: pos - 1 + len(unit) * copies] = list(unit * copies)
    contigs[chrom] = "".join(seq)


def build_scenarios() -> list[Scenario]:
    plant_specs = []
    sc = [
        Scenario("clean_snv", "chrDemo1", 5000, "A", "G", True,
                 "Strong balanced ALT support, high base quality -> likely true"),
        Scenario("lowqual_snv", "chrDemo1", 8000, "C", "T", True,
                 "Few ALT reads, low base quality, identical start positions -> likely artifact",
                 n_reads_ref=8, n_reads_alt=4, alt_baseq=6, ref_baseq=37,
                 same_start=True, mapq=12, vcf_filter="LowQual"),
        Scenario("strand_bias_snv", "chrDemo1", 11000, "G", "C", True,
                 "ALT reads all on reverse strand -> strand bias suspicion",
                 n_reads_ref=10, n_reads_alt=14, alt_strand="reverse"),
        Scenario("homopolymer_indel", "chrDemo1", 14000, "A", "AAA", True,
                 "Insertion inside poly-A run -> indel representation risk",
                 n_reads_ref=10, n_reads_alt=8, genotype="0/1"),
        Scenario("low_coverage", "chrDemo1", 17000, "T", "C", True,
                 "Almost no coverage -> insufficient evidence",
                 n_reads_ref=1, n_reads_alt=1, alt_baseq=10, vcf_filter="LowDP"),
        Scenario("caller_miss", "chrDemo1", 20000, "G", "T", False,
                 "Reads support ALT but variant is absent from the VCF -> FN candidate",
                 n_reads_ref=10, n_reads_alt=9, hotspot=True),
        Scenario("equivocal", "chrDemo1", 23000, "C", "A", True,
                 "Moderate quality + caller/pileup disagreement -> inconclusive",
                 n_reads_ref=10, n_reads_alt=8, alt_baseq=20, alt_strand="both",
                 same_start=False, genotype="0/1", vcf_filter="VQSRUninformative",
                 vcf_ad_override=(16, 2)),
        Scenario("clean_snv_chr2", "chrDemo2", 5000, "T", "A", True,
                 "Clean SNV on second contig -> contig handling check",
                 n_reads_ref=10, n_reads_alt=11),
    ]
    return sc


def _apply_alt_window(window: str, ref: str, alt: str, offset: int) -> tuple[str, list]:
    """Replace the reference allele inside a read window with the ALT allele.

    offset = index of the variant anchor inside the window (0-based).
    Returns (new_sequence, cigar) where the cigar reflects the inserted/
    deleted reference bases.
    """
    o = offset
    ref_len = len(ref)
    anchor = window[o:o + ref_len]
    if anchor.upper() != ref.upper():
        # reference does not match - fall back to a plain M alignment
        return window, [(len(window), "M")]
    inserted = alt[1:] if len(alt) > 1 else ""
    deleted = ref[1:] if len(ref) > 1 else ""
    if inserted and not deleted:                    # pure insertion
        new_seq = window[:o + 1] + inserted + window[o + 1:]
        cig = [(o + 1, "M"), (len(inserted), "I"),
               (len(window) - o - 1, "M")] if len(window) - o - 1 > 0 else [(o + 1, "M"), (len(inserted), "I")]
    elif deleted and not inserted:                  # pure deletion
        new_seq = window[:o + 1] + window[o + ref_len:]
        cig = [(o + 1, "M"), (len(deleted), "D"),
               (len(window) - o - ref_len, "M")] if len(window) - o - ref_len > 0 else [(o + 1, "M"), (len(deleted), "D")]
    elif inserted and deleted:                      # replacement (del + ins)
        new_seq = window[:o + 1] + inserted + window[o + ref_len:]
        cig = [(o + 1, "M"), (len(deleted), "D"), (len(inserted), "I"),
               (len(window) - o - ref_len, "M")]
    else:                                            # SNV / MNV
        new_seq = window[:o] + alt + window[o + ref_len:]
        cig = [(len(new_seq), "M")]
    return new_seq, cig


def simulate_reads(contigs: dict[str, str], sc: Scenario, rng: random.Random,
                   read_len: int = 100, insert_size: int = 300):
    """Generate (seq, qual, flag, mapq, start, cigar) tuples supporting a scenario.

    Reads are simulated at the LEFT-ALIGNED representation of the variant
    (as a normalised VCF would report it) so that the evidence engine sees
    the same allele representation that bcftools norm produces.
    """
    ref_seq = contigs[sc.chrom]
    reads = []
    total = sc.n_reads_ref + sc.n_reads_alt
    v_ref, v_alt = sc.norm_ref, sc.norm_alt
    v_pos = sc.norm_pos
    for i in range(total):
        is_alt = i >= sc.n_reads_ref
        if sc.same_start:
            start = max(0, v_pos - 1 - read_len // 2)
        else:
            start = max(0, v_pos - 1 - rng.randint(10, 60))
        end = start + read_len
        if end > len(ref_seq):
            start = max(0, len(ref_seq) - read_len - 1)
            end = start + read_len
        window = list(ref_seq[start:end])
        if is_alt:
            offset = v_pos - 1 - start
            if 0 <= offset < len(window):
                seq, cig = _apply_alt_window("".join(window), v_ref, v_alt, offset)
            else:
                seq = "".join(window)
                cig = [(len(seq), "M")]
        else:
            seq = "".join(window)
            cig = [(len(seq), "M")]
        baseq = sc.alt_baseq if is_alt else sc.ref_baseq
        qual = [baseq] * len(seq)
        # realistic orientation: ~50% reverse, mate on the opposite side
        reverse = (sc.alt_strand == "reverse" and is_alt) or \
                  (sc.alt_strand == "forward" and is_alt) or \
                  (sc.alt_strand == "both" and rng.random() < 0.5)
        flag = 0x1 | 0x2 | 0x20           # paired, proper pair, first in pair
        if reverse:
            flag |= 0x10                  # reverse strand
        else:
            flag |= 0x40                  # second in pair
        reads.append({"seq": seq, "qual": qual, "flag": flag, "mapq": sc.mapq,
                      "start": start, "cigar": cig, "alt": is_alt,
                      "name": f"{sc.name}_read{i:04d}"})
    return reads


def write_bam_with_pysam(path: Path, contigs: dict[str, str], reads_by_scenario,
                          sample_id: str) -> bool:
    try:
        import pysam
    except ImportError:
        return False
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(seq)} for name, seq in contigs.items()],
        "RG": [{"ID": "rg1", "SM": sample_id, "PL": "ILLUMINA", "LB": "lib1"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        all_reads = []
        for sc, reads in reads_by_scenario:
            for r in reads:
                a = pysam.AlignedSegment()
                a.query_name = r["name"]
                a.query_sequence = r["seq"]
                a.query_qualities = pysam.qualitystring_to_array("".join(
                    chr(q + 33) for q in r["qual"]))
                a.flag = r["flag"]
                a.reference_id = list(contigs).index(sc.chrom)
                a.reference_start = r["start"]
                a.mapping_quality = r["mapq"]
                a.cigarstring = "".join(f"{l}{op}" for l, op in r["cigar"])
                a.next_reference_id = a.reference_id
                a.next_reference_start = r["start"] + 250
                a.template_length = 300
                a.set_tag("RG", "rg1")
                all_reads.append((sc.chrom, r["start"], a))
        all_reads.sort(key=lambda x: (list(contigs).index(x[0]), x[1]))
        for _, _, a in all_reads:
            out.write(a)
    pysam.index(str(path))
    return True


def write_bam_purepython(path: Path, contigs: dict[str, str], reads_by_scenario,
                         sample_id: str) -> None:
    from ..io.bam import BamHeader, BamRecord, write_bam
    from ..io.fasta import write_fasta, write_fai
    header = BamHeader(
        text="@HD\tVN:1.6\tSO:coordinate\n@RG\tID:rg1\tSM:%s\tPL:ILLUMINA\tLB:lib1\n" % sample_id,
        refs=[(name, len(seq)) for name, seq in contigs.items()])
    recs = []
    ref_names = list(contigs)
    for sc, reads in reads_by_scenario:
        for r in reads:
            cig = [(l, op) for l, op in r["cigar"]]
            ref_len = sum(l for l, op in cig if op in "MDN=X")
            recs.append(BamRecord(
                qname=r["name"], flag=r["flag"],
                ref_id=ref_names.index(sc.chrom), pos=r["start"], mapq=r["mapq"],
                cigar=cig, seq=r["seq"], qual=r["qual"],
                next_ref_id=ref_names.index(sc.chrom), next_pos=r["start"] + 250, tlen=300,
            ))
    recs.sort(key=lambda r: (r.ref_id, r.pos))
    write_bam(path, header, recs)


def write_vcf(path: Path, contigs: dict[str, str], scenarios: list[Scenario],
              sample_id: str):
    lines = [
        "##fileformat=VCFv4.2",
        "##source=avrs-demo-generator",
        "##reference=demo-reference",
    ]
    for name, seq in contigs.items():
        lines.append(f"##contig=<ID={name},length={len(seq)}>")
    lines += [
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}",
    ]
    for sc in scenarios:
        if not sc.in_vcf:
            continue
        dp = sc.n_reads_ref + sc.n_reads_alt
        ad_ref, ad_alt = sc.n_reads_ref, sc.n_reads_alt
        if sc.vcf_ad_override:
            ad_ref, ad_alt = sc.vcf_ad_override
        ad = f"{ad_ref},{ad_alt}"
        lines.append(
            f"{sc.chrom}\t{sc.pos}\t.\t{sc.ref}\t{sc.alt}\t{50 + sc.n_reads_alt}\t"
            f"{sc.vcf_filter}\tDP={dp}\tGT:DP:GQ:AD\t{sc.genotype}:{dp}:{sc.vcf_gq}:{ad}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _sync_scenario_alleles(contigs: dict[str, str], scenarios: list[Scenario]) -> None:
    """Make each scenario's REF match the actual reference base(s) and
    compute the left-aligned representation used for read simulation."""
    from ..io.vcf import normalize_allele
    for sc in scenarios:
        ref_seq = contigs[sc.chrom]
        actual = ref_seq[sc.pos - 1]
        if len(sc.ref) == 1:
            sc.ref = actual
            if len(sc.alt) == 1:
                sc.alt = next(b for b in "ACGT" if b != actual)
            else:
                sc.alt = actual + sc.alt[1:]
        else:
            real = ref_seq[sc.pos - 1: sc.pos - 1 + len(sc.ref)]
            sc.ref = real
            if len(sc.alt) == 1:
                sc.alt = real[0]
            else:
                sc.alt = real[0] + sc.alt[1:]
        getter = lambda p, c=sc.chrom: ref_seq[p - 1]  # noqa: E731
        sc.norm_pos, sc.norm_ref, sc.norm_alt = normalize_allele(
            sc.pos, sc.ref, sc.alt, getter)


def generate_demo(out_dir: str | Path, sample_id: str = "DEMO_SAMPLE_001",
                  build: str = "DEMO-GRCh38") -> dict:
    """Generate the full demo dataset. Returns a manifest dict."""
    from ..io.fasta import write_fasta, write_fai
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    contigs = make_reference({"chrDemo1": 30000, "chrDemo2": 15000})
    plant_homopolymer(contigs, "chrDemo1", 14000, "A", 9)
    plant_repeat(contigs, "chrDemo2", 9000, "AT", 12)
    scenarios = build_scenarios()
    _sync_scenario_alleles(contigs, scenarios)

    fasta = out / "demo_reference.fa"
    write_fasta(fasta, contigs)
    write_fai(fasta, contigs)

    reads_by_scenario = []
    for sc in scenarios:
        reads = simulate_reads(contigs, sc, random.Random(SEED + sc.pos))
        reads_by_scenario.append((sc, reads))

    bam = out / "demo_sample.bam"
    used_pysam = write_bam_with_pysam(bam, contigs, reads_by_scenario, sample_id)
    if not used_pysam:
        write_bam_purepython(bam, contigs, reads_by_scenario, sample_id)

    vcf = out / "demo_sample.vcf"
    write_vcf(vcf, contigs, scenarios, sample_id)

    bed = out / "demo_target.bed"
    bed.write_text("chrDemo1\t1000\t25000\nchrDemo2\t500\t10000\n", encoding="utf-8")

    hotspot = out / "demo_hotspots.tsv"
    hs_rows = ["#chrom\tpos\tref\talt"]
    for sc in scenarios:
        if sc.hotspot or not sc.in_vcf:
            hs_rows.append(f"{sc.chrom}\t{sc.pos}\t{sc.ref}\t{sc.alt}")
    hotspot.write_text("\n".join(hs_rows) + "\n", encoding="utf-8")

    manifest = {
        "sample_id": sample_id, "reference_build": build,
        "fasta": str(fasta), "bam": str(bam), "bam_index": str(bam) + ".bai",
        "vcf": str(vcf), "target_bed": str(bed), "hotspots": str(hotspot),
        "bam_written_with": "pysam" if used_pysam else "purepython",
        "scenarios": [
            {"name": s.name, "chrom": s.chrom, "pos": s.pos, "ref": s.ref, "alt": s.alt,
             "in_vcf": s.in_vcf, "description": s.description} for s in scenarios],
        "disclaimer": ("Simulated data with scenario labels for pipeline testing only - "
                       "NOT clinical truth and not a truth set."),
    }
    import json
    (out / "demo_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
