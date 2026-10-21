"""Generate a small but REAL FASTQ pair + reference for the pipeline smoke test.

Produces a genuine paired-end FASTQ (R1/R2) that bwa-mem2 can align and GATK
HaplotypeCaller can call, plus the matching reference FASTA (+ .fai/.dict).
A heterozygous SNV is implanted so the caller has something real to find.
"""
from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

BASES = "ACGT"
SEED = 424242


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="smoke_data")
    ap.add_argument("--length", type=int, default=120000)
    ap.add_argument("--depth", type=int, default=30)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--insert", type=int, default=300)
    ap.add_argument("--variant-pos", type=int, default=60000)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    ref = [rng.choice(BASES) for _ in range(args.length)]
    vpos = args.variant_pos
    ref_base = ref[vpos - 1]
    alt_base = next(b for b in BASES if b != ref_base)
    # The REFERENCE keeps the REF allele; the READS carry the ALT allele, so
    # the aligner/caller sees a real difference to detect.
    ref_seq = "".join(ref)
    alt_seq = ref_seq[:vpos - 1] + alt_base + ref_seq[vpos:]
    ref_path = out / "smoke_reference.fa"
    with open(ref_path, "w") as fh:
        fh.write(">chrSmoke1 smoke test contig\n")
        for i in range(0, len(ref_seq), 60):
            fh.write(ref_seq[i:i + 60] + "\n")

    # FASTA index (samtools faidx compatible)
    with open(str(ref_path) + ".fai", "w") as fh:
        offset = len(">chrSmoke1 smoke test contig\n")
        fh.write(f"chrSmoke1\t{len(ref_seq)}\t{offset}\t60\t61\n")

    # GATK sequence dictionary (GATK looks for <basename>.dict)
    with open(str(ref_path) + ".dict", "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        fh.write("@SQ\tSN:chrSmoke1\tLN:%d\tUR:file:%s\tM5:0\t"
                 "M6:chrSmoke1\tAS:smoke\n" % (len(ref_seq), ref_path))
        fh.write("chrSmoke1\t%d\t0\t60\t61\n" % len(ref_seq))
    dict_gatk = out / "smoke_reference.dict"
    dict_gatk.write_text((str(ref_path) + ".dict").replace(".fa.dict", ".dict")
                          and open(str(ref_path) + ".dict").read(), encoding="utf-8")

    rl = args.read_len
    frag_len = args.insert
    n_frags = max(1, args.length * args.depth // (2 * rl))
    r1_path = out / "SMOKE_R1.fastq.gz"
    r2_path = out / "SMOKE_R2.fastq.gz"
    n_alt = 0
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i in range(n_frags):
            start = rng.randint(1, max(2, len(ref_seq) - frag_len - 1))
            frag = alt_seq[start - 1: start - 1 + frag_len]
            if len(frag) < frag_len:
                continue
            r1 = frag[:rl]
            r2 = revcomp(frag[-rl:])
            # quality: high everywhere, lower in a window around the variant
            def quals(s, offset):
                return "".join(chr(33 + (38 if not (vpos - 1 <= offset < vpos + 30) else 33))
                               for _ in s)
            q1 = quals(r1, start)
            q2 = quals(r2, start + frag_len - rl)
            name = f"read{i:06d}"
            f1.write(f"@{name}\n{r1}\n+\n{q1}\n")
            f2.write(f"@{name}\n{r2}\n+\n{q2}\n")
            if start - 1 <= vpos - 1 < start - 1 + frag_len:
                n_alt += 1

    print(f"reference : {ref_path} ({len(ref_seq)} bp)")
    print(f"fastq     : {r1_path} / {r2_path} ({n_frags} pairs, {n_alt} spanning the variant)")
    print(f"variant   : chrSmoke1:{vpos} {ref_base}>{alt_base} (het, 50% of spanning reads)")
    print(f"reads carry ALT={alt_base} at {vpos}; reference keeps REF={ref_base}")


if __name__ == "__main__":
    main()
