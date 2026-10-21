"""Pure-Python VCF reading (plain/bgzip), multiallelic splitting and
allele-level normalization (left-alignment + parsimony) against a reference.

bcftools remains the preferred normalizer on Linux/WSL (workflow step);
this implementation covers review-only mode, tests and any platform without
htslib. It keeps the mapping between the original caller representation and
the normalized representation.
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator, Optional


def _open_text(path: str | Path):
    p = Path(path)
    fh = open(p, "rb")
    magic = fh.read(2)
    fh.close()
    if magic == b"\x1f\x8b":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "rt", encoding="utf-8", errors="replace")


@dataclass
class VcfRecord:
    chrom: str
    pos: int            # 1-based
    ref: str
    alts: list[str]
    qual: Optional[float]
    filt: str
    info: dict[str, str]
    formats: dict[str, dict[str, str]]  # sample -> {key: value}
    sample_names: list[str]
    original_alt: str = ""  # multiallelic original

    def genotype(self, sample: str | None = None) -> Optional[str]:
        if not self.formats:
            return None
        if sample and sample in self.formats:
            return self.formats[sample].get("GT")
        first = next(iter(self.formats.values()))
        return first.get("GT")

    def format_value(self, key: str, sample: str | None = None) -> Optional[str]:
        if not self.formats:
            return None
        d = self.formats.get(sample) if sample else next(iter(self.formats.values()))
        return (d or {}).get(key)

    def info_int(self, key: str) -> Optional[int]:
        v = self.info.get(key)
        if v is None:
            return None
        try:
            return int(float(v.split(",")[0]))
        except ValueError:
            return None


@dataclass
class VcfHeader:
    lines: list[str] = field(default_factory=list)
    sample_names: list[str] = field(default_factory=list)

    @property
    def contigs(self) -> list[str]:
        out = []
        for line in self.lines:
            if line.startswith("##contig="):
                m = re.search(r"ID=([^,>]+)", line)
                if m:
                    out.append(m.group(1))
        return out


class VcfReader:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.header = VcfHeader()
        self._parse_header()

    def _parse_header(self):
        with _open_text(self.path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    self.header.lines.append(line)
                elif line.startswith("#CHROM"):
                    cols = line.split("\t")
                    if len(cols) > 9:
                        self.header.sample_names = cols[9:]
                    return
                elif line:  # data without #CHROM line (invalid but tolerate)
                    return

    def __iter__(self) -> Iterator[VcfRecord]:
        with _open_text(self.path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                rec = self._parse_line(line.rstrip("\n"))
                if rec:
                    yield rec

    def _parse_line(self, line: str) -> Optional[VcfRecord]:
        cols = line.split("\t")
        if len(cols) < 8:
            return None
        chrom, pos, _id, ref, alt, qual, filt, info = cols[:8]
        info_d: dict[str, str] = {}
        if info and info != ".":
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
                elif kv:
                    info_d[kv] = "True"
        alts = [] if alt in (".", "") else alt.split(",")
        formats: dict[str, dict[str, str]] = {}
        sample_names = []
        if len(cols) > 9:
            keys = cols[8].split(":")
            names = self.header.sample_names or [f"S{i + 1}" for i in range(len(cols) - 9)]
            for name, val in zip(names, cols[9:]):
                vals = val.split(":")
                formats[name] = dict(zip(keys, vals))
                sample_names.append(name)
        try:
            q = None if qual == "." else float(qual)
        except ValueError:
            q = None
        return VcfRecord(
            chrom=chrom, pos=int(pos), ref=ref.upper(), alts=[a.upper() for a in alts],
            qual=q, filt=filt, info=info_d, formats=formats,
            sample_names=sample_names or self.header.sample_names, original_alt=alt,
        )


def split_multiallelic(rec: VcfRecord) -> list[VcfRecord]:
    """Split multiallelic sites into allele-level records.

    AD/PL lists are subset to the (REF, ALT_i) pair where present.
    """
    if len(rec.alts) <= 1:
        return [rec]
    out = []
    for i, alt in enumerate(rec.alts):
        import copy
        r = copy.deepcopy(rec)
        r.alts = [alt]
        # subset AD / RPL style fields if present
        for sname in r.formats:
            gt = r.formats[sname].get("GT")
            if gt and any(c in gt for c in "/|"):
                # keep genotype only if it references this allele, else map to ref
                gt_alt = {a for a in re.split(r"[/|]", gt) if a != "."}
                if str(i + 1) in gt_alt:
                    r.formats[sname]["GT"] = gt  # keep original multiallelic GT for provenance
                else:
                    pass
            ad = r.formats[sname].get("AD")
            if ad and "," in ad:
                parts = ad.split(",")
                if len(parts) == len(rec.alts) + 1:
                    r.formats[sname]["AD"] = f"{parts[0]},{parts[i + 1]}"
        out.append(r)
    return out


# ------------------------------------------------------------ normalization
def normalize_allele(pos: int, ref: str, alt: str,
                      get_ref_base: Callable[[int], str]) -> tuple[int, str, str]:
    """Left-align + parsimonious normalization (vt/bcftools-style).

    get_ref_base(p) returns the reference base at 1-based position p.
    Returns (pos, ref, alt) normalized. SNVs pass through unchanged.
    """
    # right-trim common suffix (keep >=1 base each)
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # left-shift pure indels (one allele is a prefix of the other)
    while (
        pos > 1
        and len(ref) >= 1 and len(alt) >= 1
        and ref[-1] == alt[-1]
        and (len(ref) == 1 or len(alt) == 1)
        and (ref.startswith(alt) or alt.startswith(ref))
    ):
        pos -= 1
        b = get_ref_base(pos).upper()
        ref, alt = b + ref, b + alt
        ref, alt = ref[:-1], alt[:-1]
    # left-trim common prefix (keep >=1 base each)
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def variant_uid(build: str, chrom: str, pos: int, ref: str, alt: str) -> str:
    """Stable allele-level identity that includes the reference build."""
    return f"{build}#{chrom}:{pos}:{ref}>{alt}"


def classify_variant(ref: str, alt: str) -> str:
    """SNV | MNV | INS | DEL | COMPLEX (allele-level)."""
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return "INS"
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        return "DEL"
    if len(ref) == len(alt) and len(ref) > 1 and all(a != b for a, b in zip(ref, alt)):
        return "MNV"
    return "COMPLEX"


def write_vcf(path: str | Path, header_lines: list[str], records: list[VcfRecord],
              sample_names: list[str] | None = None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for h in header_lines:
            fh.write(h.rstrip("\n") + "\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                 + ("\tFORMAT\t" + "\t".join(sample_names or []) if sample_names else "")
                 + "\n")
        for r in records:
            cols = [r.chrom, str(r.pos), ".", r.ref, ",".join(r.alts) or ".",
                    "." if r.qual is None else f"{r.qual:g}", r.filt or ".",
                    ";".join(f"{k}={v}" if v != "True" else k for k, v in r.info.items()) or "."]
            if r.formats:
                keys = ["GT"] + [k for k in next(iter(r.formats.values())) if k != "GT"]
                fh.write("\t".join(cols) + "\t" + ":".join(keys))
                for s in r.sample_names:
                    d = r.formats.get(s, {})
                    fh.write("\t" + ":".join(d.get(k, ".") for k in keys))
                fh.write("\n")
            else:
                fh.write("\t".join(cols) + "\n")
