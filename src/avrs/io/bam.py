"""Pure-Python BAM reader/writer + BAI index writer/reader.

BAM structure (little-endian): magic "BAM\\x01", header text, n_ref + ref
dict, then records: [block_size][refID][pos][l_read_name u8][mapq u8][bin
u16][n_cigar_op u16][flag u16][l_seq i32][next_refID][next_pos][tlen i32]
read_name(\\0-term) cigar(n*uint32) seq(nibbles) qual(l_seq) tags.

BAI: bin index (reg2bin) + linear index (16kb windows, min virtual offset).

Used by the pure-python evidence backend (Windows/demo). On Linux/WSL the
pysam backend is preferred; this implementation is validated by round-trip
unit tests. Cross-validation against samtools-indexed files is documented in
VALIDATION.md.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

from .bgzf import BGZFReader, BGZFWriter

SEQ_NT16 = "=ACMGRSVTWYHKDBN"  # index by 4-bit code
_NT16_CODE = {c: i for i, c in enumerate(SEQ_NT16)}
CIGAR_OPS = "MIDNSHP=X"

BAM_FPAIRED = 1
BAM_FPROPER_PAIR = 2
BAM_FUNMAP = 4
BAM_FMUNMAP = 8
BAM_FREVERSE = 16
BAM_FMREVERSE = 32
BAM_FREAD1 = 64
BAM_FREAD2 = 128
BAM_FSECONDARY = 256
BAM_FQCFAIL = 512
BAM_FDUP = 1024
BAM_FSUPPLEMENTARY = 2048


@dataclass
class BamHeader:
    text: str = ""
    refs: list[tuple[str, int]] = field(default_factory=list)

    @property
    def ref_names(self) -> list[str]:
        return [r[0] for r in self.refs]


@dataclass
class BamRecord:
    qname: str
    flag: int
    ref_id: int
    pos: int               # 0-based leftmost aligned position (-1 if unmapped)
    mapq: int
    cigar: list[tuple[int, int]]  # [(op_len, op_char)]
    seq: str
    qual: list[int]        # Phred, per base; [] means '*'
    next_ref_id: int = -1
    next_pos: int = -1
    tlen: int = 0
    tags: bytes = b""
    virtual_offset: int = -1  # filled during sequential iteration

    @property
    def chrom(self) -> Optional[str]:
        return None if self.ref_id < 0 else None  # resolved via header at call sites

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & BAM_FREVERSE)

    @property
    def is_dup(self) -> bool:
        return bool(self.flag & BAM_FDUP)

    @property
    def is_secondary(self) -> bool:
        return bool(self.flag & BAM_FSECONDARY)

    @property
    def is_supplementary(self) -> bool:
        return bool(self.flag & BAM_FSUPPLEMENTARY)

    @property
    def is_qcfail(self) -> bool:
        return bool(self.flag & BAM_FQCFAIL)

    @property
    def is_unmapped(self) -> bool:
        return bool(self.flag & BAM_FUNMAP)

    def reference_length(self) -> int:
        return sum(l for l, op in self.cigar if op in "MDN=X")

    def end_pos(self) -> int:
        """0-based exclusive end of aligned span."""
        return self.pos + self.reference_length()

    def aligned_pairs(self) -> list[tuple[int, int]]:
        """Walk CIGAR: list of (ref_pos0, query_pos0) for M/=/X bases."""
        out = []
        ref = self.pos
        qry = 0
        for op_len, op in self.cigar:
            if op in "M=X":
                for k in range(op_len):
                    out.append((ref + k, qry + k))
                ref += op_len
                qry += op_len
            elif op == "I" or op == "S":
                qry += op_len
            elif op == "D" or op == "N":
                ref += op_len
        return out

    def base_at(self, ref_pos0: int) -> tuple[Optional[str], Optional[int]]:
        """Base and quality aligned at a reference position (0-based), if covered."""
        for r, q in self.aligned_pairs():
            if r == ref_pos0:
                idx = q
                base = self.seq[idx] if idx < len(self.seq) else None
                qual = self.qual[idx] if idx < len(self.qual) else None
                return base, qual
        return None, None


# ---------------------------------------------------------------- read BAM
def read_header(reader: BGZFReader) -> BamHeader:
    magic = reader.read(4)
    if magic != b"BAM\x01":
        raise ValueError("Not a BAM file (bad magic)")
    l_text = struct.unpack("<i", reader.read(4))[0]
    text = reader.read(l_text).decode("utf-8", errors="replace")
    n_ref = struct.unpack("<i", reader.read(4))[0]
    refs = []
    for _ in range(n_ref):
        l_name = struct.unpack("<i", reader.read(4))[0]
        name = reader.read(l_name).rstrip(b"\x00").decode("utf-8", errors="replace")
        l_ref = struct.unpack("<i", reader.read(4))[0]
        refs.append((name, l_ref))
    return BamHeader(text=text, refs=refs)


def _unpack_record(buf: bytes, virtual_offset: int = -1) -> BamRecord:
    (ref_id, pos, l_read_name, mapq, _bin, n_cigar, flag,
     l_seq, next_ref_id, next_pos, tlen) = struct.unpack("<iiBBHHHiiii", buf[:32])
    o = 32
    qname = buf[o:o + l_read_name].rstrip(b"\x00").decode("utf-8", errors="replace")
    o += l_read_name
    cigar_u32 = struct.unpack(f"<{n_cigar}I", buf[o:o + 4 * n_cigar]) if n_cigar else ()
    o += 4 * n_cigar
    n_seq_bytes = (l_seq + 1) // 2
    seq_bytes = buf[o:o + n_seq_bytes]
    o += n_seq_bytes
    qual_bytes = buf[o:o + l_seq]
    o += l_seq
    tags = buf[o:]

    seq = []
    for b in seq_bytes:
        seq.append(SEQ_NT16[b >> 4])
        seq.append(SEQ_NT16[b & 0xF])
    seq = "".join(seq[:l_seq])
    if qual_bytes and qual_bytes[0] == 0xFF and len(set(qual_bytes)) == 1:
        quals: list[int] = []          # '*' - quality not stored
    else:
        # BAM stores raw Phred values (0-93), NOT ASCII-offset scores
        quals = list(qual_bytes)
    cigar = [(c >> 4, CIGAR_OPS[c & 0xF]) for c in cigar_u32]
    return BamRecord(
        qname=qname, flag=flag, ref_id=ref_id, pos=pos, mapq=mapq, cigar=cigar,
        seq=seq, qual=quals, next_ref_id=next_ref_id, next_pos=next_pos,
        tlen=tlen, tags=tags, virtual_offset=virtual_offset,
    )


def iter_records(path: str | Path, header: BamHeader | None = None) -> Iterator[tuple[BamHeader, BamRecord]]:
    """Sequential iteration; yields (header, record) with virtual offsets set."""
    with BGZFReader(path) as reader:
        hdr = read_header(reader)
        if header is not None and header.refs != hdr.refs:
            raise ValueError("BAM header mismatch")
        while True:
            voff = reader.tell_virtual()
            size_bytes = reader.read(4)
            if len(size_bytes) < 4:
                return
            (block_size,) = struct.unpack("<i", size_bytes)
            buf = reader.read(block_size)
            if len(buf) < block_size:
                return
            yield hdr, _unpack_record(buf, voff)


def fetch_region(path: str | Path, chrom: str, start0: int, end0: int,
                 require_index: bool = False) -> list[BamRecord]:
    """Return records overlapping [start0, end0) using the BAI linear index if
    available, else sequential scan (fine for demo-scale BAMs)."""
    # try indexed access
    bai_path = Path(str(path) + ".bai")
    if require_index and not bai_path.exists():
        raise FileNotFoundError(f"BAM index missing: {bai_path}")
    if bai_path.exists():
        try:
            return _fetch_indexed(path, bai_path, chrom, start0, end0)
        except Exception:
            pass  # fall back to scan; index parse failures must not silently corrupt results
    out = []
    for hdr, rec in iter_records(path):
        if rec.ref_id < 0:
            continue
        name = hdr.refs[rec.ref_id][0]
        if name == chrom and rec.pos < end0 and rec.end_pos() > start0:
            out.append(rec)
    return out


# ---------------------------------------------------------------- BAI
def reg2bin(beg: int, end: int) -> int:
    end -= 1
    if beg >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (beg >> 26)
    if beg >> 29 == end >> 29:
        return ((1 << 0) - 1) // 7 + (beg >> 29)
    return 0


def _parse_bai(bai_path: Path):
    """Return (ref_names_placeholder=None, n_ref, bins, linear) structures.

    bins: list over refs of dict bin -> list[(vbeg, vend)]
    linear: list over refs of list of ioffsets (per 16kb window)
    """
    data = open(bai_path, "rb").read()
    if data[:4] != b"BAI\x01":
        raise ValueError("Not a BAI index")
    n_ref = struct.unpack_from("<i", data, 4)[0]
    off = 8
    bins = []
    linear = []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", data, off)[0]
        off += 4
        bd = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", data, off)
            off += 8
            chunks = []
            for _ in range(n_chunk):
                vb, ve = struct.unpack_from("<QQ", data, off)
                off += 16
                chunks.append((vb, ve))
            bd[bin_id] = chunks
        n_intv = struct.unpack_from("<i", data, off)[0]
        off += 4
        ints = list(struct.unpack_from(f"<{n_intv}Q", data, off)) if n_intv else []
        off += 8 * n_intv
        bins.append(bd)
        linear.append(ints)
    return n_ref, bins, linear


def _fetch_indexed(bam_path: Path, bai_path: Path, chrom: str, start0: int, end0: int) -> list[BamRecord]:
    hdr = None
    with BGZFReader(bam_path) as reader:
        hdr = read_header(reader)
    ref_id = next((i for i, (n, _) in enumerate(hdr.refs) if n == chrom), -1)
    if ref_id < 0:
        return []
    _, bins, linear = _parse_bai(bai_path)
    if ref_id >= len(bins):
        return []
    # linear index: min virtual offset for each 16kb window overlapping query
    win_beg = start0 >> 14
    win_end = end0 >> 14
    if linear[ref_id]:
        start_voff = linear[ref_id][win_beg] if win_beg < len(linear[ref_id]) else linear[ref_id][-1]
        if start_voff == 0:
            start_voff = linear[ref_id][max(0, win_beg - 1)] if linear[ref_id] else 0
    else:
        start_voff = 0
    # candidate bins
    wanted = {reg2bin(s, min(s + 16384, end0)) for s in range(start0, end0, 16384)}
    wanted.add(0)  # include reads spanning windows; chunk ranges also filter below
    offsets: list[tuple[int, int]] = []
    for bin_id in wanted:
        for (vb, ve) in bins[ref_id].get(bin_id, []):
            offsets.append((vb, ve))
    out = []
    with BGZFReader(bam_path) as reader:
        seen = set()
        for vb, ve in sorted(offsets):
            begin = max(vb, start_voff)
            reader.seek_virtual(begin)
            limit = ve
            while reader.tell_virtual() < limit:
                voff = reader.tell_virtual()
                size_bytes = reader.read(4)
                if len(size_bytes) < 4:
                    break
                (block_size,) = struct.unpack("<i", size_bytes)
                buf = reader.read(block_size)
                if len(buf) < block_size:
                    break
                if voff in seen:
                    continue
                seen.add(voff)
                rec = _unpack_record(buf, voff)
                if rec.ref_id == ref_id and rec.pos < end0 and rec.end_pos() > start0:
                    out.append(rec)
    return out


# ---------------------------------------------------------------- write BAM
def encode_record(rec: BamRecord, ref_id: int) -> bytes:
    l_read_name = len(rec.qname) + 1
    n_cigar = len(rec.cigar)
    cigar_u32 = [(l << 4) | CIGAR_OPS.index(op) for l, op in rec.cigar]
    seq = rec.seq.upper()
    l_seq = len(seq)
    nib = bytearray((l_seq + 1) // 2)
    for i, c in enumerate(seq):
        code = _NT16_CODE.get(c, 15)
        if i % 2 == 0:
            nib[i // 2] |= code << 4
        else:
            nib[i // 2] |= code
    if rec.qual:
        qual_bytes = bytes(min(93, max(0, q)) for q in rec.qual)
    else:
        qual_bytes = b"\xff" * l_seq
    head = struct.pack(
        "<iiBBHHHiiii", ref_id, rec.pos, l_read_name, rec.mapq, reg2bin(rec.pos, rec.pos + rec.reference_length()),
        n_cigar, rec.flag, l_seq, rec.next_ref_id, rec.next_pos, rec.tlen,
    )
    body = (
        head
        + rec.qname.encode() + b"\x00"
        + (struct.pack(f"<{n_cigar}I", *cigar_u32) if n_cigar else b"")
        + bytes(nib) + qual_bytes + rec.tags
    )
    return struct.pack("<i", len(body)) + body


def write_bam(path: str | Path, header: BamHeader, records: list[BamRecord]) -> None:
    """Write a coordinate-sorted BAM with BAI index (records must be sorted)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        w = BGZFWriter(fh)
        text = header.text if header.text else "@HD\tVN:1.6\tSO:coordinate\n"
        ref_blob = struct.pack("<i", len(header.refs))
        for name, l_ref in header.refs:
            nb = name.encode() + b"\x00"
            ref_blob += struct.pack("<i", len(nb)) + nb + struct.pack("<i", l_ref)
        w.write(b"BAM\x01" + struct.pack("<i", len(text)) + text.encode() + ref_blob)
        for i, rec in enumerate(records):
            rid = rec.ref_id if rec.ref_id >= 0 else -1
            w.write(encode_record(rec, rid))
        w.close()
    write_bai(path, header, records)


def write_bai(bam_path: Path, header: BamHeader, records: list[BamRecord]) -> None:
    """Generate a BAI index from the BAM we just wrote.

    BAI layout: per reference -> bin index (reg2bin) + linear index (16 kb
    windows). A bin chunk spans a *range of records*, so chunk end must be
    the virtual offset of the record following the last one in the run -
    not the last record's own offset.
    """
    # single sequential pass: (ref_id, pos, end, voff, next_voff)
    entries: list[tuple[int, int, int, int, int]] = []
    with BGZFReader(bam_path) as reader:
        read_header(reader)
        while True:
            voff = reader.tell_virtual()
            size_bytes = reader.read(4)
            if len(size_bytes) < 4:
                break
            (block_size,) = struct.unpack("<i", size_bytes)
            buf = reader.read(block_size)
            if len(buf) < block_size:
                break
            rec = _unpack_record(buf, voff)
            next_voff = voff + 4 + block_size
            if rec.ref_id >= 0:
                entries.append((rec.ref_id, rec.pos, rec.end_pos(), voff, next_voff))

    per_ref: list[list[tuple[int, int, int, int, int]]] = [[] for _ in header.refs]
    for e in entries:
        if e[0] < len(per_ref):
            per_ref[e[0]].append(e)

    bins_out: list[dict[int, list[tuple[int, int]]]] = []
    linear_out: list[list[int]] = []
    for ref_entries in per_ref:
        bd: dict[int, list[tuple[int, int]]] = {}
        li: dict[int, int] = {}
        # group consecutive records sharing a bin into chunks
        run_bin: Optional[int] = None
        run_start = 0
        run_end = 0
        for (rid, pos, end, voff, next_voff) in ref_entries:
            bin_id = reg2bin(pos, end)
            for wstart in range(pos >> 14, ((end - 1) >> 14) + 1):
                if wstart not in li or voff < li[wstart]:
                    li[wstart] = voff
            if run_bin is None:
                run_bin, run_start, run_end = bin_id, voff, next_voff
            elif bin_id == run_bin:
                run_end = next_voff
            else:
                bd.setdefault(run_bin, []).append((run_start, run_end))
                run_bin, run_start, run_end = bin_id, voff, next_voff
        if run_bin is not None:
            bd.setdefault(run_bin, []).append((run_start, run_end))
        bins_out.append(bd)
        n_windows = (max(li) + 1) if li else 0
        linear = [li.get(i, 0) for i in range(n_windows)]
        last = 0
        for i in range(len(linear)):
            if linear[i] == 0:
                linear[i] = last
            last = linear[i]
        linear_out.append(linear)

    buf = bytearray(b"BAI\x01" + struct.pack("<i", len(header.refs)))
    for bd, linear in zip(bins_out, linear_out):
        buf += struct.pack("<i", len(bd))
        for bin_id in sorted(bd):
            chunks = bd[bin_id]
            buf += struct.pack("<Ii", bin_id, len(chunks))
            for vb, ve in chunks:
                buf += struct.pack("<QQ", vb, ve)
        buf += struct.pack("<i", len(linear))
        for v in linear:
            buf += struct.pack("<Q", v)
    with open(str(bam_path) + ".bai", "wb") as fh:
        fh.write(bytes(buf))
