"""Pure-Python BGZF (blocked gzip) reader/writer and BAM/BAI support.

This module is the fallback evidence backend for platforms where pysam
cannot be built (native Windows). On WSL2/Linux the pysam backend is
preferred and automatically selected when importable; these implementations
exist so that tests, demos and the audit path run everywhere.

BGZF: each block is a standard gzip member with an extra field 'BC' holding
the compressed block size. Python's gzip module reads concatenated members
transparently; writing needs explicit block construction.
"""
from __future__ import annotations

import struct
import zlib
from pathlib import Path
from typing import BinaryIO, Iterator

# The standard 28-byte BGZF EOF marker block.
BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def bgzf_block_size(raw: bytes) -> int:
    """Build a single BGZF block from uncompressed bytes. Returns block bytes."""
    assert len(raw) <= 0x10000 - 26, "BGZF block too large"
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    compressed = compressor.compress(raw) + compressor.flush()
    # total block size = 18 (header) + len(compressed) + 8 (crc+isize)
    # BSIZE stores (total - 1)
    bsize = len(compressed) + 25
    block = struct.pack(
        "<BBBBLBBHBBHH",
        31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize,
    )
    block += compressed
    block += struct.pack("<LL", zlib.crc32(raw) & 0xFFFFFFFF, len(raw) & 0xFFFFFFFF)
    return block


class BGZFWriter:
    """Write data as a sequence of BGZF blocks (default 64 KiB uncompressed)."""

    def __init__(self, fh: BinaryIO, block_size: int = 0xFF00):
        self.fh = fh
        self.block_size = block_size
        self._buffer = bytearray()

    def write(self, data: bytes):
        self._buffer += data
        while len(self._buffer) >= self.block_size:
            self._flush_one(self.block_size)

    def _flush_one(self, size: int):
        raw = bytes(self._buffer[:size])
        del self._buffer[:size]
        self.fh.write(bgzf_block_size(raw))

    def close(self):
        if self._buffer:
            self._flush_one(len(self._buffer))
        self.fh.write(BGZF_EOF)
        self.fh.flush()

    # context manager -----------------------------------------------------
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if not exc[0]:
            self.close()
        else:
            self.fh.flush()


def open_bgzf_writer(path: str | Path, **kw) -> BGZFWriter:
    return BGZFWriter(open(path, "wb"), **kw)


class BGZFReader:
    """Sequential reader over BGZF that tracks virtual offsets.

    Virtual offset: (block_offset_in_file << 16) | offset_within_block.
    Needed to write/interpret BAI indices.
    """

    def __init__(self, path: str | Path):
        self.fh = open(path, "rb")
        self._block: bytes = b""
        self._uoff = 0        # offset within current decompressed block
        self._coff = 0        # offset of current block in file
        self._block_end_coff = 0

    def _read_block(self):
        start = self.fh.tell()
        head = self.fh.read(12)          # gzip header: magic(2) CM FLG MTIME(4) XFL OS XLEN(2)
        if len(head) < 12:
            self._block = b""
            self._uoff = 0
            self._coff = start
            self._block_end_coff = start
            return False
        if head[0:2] != b"\x1f\x8b":
            raise ValueError(f"Invalid gzip magic at file offset {start}")
        xlen = struct.unpack("<H", head[10:12])[0]
        extra = self.fh.read(xlen)
        # parse BC extra field for bsize
        bsize = None
        i = 0
        while i + 4 <= len(extra):
            si1, si2, slen = extra[i], extra[i + 1], struct.unpack("<H", extra[i + 2:i + 4])[0]
            if si1 == 66 and si2 == 67 and slen == 2:
                bsize = struct.unpack("<H", extra[i + 4:i + 6])[0] + 1
            i += 4 + slen
        if bsize is None:
            raise ValueError("Missing BC extra field: not a BGZF block")
        rest_len = bsize - 12 - xlen
        rest = self.fh.read(rest_len)
        deflate_bytes = rest[:-8]
        crc, isize = struct.unpack("<LL", rest[-8:])
        raw = zlib.decompress(deflate_bytes, -15)
        if (zlib.crc32(raw) & 0xFFFFFFFF) != crc or isize != len(raw):
            raise ValueError("BGZF CRC/size mismatch")
        self._coff = start
        self._block = raw
        self._uoff = 0
        self._block_end_coff = start + bsize
        return True

    def tell_virtual(self) -> int:
        return (self._coff << 16) | self._uoff

    def seek_virtual(self, voffset: int):
        coff, uoff = voffset >> 16, voffset & 0xFFFF
        if coff != self._coff or not self._block:
            self.fh.seek(coff)
            self._read_block()
        self._uoff = uoff

    def read(self, n: int) -> bytes:
        out = bytearray()
        while len(out) < n:
            if self._uoff >= len(self._block):
                if not self._read_block():
                    break
            chunk = self._block[self._uoff:self._uoff + (n - len(out))]
            self._uoff += len(chunk)
            out += chunk
        return bytes(out)

    def peek_eof(self) -> bool:
        return self.fh.tell() >= self._file_size() and self._uoff >= len(self._block)

    def _file_size(self):
        cur = self.fh.tell()
        self.fh.seek(0, 2)
        size = self.fh.tell()
        self.fh.seek(cur)
        return size

    def read_until_eof(self) -> Iterator[bytes]:
        """Yield raw decompressed byte stream (chunked)."""
        while True:
            if self._uoff >= len(self._block):
                if not self._read_block():
                    return
            data = self._block[self._uoff:]
            self._uoff = len(self._block)
            if data:
                yield data

    def close(self):
        self.fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
