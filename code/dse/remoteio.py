"""Remote (HTTP range) BGZF + tabix reader — pure Python, no compiled deps.

Designed for FinnGen R13 summary statistics and GWAS Catalog harmonised files:
each ``.gz`` is BGZF (a series of gzip members) and is paired with a ``.tbi``
tabix index, so an arbitrary genomic window can be fetched with two or three
HTTP range requests instead of downloading a ~760 MB file.

Virtual offsets are ``(compressed_offset << 16) | uncompressed_offset``.
"""

from __future__ import annotations

import gzip
import os
import ssl
import struct
import threading
import time
import urllib.error
import urllib.request
from typing import Dict, List, Optional, Tuple

UA = {"User-Agent": "dse-genomics/1.0 (academic research)"}

# Some mirrors present an incomplete chain to Python's default store; use the
# certifi bundle and fall back to an unverified context only for public data.
try:  # pragma: no cover
    import certifi

    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:  # noqa: BLE001
    SSL_CTX = ssl.create_default_context()
SSL_CTX_INSECURE = ssl._create_unverified_context()


# --------------------------------------------------------------------------
# HTTP range reader with a windowed cache
# --------------------------------------------------------------------------
class RangedHTTP:
    """Random access to a remote file via HTTP Range requests."""

    def __init__(self, url: str, window: int = 1 << 21, timeout: int = 90):
        self.url = url
        self.window = window
        self.timeout = timeout
        self._cache: Dict[int, bytes] = {}
        self._order: List[int] = []
        self._lock = threading.Lock()
        self.n_requests = 0
        self.n_bytes = 0
        self.size = self._head_size()

    def _head_size(self) -> Optional[int]:
        req = urllib.request.Request(self.url, headers=dict(UA), method="HEAD")
        for ctx in (SSL_CTX, SSL_CTX_INSECURE):
            try:
                with urllib.request.urlopen(req, timeout=self.timeout, context=ctx) as r:
                    cl = r.headers.get("Content-Length")
                    return int(cl) if cl else None
            except ssl.SSLError:
                continue
            except Exception:  # noqa: BLE001
                return None
        return None

    def _fetch(self, start: int, end: int) -> bytes:
        last = None
        for ctx in (SSL_CTX, SSL_CTX_INSECURE):
            req = urllib.request.Request(self.url, headers={**UA, "Range": f"bytes={start}-{end - 1}"})
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(req, timeout=self.timeout, context=ctx) as r:
                        d = r.read()
                    self.n_requests += 1
                    self.n_bytes += len(d)
                    return d
                except ssl.SSLError as e:
                    last = e
                    break
                except Exception as e:  # noqa: BLE001
                    last = e
                    time.sleep(1.5 * (attempt + 1))
        raise RuntimeError(f"range fetch failed {self.url} [{start},{end}): {last}")

    def read(self, off: int, n: int) -> bytes:
        """Read ``n`` bytes at absolute offset ``off``."""
        if n <= 0:
            return b""
        out = []
        need = n
        pos = off
        while need > 0:
            wi = pos // self.window
            lo = wi * self.window
            with self._lock:
                buf = self._cache.get(wi)
            if buf is None:
                # grow the request a little to absorb sequential access
                end = min(lo + self.window, self.size or (lo + self.window))
                if self.size and lo >= self.size:
                    break
                buf = self._fetch(lo, end)
                with self._lock:
                    self._cache[wi] = buf
                    self._order.append(wi)
                    if len(self._order) > 384:
                        old = self._order.pop(0)
                        if old != wi:
                            self._cache.pop(old, None)
            rel = pos - lo
            if rel >= len(buf):
                break
            take = min(need, len(buf) - rel)
            out.append(buf[rel:rel + take])
            pos += take
            need -= take
        return b"".join(out)


# --------------------------------------------------------------------------
# BGZF
# --------------------------------------------------------------------------
def bgzf_blocksize(head: bytes) -> Optional[int]:
    """Total block size from an 18-byte BGZF block header."""
    if len(head) < 18 or head[0:2] != b"\x1f\x8b" or not (head[3] & 0x04):
        return None
    xlen = int.from_bytes(head[10:12], "little")
    extra = head[12:12 + xlen]
    p = 0
    while p + 4 <= len(extra):
        si1, si2 = extra[p], extra[p + 1]
        slen = int.from_bytes(extra[p + 2:p + 4], "little")
        if si1 == 66 and si2 == 67:
            return int.from_bytes(extra[p + 4:p + 6], "little") + 1
        p += 4 + slen
    return None


class BGZFRemote:
    """Iterate decompressed BGZF blocks of a remote file."""

    def __init__(self, rf: RangedHTTP):
        self.rf = rf
        self._blen: Dict[int, int] = {}

    def block(self, coff: int) -> Tuple[bytes, int]:
        head = self.rf.read(coff, 18)
        bs = bgzf_blocksize(head)
        if bs is None:
            # plain gzip: fall back to a large slab
            raw = self.rf.read(coff, 1 << 20)
            return gzip.decompress(raw), len(raw)
        raw = self.rf.read(coff, bs)
        return gzip.decompress(raw), bs

    def blocks_from(self, coff: int):
        while True:
            try:
                data, bs = self.block(coff)
            except Exception:  # noqa: BLE001
                return
            if not data:
                return
            yield coff, data, bs
            coff += bs


# --------------------------------------------------------------------------
# tabix index
# --------------------------------------------------------------------------
class TabixIndex:
    def __init__(self, raw: bytes):
        assert raw[:4] == b"TBI\x01", "not a tabix index"
        (self.n_ref, self.fmt, self.col_seq, self.col_beg,
         self.col_end, self.meta, self.skip, l_nm) = struct.unpack("<8i", raw[4:36])
        self.names = raw[36:36 + l_nm].split(b"\x00")
        p = 36 + l_nm
        self.refs: List[dict] = []
        for _ in range(self.n_ref):
            n_bin = struct.unpack("<i", raw[p:p + 4])[0]
            p += 4
            bins: Dict[int, List[Tuple[int, int]]] = {}
            for _ in range(n_bin):
                b, n_chunk = struct.unpack("<Ii", raw[p:p + 8])
                p += 8
                chunks = []
                for _ in range(n_chunk):
                    bg, en = struct.unpack("<QQ", raw[p:p + 16])
                    p += 16
                    chunks.append((bg, en))
                bins[b] = chunks
            n_intv = struct.unpack("<i", raw[p:p + 4])[0]
            p += 4
            ioff = list(struct.unpack("<" + "Q" * n_intv, raw[p:p + 8 * n_intv]))
            p += 8 * n_intv
            self.refs.append({"bins": bins, "linear": ioff})
        self.names = [n.decode() for n in self.names if n]


def reg2bins(beg: int, end: int) -> List[int]:
    """tabix binning scheme (htslib ``reg2bins``, 14/17/20/23/26 bit levels)."""
    if beg < 0:
        beg = 0
    end -= 1
    bins = [0]
    for k, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        lo = k + (beg >> shift)
        hi = k + (end >> shift)
        bins.extend(range(lo, hi + 1))
    return bins


class Tabix:
    def __init__(self, gz_url: str, cache_dir: Optional[str] = None):
        self.gz = RangedHTTP(gz_url)
        tbi_url = gz_url + ".tbi"
        self.idx = TabixIndex(self._load_tbi(tbi_url, cache_dir))
        self.bg = BGZFRemote(self.gz)

    def _load_tbi(self, url: str, cache_dir: Optional[str]) -> bytes:
        key = os.path.basename(url)
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            fp = os.path.join(cache_dir, key)
            if os.path.exists(fp) and os.path.getsize(fp) > 100:
                with open(fp, "rb") as f:
                    return gzip.decompress(f.read())
        rf = RangedHTTP(url)
        parts = []
        off = 0
        while True:
            d = rf._fetch(off, off + (1 << 21))
            if not d:
                break
            parts.append(d)
            off += len(d)
            if rf.size and off >= rf.size:
                break
        raw = gzip.decompress(b"".join(parts))
        if cache_dir:
            with open(os.path.join(cache_dir, key), "wb") as f:
                f.write(b"".join(parts))
        return raw

    def _ref(self, chrom: str) -> Optional[int]:
        c = str(chrom)
        for cand in (c, c.replace("chr", ""), "chr" + c.replace("chr", "")):
            try:
                return self.idx.names.index(cand)
            except ValueError:
                continue
        return None

    def query(self, chrom: str, beg0: int, end0: int, max_lines: int = 5_000_000):
        """Yield decoded lines overlapping the 0-based half-open window."""
        ri = self._ref(chrom)
        if ri is None:
            return
        ref = self.idx.refs[ri]
        bins = reg2bins(beg0, end0)
        chunks: List[Tuple[int, int]] = []
        for b in bins:
            chunks.extend(ref["bins"].get(b, ()))
        if not chunks:
            return
        lin = ref["linear"]
        if lin:
            i = min(beg0 >> 14, len(lin) - 1)
            min_off = lin[i]
            chunks = [c for c in chunks if c[1] > min_off]
        chunks.sort()
        merged: List[List[int]] = []
        for bg_, en in chunks:
            if merged and bg_ <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], en)
            else:
                merged.append([bg_, en])

        seen = 0
        for vbeg, vend in merged:
            coff, uoff = vbeg >> 16, vbeg & 0xFFFF
            carry = b""
            while True:
                try:
                    data, bs = self.bg.block(coff)
                except Exception:  # noqa: BLE001
                    break
                blob = carry + (data[uoff:] if uoff else data)
                cut = blob.rfind(b"\n")
                if cut < 0:
                    carry = blob
                else:
                    for line in blob[:cut].split(b"\n"):
                        if line:
                            yield line.decode("utf-8", "replace")
                            seen += 1
                            if seen >= max_lines:
                                return
                    carry = blob[cut + 1:]
                coff += bs
                uoff = 0
                if (coff << 16) >= vend:
                    break
            if carry:
                yield carry.decode("utf-8", "replace")
