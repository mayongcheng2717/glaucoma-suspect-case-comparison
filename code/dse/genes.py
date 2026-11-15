"""Gene annotation via Ensembl REST (cached)."""

from __future__ import annotations

import json
import os
import time
import urllib.parse
import urllib.request
from typing import Dict, List, Optional

CACHE = r"D:\eye100\cache\ensembl"
BASE = "https://rest.ensembl.org"
UA = {"User-Agent": "dse/1.0", "Content-Type": "application/json"}


def _get(url: str, tries: int = 3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=45) as r:
                return json.load(r)
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (k + 1))
    return None


def gene(symbol: str) -> Optional[dict]:
    os.makedirs(CACHE, exist_ok=True)
    fp = os.path.join(CACHE, f"gene_{symbol.replace('/', '_')}.json")
    if os.path.exists(fp):
        with open(fp, encoding="utf-8") as f:
            return json.load(f)
    d = _get(f"{BASE}/lookup/symbol/homo_sapiens/{urllib.parse.quote(symbol)}?expand=0")
    if d:
        with open(fp, "w", encoding="utf-8") as f:
            json.dump(d, f)
    return d


def genes_in(chrom: str, start: int, end: int) -> List[dict]:
    chrom = str(chrom).replace("chr", "")
    os.makedirs(CACHE, exist_ok=True)
    fp = os.path.join(CACHE, f"region_{chrom}_{start}_{end}.json")
    if os.path.exists(fp):
        with open(fp, encoding="utf-8") as f:
            return json.load(f)
    d = _get(f"{BASE}/overlap/region/human/{chrom}:{start}-{end}?feature=gene;content-type=application/json")
    d = d or []
    keep = [{"symbol": g.get("external_name") or g.get("id"),
             "ensembl": g.get("id"), "start": g.get("start"), "end": g.get("end"),
             "strand": g.get("strand"), "biotype": g.get("biotype")}
            for g in d if g.get("biotype") == "protein_coding"]
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(keep, f)
    return keep


def window(symbol: str, flank: int = 200_000):
    g = gene(symbol)
    if not g:
        return None
    return str(g["seq_region_name"]), int(g["start"]) - flank, int(g["end"]) + flank, g
