"""Shared HTTP helpers: requests session with retries and TLS fallback."""

from __future__ import annotations

import time
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

_session: Optional[requests.Session] = None
_session_noverify: Optional[requests.Session] = None


def session(verify: bool = True) -> requests.Session:
    global _session, _session_noverify
    if verify:
        if _session is None:
            s = requests.Session()
            r = Retry(total=4, backoff_factor=1.2, status_forcelist=[429, 500, 502, 503, 504],
                      allowed_methods=["GET", "HEAD"])
            s.mount("https://", HTTPAdapter(max_retries=r))
            s.headers.update({"User-Agent": "dse-genomics/1.0 (academic)"})
            _session = s
        return _session
    if _session_noverify is None:
        s = requests.Session()
        r = Retry(total=3, backoff_factor=1.0, status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET", "HEAD"])
        s.mount("https://", HTTPAdapter(max_retries=r))
        s.headers.update({"User-Agent": "dse-genomics/1.0 (academic)"})
        _session_noverify = s
    return _session_noverify


def get(url: str, timeout: int = 90, params: dict | None = None, verify: bool = True):
    """GET returning (status, text). Falls back to unverified TLS once."""
    for attempt in range(3):
        try:
            r = session(verify).get(url, timeout=timeout, params=params)
            return r.status_code, r.text
        except requests.exceptions.SSLError:
            try:
                r = session(False).get(url, timeout=timeout, params=params)
                return r.status_code, r.text
            except Exception:  # noqa: BLE001
                time.sleep(1.5)
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
    return -1, ""


def get_json(url: str, timeout: int = 90, params: dict | None = None):
    st, txt = get(url, timeout=timeout, params=params)
    if st != 200 or not txt:
        return None
    try:
        return requests.models.complexjson.loads(txt) if hasattr(requests.models, "complexjson") else __import__("json").loads(txt)
    except Exception:  # noqa: BLE001
        return None
