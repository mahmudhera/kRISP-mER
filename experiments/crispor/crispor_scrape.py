#!/usr/bin/env python3
from __future__ import annotations

import argparse
import random
import re
import sys
import time
from pathlib import Path
from typing import Optional
from urllib.parse import urljoin, urlparse, parse_qs

import requests
from bs4 import BeautifulSoup
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://crispor.gi.ucsc.edu/"
CRISPOR = urljoin(BASE, "crispor.py")


def die(msg: str, code: int = 2) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    raise SystemExit(code)


def normalize(s: str) -> str:
    return " ".join(s.strip().lower().split())


def get_batch_id_from_url(url: str) -> Optional[str]:
    qs = parse_qs(urlparse(url).query)
    return qs.get("batchId", [None])[0]


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "Mozilla/5.0 crispor-web-scraper/1.3"})

    retry = Retry(
        total=10,
        backoff_factor=1.2,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST", "HEAD"),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


def safe_request(session: requests.Session, method: str, url: str, **kwargs) -> requests.Response:
    tries = kwargs.pop("tries", 8)
    base_sleep = kwargs.pop("base_sleep", 1.5)

    last_exc = None
    for i in range(tries):
        try:
            return session.request(method, url, **kwargs)
        except (requests.exceptions.ConnectionError,
                requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ReadTimeout) as e:
            last_exc = e
            sleep_s = base_sleep * (2 ** i) + random.uniform(0, 0.6)
            time.sleep(min(sleep_s, 20))
    raise requests.exceptions.ConnectionError(f"Failed after retries: {last_exc}")


def fetch_homepage(session: requests.Session) -> str:
    r = safe_request(session, "GET", CRISPOR, timeout=60)
    r.raise_for_status()
    return r.text


def extract_org_and_pam_values(home_html: str, genome_query: str, pam_query: str) -> tuple[str, str]:
    soup = BeautifulSoup(home_html, "html.parser")

    org_select = soup.find("select", attrs={"name": "org"}) or soup.find("select", attrs={"name": "genome"})
    if org_select is None:
        die("Could not find genome dropdown (<select name='org'>).")

    pam_select = soup.find("select", attrs={"name": "pam"})
    if pam_select is None:
        die("Could not find PAM dropdown (<select name='pam'>).")

    gq = normalize(genome_query)
    pq = normalize(pam_query)

    org_value = None
    for opt in org_select.find_all("option"):
        label = normalize(opt.get_text(" ", strip=True))
        if gq in label or gq == normalize(opt.get("value") or ""):
            org_value = opt.get("value")
            break
    if not org_value:
        die("Could not match genome label. Use a more unique substring or the GCA_... id.")

    pam_value = None
    for opt in pam_select.find_all("option"):
        label = normalize(opt.get_text(" ", strip=True))
        val = normalize(opt.get("value") or "")
        if pq in label or pq == val:
            pam_value = opt.get("value")
            break
    if not pam_value:
        die("Could not match PAM. Try 'NGG' or part of label like 'SpCas9'.")

    return org_value, pam_value


def submit_job(session: requests.Session, sequence: str, org_value: str, pam_value: str) -> str:
    payload = {"seq": sequence, "org": org_value, "pam": pam_value, "submit": "SUBMIT"}
    r = safe_request(session, "POST", CRISPOR, data=payload, allow_redirects=True, timeout=120)
    r.raise_for_status()

    bid = get_batch_id_from_url(r.url)
    if bid:
        return bid
    if "batchId" in session.cookies:
        return session.cookies.get("batchId")

    m = re.search(r"batchId=([A-Za-z0-9]+)", r.text)
    if m:
        return m.group(1)

    die("Could not find batchId after submission.")
    raise AssertionError


def looks_like_ready_file(r: requests.Response) -> bool:
    ctype = (r.headers.get("Content-Type") or "").lower()
    if "text/html" in ctype:
        return False
    if r.status_code != 200:
        return False
    return True


def wait_until_ready_by_download(session: requests.Session, batch_id: str, timeout_s: int = 1200, poll_s: int = 8) -> None:
    probe_url = f"{CRISPOR}?batchId={batch_id}&download=guides"
    t0 = time.time()
    while True:
        r = safe_request(session, "GET", probe_url, timeout=120)
        if looks_like_ready_file(r):
            return
        if time.time() - t0 > timeout_s:
            die(f"Timed out waiting for CRISPOR job to finish (batchId={batch_id})")
        time.sleep(poll_s)


def get_download_links_from_results(session: requests.Session, batch_id: str) -> dict[str, str]:
    """
    Fetch results page ONCE (after ready), extract the *actual* download URLs
    for: Guides, Guides, all scores, Off-targets.
    """
    url = f"{CRISPOR}?batchId={batch_id}"
    r = safe_request(session, "GET", url, timeout=180)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    links = {}
    for a in soup.select("a[href]"):
        text = normalize(a.get_text(" ", strip=True))
        href = a["href"]

        if text == "guides":
            links["guides"] = urljoin(BASE, href)
        elif text in ("guides, all scores", "guides, all score", "all scores"):
            links["all_scores"] = urljoin(BASE, href)
        elif text in ("off-targets", "off targets", "off-target"):
            links["off_targets"] = urljoin(BASE, href)

    missing = {"guides", "all_scores", "off_targets"} - set(links.keys())
    if missing:
        die(f"Could not find download link(s) on results page: {missing}")

    return links


def download_url(session: requests.Session, url: str, outpath: Path) -> None:
    r = safe_request(session, "GET", url, stream=True, timeout=300)
    r.raise_for_status()

    ctype = (r.headers.get("Content-Type") or "").lower()
    if "text/html" in ctype:
        txt = r.text[:300].replace("\n", " ")
        die(f"Expected file, got HTML instead. Snippet: {txt}")

    with outpath.open("wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 128):
            if chunk:
                f.write(chunk)


def read_sequence(sequence: Optional[str], sequence_file: Optional[str]) -> str:
    if sequence_file:
        p = Path(sequence_file)
        if not p.exists():
            die(f"Sequence file not found: {p}")
        seq = p.read_text().strip()
    elif sequence:
        seq = sequence.strip()
    else:
        die("Provide --sequence or --sequence-file (unless using --batch-id).")

    seq = re.sub(r"\s+", "", seq)
    if len(seq) < 20:
        die("Sequence looks too short.")
    return seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sequence", help="Raw DNA sequence.")
    ap.add_argument("--sequence-file", help="File containing DNA sequence.")
    ap.add_argument("--genome", help="Genome label substring or GCA_... id.")
    ap.add_argument("--pam", default="NGG", help="PAM value/label substring (default: NGG).")
    ap.add_argument("--batch-id", help="If provided, skips submission.")
    ap.add_argument("--outdir", default="crispor_results")
    ap.add_argument("--timeout", type=int, default=1200)
    ap.add_argument("--poll", type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    session = make_session()
    try:
        if args.batch_id:
            batch_id = args.batch_id.strip()
        else:
            if not args.genome:
                die("Provide --genome if submitting a new job.")
            seq = read_sequence(args.sequence, args.sequence_file)
            home = fetch_homepage(session)
            org_value, pam_value = extract_org_and_pam_values(home, args.genome, args.pam)
            batch_id = submit_job(session, seq, org_value, pam_value)

        print(f"[INFO] batchId = {batch_id}")
        print("[INFO] Waiting until job is ready...")
        wait_until_ready_by_download(session, batch_id, timeout_s=args.timeout, poll_s=args.poll)

        print("[INFO] Ready. Fetching real download links from results page...")
        dl = get_download_links_from_results(session, batch_id)
        print("[INFO] Download URLs:")
        for k, v in dl.items():
            print(f"  - {k}: {v}")

        print("[INFO] Downloading tables...")
        download_url(session, dl["guides"], outdir / f"{batch_id}.guides.xls")
        download_url(session, dl["all_scores"], outdir / f"{batch_id}.guides_all_scores.xls")
        download_url(session, dl["off_targets"], outdir / f"{batch_id}.offtargets.xls")

        print("[INFO] Done. Files saved to:", outdir.resolve())

    finally:
        session.close()


if __name__ == "__main__":
    main()
