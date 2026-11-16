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
    """Used for fatal errors before the batch loop starts."""
    print(f"[ERROR] {msg}", file=sys.stderr)
    raise SystemExit(code)


def normalize(s: str) -> str:
    return " ".join(s.strip().lower().split())


def get_batch_id_from_url(url: str) -> Optional[str]:
    qs = parse_qs(urlparse(url).query)
    return qs.get("batchId", [None])[0]


def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "Mozilla/5.0 crispor-batch/1.0"})

    retry_kwargs = dict(
        total=10,
        backoff_factor=1.2,
        status_forcelist=(429, 500, 502, 503, 504),
        raise_on_status=False,
    )

    try:
        retry = Retry(allowed_methods=("GET", "POST", "HEAD"), **retry_kwargs)
    except TypeError:
        retry = Retry(method_whitelist=("GET", "POST", "HEAD"), **retry_kwargs)

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
            time.sleep(min(sleep_s, 25))
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
        val = normalize(opt.get("value") or "")
        if gq in label or gq == val:
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


def read_fasta_sequence(fasta_path: Path) -> tuple[str, str]:
    text = fasta_path.read_text().strip().splitlines()
    if not text or not text[0].startswith(">"):
        raise ValueError(f"{fasta_path} is not a valid FASTA (missing header line).")

    name = text[0][1:].strip()
    seq = "".join(line.strip() for line in text[1:] if line.strip())
    seq = re.sub(r"\s+", "", seq).upper()
    if len(seq) < 20:
        raise ValueError(f"{fasta_path} sequence too short.")
    return name, seq


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

    raise RuntimeError("Could not find batchId after submission.")


def looks_like_ready_file(r: requests.Response) -> bool:
    ctype = (r.headers.get("Content-Type") or "").lower()
    if "text/html" in ctype:
        return False
    return r.status_code == 200


def wait_until_ready_by_probe(session: requests.Session, batch_id: str, timeout_s: int = 1200, poll_s: int = 8) -> None:
    probe_url = f"{CRISPOR}?batchId={batch_id}&download=guides"
    t0 = time.time()
    while True:
        r = safe_request(session, "GET", probe_url, timeout=120)
        if looks_like_ready_file(r):
            return
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f"Timed out waiting for CRISPOR job to finish (batchId={batch_id})")
        time.sleep(poll_s)


def get_download_links_from_results(session: requests.Session, batch_id: str) -> dict[str, str]:
    url = f"{CRISPOR}?batchId={batch_id}"
    r = safe_request(session, "GET", url, timeout=180)
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "html.parser")

    links: dict[str, str] = {}
    for a in soup.select("a[href]"):
        txt = normalize(a.get_text(" ", strip=True))
        href = a["href"]
        if txt == "guides":
            links["guides"] = urljoin(BASE, href)
        elif txt in ("guides, all scores", "guides, all score", "all scores"):
            links["all_scores"] = urljoin(BASE, href)
        elif txt in ("off-targets", "off targets", "off-target"):
            links["off_targets"] = urljoin(BASE, href)

    missing = {"guides", "all_scores", "off_targets"} - set(links.keys())
    if missing:
        # Raise standard Exception so main loop catches it without halting the program
        raise RuntimeError(f"Could not find download link(s) on results page: {missing}. Sequence likely has 0 valid guides.")
    return links


def download_url(session: requests.Session, url: str, outpath: Path) -> None:
    r = safe_request(session, "GET", url, stream=True, timeout=300)
    r.raise_for_status()

    ctype = (r.headers.get("Content-Type") or "").lower()
    if "text/html" in ctype:
        txt = r.text[:300].replace("\n", " ")
        raise RuntimeError(f"Expected file, got HTML instead. Snippet: {txt}")

    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("wb") as f:
        for chunk in r.iter_content(chunk_size=1024 * 128):
            if chunk:
                f.write(chunk)


def run_one_target(
    session: requests.Session,
    fasta_path: Path,
    org_value: str,
    pam_value: str,
    out_root: Path,
    timeout: int,
    poll: int,
    sleep_between: float,
) -> bool:
    target_folder = out_root / fasta_path.stem
    guides_fp = target_folder / "guides.xls"
    all_scores_fp = target_folder / "guides_all_scores.xls"
    offt_fp = target_folder / "offtargets.xls"

    if all(p.exists() and p.stat().st_size > 0 for p in (guides_fp, all_scores_fp, offt_fp)):
        print(f"[SKIP] {fasta_path.name} (already downloaded)")
        return True

    name, seq = read_fasta_sequence(fasta_path)
    print(f"[RUN] {fasta_path.name} | header='{name}' | len={len(seq)}")

    batch_id = submit_job(session, seq, org_value, pam_value)
    print(f"  [INFO] batchId={batch_id}")

    wait_until_ready_by_probe(session, batch_id, timeout_s=timeout, poll_s=poll)
    dl = get_download_links_from_results(session, batch_id)

    download_url(session, dl["guides"], guides_fp)
    download_url(session, dl["all_scores"], all_scores_fp)
    download_url(session, dl["off_targets"], offt_fp)

    print(f"  [OK] saved -> {target_folder}")
    time.sleep(max(0.0, sleep_between))
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets-dir", required=True, help="Folder containing *.fasta targets")
    ap.add_argument("--outdir", default="crispor_outputs", help="Output root folder")
    ap.add_argument("--genome", required=True, help="Genome label substring or GCA_... id")
    ap.add_argument("--pam", default="NGG", help="PAM, default NGG")
    ap.add_argument("--timeout", type=int, default=1200, help="Max seconds to wait per job")
    ap.add_argument("--poll", type=int, default=8, help="Polling interval seconds")
    ap.add_argument("--sleep-between", type=float, default=2.0, help="Sleep seconds between targets")
    ap.add_argument("--limit", type=int, default=0, help="Process only first N targets (0 = all)")
    args = ap.parse_args()

    targets_dir = Path(args.targets_dir)
    if not targets_dir.exists():
        die(f"targets-dir not found: {targets_dir}")

    out_root = Path(args.outdir)
    out_root.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(targets_dir.glob("*.fasta"))
    if not fasta_files:
        die(f"No .fasta files found in {targets_dir}")

    if args.limit and args.limit > 0:
        fasta_files = fasta_files[:args.limit]

    session = make_session()
    try:
        home = fetch_homepage(session)
        org_value, pam_value = extract_org_and_pam_values(home, args.genome, args.pam)
        print(f"[INFO] Using org={org_value} pam={pam_value}")
        print(f"[INFO] Targets: {len(fasta_files)} | outdir={out_root.resolve()}")

        ok = 0
        fail = 0
        for fp in fasta_files:
            try:
                run_one_target(
                    session=session,
                    fasta_path=fp,
                    org_value=org_value,
                    pam_value=pam_value,
                    out_root=out_root,
                    timeout=args.timeout,
                    poll=args.poll,
                    sleep_between=args.sleep_between,
                )
                ok += 1
            except Exception as e:
                fail += 1
                # Prints the error and continues the loop!
                print(f"[FAIL] {fp.name}: {e}", file=sys.stderr)
                time.sleep(5)

        print(f"[DONE] ok={ok} fail={fail}")

    finally:
        session.close()


if __name__ == "__main__":
    main()