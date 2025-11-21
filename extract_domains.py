#!/usr/bin/env python3
"""
extract_domains.py
Usage:
  python3 extract_domains.py --ids ids.txt --outfasta out.fasta --outtsv out.tsv --mode by-domains \
    --domain1 "La motif" --domain2 "RRM" --workers 5
or
  python3 extract_domains.py --ids ids.txt --outfasta out.fasta --outtsv out.tsv --mode by-order
"""
import requests, time, sys, argparse, json
from concurrent.futures import ThreadPoolExecutor, as_completed

API_BASE = "https://rest.uniprot.org/uniprotkb/{}.json"
HEADERS = {"Accept": "application/json"}

def fetch_json(uid, max_retries=3, backoff=1.0):
    url = API_BASE.format(uid)
    for i in range(max_retries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            if r.status_code == 200:
                return r.json()
            else:
                # handle 404 or rate limit
                if r.status_code == 404:
                    return {"_error": f"404 not found for {uid}"}
                # other: retry
                time.sleep(backoff*(i+1))
        except Exception as e:
            time.sleep(backoff*(i+1))
    return {"_error": f"failed after retries for {uid}"}

def extract_by_order(j):
    seq = j.get("sequence", {}).get("value")
    if not seq:
        return None, "no_sequence"
    domains = [f for f in j.get("features", []) if f.get("type") == "Domain"]
    if len(domains) < 2:
        return None, f"less_than_2_domains:{len(domains)}"
    # take start of first, end of second
    s = domains[0]["location"]["start"]["value"]
    e = domains[1]["location"]["end"]["value"]
    return seq[int(s)-1:int(e)], f"{s}-{e}"

def matches_domain(feature, pattern):
    """pattern can be substring of description or pfam id in dbReferences"""
    desc = feature.get("description", "") or ""
    if pattern.lower() in desc.lower():
        return True
    # check dbReferences for Pfam
    for db in feature.get("evidences", []) + feature.get("dbReferences", []):
        pass
    # Pfam info is usually in feature['dbReferences'] where db=='Pfam' and id like 'PF01250'
    for dbr in feature.get("dbReferences", []):
        if dbr.get("type","").lower() == "pfam":
            if pattern.upper() in dbr.get("id","").upper() or pattern.lower() in dbr.get("id","").lower():
                return True
    return False

def extract_by_domain_names(j, domain1_pat, domain2_pat):
    seq = j.get("sequence", {}).get("value")
    if not seq:
        return None, "no_sequence"
    domains = [f for f in j.get("features", []) if f.get("type") == "Domain"]
    d1 = None; d2 = None
    for f in domains:
        if d1 is None and matches_domain(f, domain1_pat):
            d1 = f
            continue
        if d2 is None and matches_domain(f, domain2_pat):
            d2 = f
            if d1 is not None:
                break
    if d1 and d2:
        s = d1["location"]["start"]["value"]
        e = d2["location"]["end"]["value"]
        return seq[int(s)-1:int(e)], f"{s}-{e}"
    return None, f"pattern_not_found d1:{bool(d1)} d2:{bool(d2)}"

def worker(uid, mode, domain1, domain2, delay=0.2):
    j = fetch_json(uid)
    if "_error" in j:
        return uid, None, j["_error"]
    if mode == "by-order":
        seqfrag, status = extract_by_order(j)
    else:
        seqfrag, status = extract_by_domain_names(j, domain1, domain2)
    # be polite
    time.sleep(delay)
    return uid, seqfrag, status

def main(args):
    ids = [line.strip() for line in open(args.ids) if line.strip()]
    outf = open(args.outfasta, "w")
    outt = open(args.outtsv, "w")
    outt.write("id\tstatus\n")
    workers = min(args.workers, max(1, len(ids)))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(worker, uid, args.mode, args.domain1, args.domain2, args.delay): uid for uid in ids}
        for fut in as_completed(futs):
            uid = futs[fut]
            try:
                uid, seqfrag, status = fut.result()
            except Exception as e:
                outt.write(f"{uid}\texception:{e}\n")
                continue
            if seqfrag:
                # wrap to 60
                outf.write(f">{uid}|{status}\n")
                for i in range(0, len(seqfrag), 60):
                    outf.write(seqfrag[i:i+60] + "\n")
                outt.write(f"{uid}\tok:{status}\n")
            else:
                outt.write(f"{uid}\t{status}\n")
    outf.close(); outt.close()

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--ids", required=True, help="file with UniProt accessions, one per line")
    p.add_argument("--outfasta", default="out.fasta")
    p.add_argument("--outtsv", default="out.tsv")
    p.add_argument("--mode", choices=["by-order","by-domains"], default="by-order")
    p.add_argument("--domain1", default="HTH La-type RNA-binding", help="pattern for domain1 matching (used with by-domains)")
    p.add_argument("--domain2", default="RRM", help="pattern for domain2 matching (used with by-domains)")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--delay", type=float, default=0.2, help="per-request delay (seconds)")
    args = p.parse_args()
    main(args)

