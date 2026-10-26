#!/usr/bin/env python3
"""01d_redownload_verified.py — verified re-download of genomes from BV-BRC.

WHY THIS EXISTS
The original FTPS download (scripts/01_download_bvbrc.sh, wget -q) silently
produced truncated .fna files for about 47% of the collection. The transfer stops
mid-sequence and wget still exits 0. A truncated FASTA still parses, so neither
the structural check nor the BV-BRC-deposited QC metrics (which describe the full
genome) detected it. Kleborate then typed the truncated files, so genes on the
missing contigs (often plasmid-borne carbapenemase and ESBL genes) were missed.

WHAT IT DOES
Downloads each genome from the BV-BRC HTTPS data API and ACCEPTS it only if
  total sequence length == BV-BRC genome_length
  number of records     == BV-BRC contigs
  the file ends with a newline (no mid-sequence cut).
Accepted files are written atomically (<id>.fna.part -> <id>.fna); a genome that
fails is retried in later rounds. Original files are never touched.

USAGE (from the project root)
  python3 scripts/01d_redownload_verified.py IDS.txt OUTDIR
      [--workers 2] [--rounds 10] [--mode single|bulk] [--batch 8] [--budget SECONDS]

  IDS.txt  one BV-BRC genome id per line (e.g. results/qc_concordance/truncated_ids.txt)
  OUTDIR   e.g. data/raw/genomes_redownload

  --workers  parallel connections in single mode (default 2; use 1 if transfers stall)
  --rounds   retry rounds over still-missing genomes (default 10)
  --mode     single = one genome per request (default, most robust)
             bulk   = --batch genomes per request (faster on a good connection)
  --budget   stop after this many seconds (default: run to completion)

OUTPUTS (in OUTDIR)
  <id>.fna         verified complete genome
  manifest.tsv     one row per attempt: genome_id, status, expected and received length/contigs
  failed_ids.txt   genomes still not verified at the end (empty when everything succeeded)
"""
import argparse, json, os, re, shutil, subprocess, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

API = "https://www.bv-brc.org/api"
UA = "kp-discordance/1.0 (verified redownload)"
HAVE_CURL = shutil.which("curl") is not None
HDR = re.compile(rb"\|\s*(\d+\.\d+)\]\s*$")
ACC = re.compile(rb"^accn\|(\d+\.\d+)\.con\.")


def http_get(url, accept, timeout):
    """Return response bytes (b'' on failure). Uses curl when available."""
    if HAVE_CURL:
        r = subprocess.run(["curl", "-sS", "-m", str(timeout), "-H", f"Accept: {accept}",
                            "-A", UA, url], capture_output=True)
        return r.stdout if r.returncode == 0 else b""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": accept})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except Exception:
        return b""


def fetch_meta(ids):
    """BV-BRC genome_length and contigs for each id."""
    meta = {}
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        q = (f"in(genome_id,({','.join(chunk)}))"
             f"&select(genome_id,genome_length,contigs,genome_status)&limit(1000)")
        for attempt in range(4):
            raw = http_get(f"{API}/genome/?{q}", "application/json", 120)
            try:
                for row in json.loads(raw):
                    meta[row["genome_id"]] = row
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
    return meta


def parse(data):
    recs = data.split(b">")[1:]
    length = sum(len(r[r.find(b"\n") + 1:].replace(b"\n", b"").replace(b"\r", b"")) for r in recs)
    return len(recs), length


def accept(gid, data, meta, outdir):
    m = meta.get(gid)
    if not m:
        return (gid, "no_bvbrc_record", "", "", "", "")
    exp_len, exp_ctg = int(m.get("genome_length") or 0), int(m.get("contigs") or 0)
    n, length = parse(data) if data else (0, 0)
    if data and length == exp_len and n == exp_ctg and data.endswith(b"\n"):
        part = os.path.join(outdir, gid + ".fna.part")
        with open(part, "wb") as fh:
            fh.write(data)
        os.replace(part, os.path.join(outdir, gid + ".fna"))
        return (gid, "ok", exp_len, exp_ctg, length, n)
    return (gid, f"mismatch:{n}ctg/{length}bp" if data else "no_data", exp_len, exp_ctg, length, n)


def one(gid, meta, outdir):
    url = f"{API}/genome_sequence/?eq(genome_id,{gid})&limit(25000)&http_accept=application/dna+fasta"
    row = (gid, "no_data", "", "", "", "")
    for attempt in range(3):
        row = accept(gid, http_get(url, "application/dna+fasta", 60), meta, outdir)
        if row[1] == "ok":
            break
        time.sleep(2 * (attempt + 1))
    return [row]


def split_by_genome(data):
    out = {}
    for rec in data.split(b">")[1:]:
        nl = rec.find(b"\n")
        head = rec[:nl] if nl >= 0 else rec
        m = HDR.search(head) or ACC.search(head)
        if m:
            out.setdefault(m.group(1).decode(), []).append(b">" + rec)
    return {g: b"".join(v) for g, v in out.items()}


def bulk(gids, meta, outdir):
    url = (f"{API}/genome_sequence/?in(genome_id,({','.join(gids)}))"
           f"&limit(25000)&http_accept=application/dna+fasta")
    per = split_by_genome(http_get(url, "application/dna+fasta", 180))
    return [accept(g, per.get(g, b""), meta, outdir) for g in gids]


def run_step(a, ids, have, man):
    """Time-safe single pass for short shells: sequential bulk requests, each
    with a timeout that cannot overrun the budget. Re-run until FINISHED."""
    t0 = time.time()
    todo = [g for g in ids if not have(g)]
    if not todo:
        print(f"FINISHED: {len(ids)}/{len(ids)} verified complete; 0 still missing", flush=True)
        open(os.path.join(a.outdir, "failed_ids.txt"), "w").close()
        return 0
    cand = todo[:a.max_per_call]
    meta = fetch_meta(cand)
    ok = bad = 0
    i = 0
    with open(man, "a") as mf:
        while i < len(cand):
            left = a.budget - (time.time() - t0)
            if left < 15:
                break
            batch = cand[i:i + a.batch]; i += len(batch)
            url = (f"{API}/genome_sequence/?in(genome_id,({','.join(batch)}))"
                   f"&limit(25000)&http_accept=application/dna+fasta")
            data = http_get(url, "application/dna+fasta", int(min(a.timeout, left - 5)))
            per = split_by_genome(data)
            for g in batch:
                row = accept(g, per.get(g, b""), meta, a.outdir)
                mf.write("\t".join(str(x) for x in row) + "\n")
                ok += row[1] == "ok"; bad += row[1] != "ok"
            mf.flush()
    total = sum(map(have, ids))
    print(f"step: +{ok} verified, {bad} retry-later; {total}/{len(ids)} verified; "
          f"{time.time() - t0:.0f}s", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ids"); ap.add_argument("outdir")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--rounds", type=int, default=10)
    ap.add_argument("--mode", choices=["single", "bulk"], default="single")
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--budget", type=float, default=0)
    ap.add_argument("--timeout", type=int, default=40, help="seconds per request in step mode")
    ap.add_argument("--max-per-call", dest="max_per_call", type=int, default=250)
    ap.add_argument("--step", action="store_true",
                    help="one time-safe pass of --budget seconds (for short-lived shells)")
    a = ap.parse_args()

    t0 = time.time()
    os.makedirs(a.outdir, exist_ok=True)
    ids = [l.strip() for l in open(a.ids) if l.strip()]
    man = os.path.join(a.outdir, "manifest.tsv")
    if not os.path.exists(man):
        with open(man, "w") as fh:
            fh.write("genome_id\tstatus\tbvbrc_length\tbvbrc_contigs\tgot_length\tgot_contigs\n")
    have = lambda g: os.path.exists(os.path.join(a.outdir, g + ".fna"))
    out_of_time = lambda: a.budget and time.time() - t0 > a.budget
    if a.step:
        sys.exit(run_step(a, ids, have, man))

    print(f"{len(ids)} genomes requested; {sum(map(have, ids))} already verified; "
          f"mode={a.mode} workers={a.workers} curl={'yes' if HAVE_CURL else 'no (urllib)'}", flush=True)
    for rnd in range(1, a.rounds + 1):
        todo = [g for g in ids if not have(g)]
        if not todo or out_of_time():
            break
        meta = fetch_meta(todo)
        print(f"round {rnd}: {len(todo)} to fetch (BV-BRC records found for {len(meta)})", flush=True)
        if a.mode == "single":
            jobs = [[g] for g in todo]
            worker, nthreads = (lambda grp: one(grp[0], meta, a.outdir)), a.workers
        else:
            jobs = [todo[i:i + a.batch] for i in range(0, len(todo), a.batch)]
            worker, nthreads = (lambda grp: bulk(grp, meta, a.outdir)), 1
        ok_round, done_n = 0, 0
        with open(man, "a") as mf, ThreadPoolExecutor(nthreads) as ex:
            futs = []
            for job in jobs:
                futs.append(ex.submit(worker, job))
            for f in as_completed(futs):
                for row in f.result():
                    mf.write("\t".join(str(x) for x in row) + "\n")
                    ok_round += row[1] == "ok"; done_n += 1
                mf.flush()
                if done_n % 100 < len(jobs[0]):
                    print(f"   {done_n}/{len(todo)} attempted, {ok_round} verified this round", flush=True)
                if out_of_time():
                    for x in futs:
                        x.cancel()
                    break
        print(f"round {rnd} done: {ok_round} verified", flush=True)
        if ok_round == 0:
            time.sleep(60)  # back off if the server is refusing, then try again

    missing = [g for g in ids if not have(g)]
    with open(os.path.join(a.outdir, "failed_ids.txt"), "w") as fh:
        fh.write("\n".join(missing) + ("\n" if missing else ""))
    print(f"FINISHED: {len(ids) - len(missing)}/{len(ids)} verified complete; "
          f"{len(missing)} still missing (see failed_ids.txt); {time.time() - t0:.0f}s", flush=True)
    sys.exit(0 if not missing else 2)


if __name__ == "__main__":
    main()
