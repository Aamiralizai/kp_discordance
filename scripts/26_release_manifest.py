#!/usr/bin/env python3
"""26_release_manifest.py — one audit trail for a release.

WHY THIS EXISTS
This answers the question a reader of the paper will ask: which code, run
against which inputs, produced the numbers in the manuscript? The answer was
spread across several logs from different days, and nothing tied them
together. This writes that link as a single file: the environment, the exact
version of every script, the checksum of every frozen input, the checksum of
every generated table and figure, and which stages the run did not perform.

It records what happened. It does not certify that a run was complete — a run
that skipped a stage is recorded as having skipped it, and run_all.sh's exit
status is the authority on that.

USAGE (from the project root, normally as the last step of run_all.sh)
    python3 scripts/26_release_manifest.py [--label "a name for this release"]

OUTPUT
    results/RELEASE_MANIFEST.txt
    results/RELEASE_MANIFEST.sha256   checksums alone, for sha256sum -c
"""
import argparse
import hashlib
import os
import platform
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRACKED = [
    ("code", ["scripts", "tests", "tools", "config", "assets"],
     (".py", ".sh", ".csv", ".yaml", ".yml", ".svg")),
    ("frozen inputs", ["data/processed", "data/interim", "data/raw/phenotypes",
                       "data/raw/genomes_redownload"], (".tsv", ".csv")),
    ("generated results", ["results/tables", "results/supplementary",
                           "results/figures", "results/qc_concordance",
                           "results/allele_search"],
     (".tsv", ".csv", ".txt", ".json", ".jsonl", ".svg", ".xlsx", ".fa", ".sha256",
      ".md")),
]
SKIP_DIRS = {"__pycache__", ".git", ".ipynb_checkpoints", ".venv"}


def sha256(path, cap=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def walk(rel, exts):
    base = os.path.join(ROOT, rel)
    if not os.path.isdir(base):
        return
    for dirpath, dirnames, files in os.walk(base):
        dirnames[:] = [d for d in sorted(dirnames) if d not in SKIP_DIRS]
        for f in sorted(files):
            if exts and not f.endswith(exts):
                continue
            p = os.path.join(dirpath, f)
            yield os.path.relpath(p, ROOT), p


def git(*args):
    try:
        out = subprocess.run(["git", "-C", ROOT, *args], capture_output=True,
                             text=True, timeout=20)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:
        return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", default="")
    ap.add_argument("--allow-env-drift", action="store_true",
                    help="write the manifest even if the installed packages do "
                         "not match the requirements.txt pins. The mismatch is "
                         "recorded in the manifest either way.")
    a = ap.parse_args()

    stamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    lines = ["kp-multiscale-amr — release manifest", ""]
    if a.label:
        lines.append(f"label        : {a.label}")
    lines += [
        f"written      : {stamp}",
        f"host         : {platform.node()}",
        f"platform     : {platform.platform()}",
        f"python       : {sys.version.split()[0]}",
    ]
    commit, status, branch = git("rev-parse", "HEAD"), git("status", "--porcelain"), git("rev-parse", "--abbrev-ref", "HEAD")
    if commit:
        lines += [f"git commit   : {commit}",
                  f"git branch   : {branch or '(detached)'}",
                  f"git state    : {'clean' if not status else 'MODIFIED, see below'}"]
        if status:
            lines += ["             : " + l for l in status.splitlines()[:20]]
    elif os.path.isdir(os.path.join(ROOT, ".git")):
        lines.append("git commit   : no commit yet, this manifest is written "
                     "into the snapshot it describes")
    else:
        lines.append("git commit   : not a git working tree")

    for name in ("kleborate", "python"):
        v = ""
        try:
            v = subprocess.run([name, "--version"], capture_output=True,
                               text=True, timeout=20).stdout.strip()
        except Exception:
            pass
        if v:
            lines.append(f"{name:<13}: {v}")
    # Probe whatever requirements.txt pins, not a list fixed in this file. A
    # list fixed here goes stale the moment a package is added to the pins, and
    # the missing probe then reads as "not importable", which looks like drift
    # and blocks the release for a package that is in fact installed.
    import_name = {"scikit-learn": "sklearn", "pyyaml": "yaml"}
    req = os.path.join(ROOT, "requirements.txt")
    pinned = []
    if os.path.exists(req):
        with open(req, errors="replace") as fh:
            for raw in fh:
                line = raw.split("#")[0].strip()
                if "==" in line:
                    name, want = (x.strip() for x in line.split("==", 1))
                    pinned.append((name, want, import_name.get(name.lower(),
                                                               name.lower())))
    installed = {}
    import importlib
    for mod_name in {m for _, _, m in pinned} | {"pandas", "numpy", "scipy",
                                                 "statsmodels", "matplotlib",
                                                 "sklearn"}:
        try:
            mod = importlib.import_module(mod_name)
            installed[mod_name] = (getattr(mod, "__version__", None)
                                   or mod.etree.__version__)
        except Exception:
            pass
    shown = ("pandas", "numpy", "scipy", "statsmodels", "matplotlib", "sklearn",
             "pyrodigal", "pyswrd")
    lines.append("packages     : " + ", ".join(f"{m} {installed[m]}"
                                               for m in shown if m in installed))

    # The manifest describes the environment the shipped results came from, so
    # it must not be written from an interpreter that has drifted off the pins.
    # Regenerating the manifest later in an upgraded environment would otherwise
    # rewrite the environment claim without touching a single result.
    drift = []
    if pinned:
        for name, want, mod in pinned:
            have = installed.get(mod)
            if have is None:
                drift.append(f"{name}: pinned {want}, not importable")
            elif have != want:
                drift.append(f"{name}: pinned {want}, installed {have}")
        lines.append("pins         : " + ("installed packages match requirements.txt"
                                          if not drift else
                                          f"MISMATCH against requirements.txt ({len(drift)})"))
        lines += ["             : " + d for d in drift]
        if drift and a.allow_env_drift:
            lines.append("             : written anyway with --allow-env-drift, so the "
                         "packages line above is NOT the environment that produced "
                         "the results")
    else:
        lines.append("pins         : requirements.txt not found, pins unchecked")

    log = os.path.join(ROOT, "results", "logs", "latest_summary.txt")
    if os.path.exists(log):
        lines += ["", "--- step outcomes from the run that produced these results ---"]
        with open(log, errors="replace") as fh:
            keep = False
            for line in fh:
                if line.startswith("--- step outcomes"):
                    keep = True
                    continue
                if keep and line.startswith("---"):
                    break
                if keep and line.strip():
                    lines.append("  " + line.rstrip())

    checksums = []
    for section, dirs, exts in TRACKED:
        lines += ["", f"--- {section} ---"]
        n = 0
        for rel, p in ((r, q) for d in dirs for r, q in walk(d, exts)):
            digest = sha256(p)
            size = os.path.getsize(p)
            lines.append(f"  {digest[:16]}  {size:>12,}  {rel}")
            checksums.append(f"{digest}  {rel}")
            n += 1
        lines.append(f"  ({n} file(s))")

    if drift and not a.allow_env_drift:
        sys.stderr.write(
            "26_release_manifest.py: refusing to write. The installed packages do "
            "not match requirements.txt, so the manifest would describe an "
            "environment that did not produce these results:\n")
        for d in drift:
            sys.stderr.write(f"    {d}\n")
        sys.stderr.write(
            "  Install the pins (pip install -r requirements.txt) and rerun, or "
            "pass --allow-env-drift to record the mismatch in the manifest.\n")
        return 1

    out = os.path.join(ROOT, "results", "RELEASE_MANIFEST.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    sums = os.path.join(ROOT, "results", "RELEASE_MANIFEST.sha256")
    with open(sums, "w") as fh:
        fh.write("\n".join(checksums) + "\n")
    print(f"  {len(checksums):,} files recorded")
    print(f"  wrote {out}")
    print(f"  wrote {sums}  (verify with: sha256sum -c RELEASE_MANIFEST.sha256)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
