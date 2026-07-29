#!/usr/bin/env python
"""Regression tests for breakpoint classification and genotype parsing.

Each test encodes an error that was made and found in review. Run before and
after any change to classification, parsing or table I/O:

    python tests/test_classification.py
"""
import importlib.util
import io
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

spec = importlib.util.spec_from_file_location(
    "rl", ROOT / "scripts" / "08b_discordance_recordlevel.py")
rl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rl)
from kp_defs import carb_family, porin_status  # noqa: E402

ri, cl, ix = rl.record_interval, rl.classify_interval, rl.intersect
fails = []


def check(name, got, want):
    ok = got == want
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {got!r}"
          + ("" if ok else f"  (expected {want!r})"))
    if not ok:
        fails.append(name)


print("EUCAST meropenem  S <= 2, R > 8")
check(">8  -> R", cl(ri(8, ">"), 2, 8, "<=", ">"), "R")
check(">=8 -> ambiguous", cl(ri(8, ">="), 2, 8, "<=", ">"), "ambiguous")
check("exact 8 -> I", cl(ri(8, ""), 2, 8, "<=", ">"), "I")
check("exact 16 -> R", cl(ri(16, ""), 2, 8, "<=", ">"), "R")
check("<=2 -> S", cl(ri(2, "<="), 2, 8, "<=", ">"), "S")
check("<=4 -> ambiguous", cl(ri(4, "<="), 2, 8, "<=", ">"), "ambiguous")

print("\nEUCAST imipenem  S <= 2, R > 4")
check(">4  -> R", cl(ri(4, ">"), 2, 4, "<=", ">"), "R")
check("exact 4 -> I", cl(ri(4, ""), 2, 4, "<=", ">"), "I")
check("exact 8 -> R", cl(ri(8, ""), 2, 4, "<=", ">"), "R")

print("\nCLSI  S <= 1, R >= 4  (comparator differs from EUCAST)")
check("exact 4 -> R", cl(ri(4, ""), 1, 4, "<=", ">="), "R")
check("exact 2 -> I", cl(ri(2, ""), 1, 4, "<=", ">="), "I")
check("<=1 -> S", cl(ri(1, "<="), 1, 4, "<=", ">="), "S")

print("\ninterval intersection")
check("<=1 & exact 4 -> contradictory", ix(ri(1, "<="), ri(4, "")), None)
check("<=8 & >=2 -> [2,8]", ix(ri(8, "<="), ri(2, ">=")), (2.0, 8.0, True, True))

print("\nboolean round-trip through TSV (dtype=str)")
d = pd.read_csv(io.StringIO("lo_inc\nFalse\n"), sep="\t", dtype=str)
check("raw string is truthy (the bug)", bool(d["lo_inc"][0]), True)
d["lo_inc"] = d["lo_inc"].astype(str).str.strip().str.lower().isin({"true", "1"})
check("after conversion", bool(d["lo_inc"][0]), False)

print("\nporin HGVS parsing")
check("loop-3 GD", porin_status("OmpK36:p.134_135insGlyAsp"), "ompK36_loop3")
check("loop-3 TD", porin_status("OmpK36:p.136_137insThrAsp"), "ompK36_loop3")
check("non-loop-3 insertion", porin_status("OmpK36:p.220_221insAlaGly"),
      "ompK36_other_variant")
check("c.25C>T is reduced expression, NOT a stop",
      porin_status("OmpK36:c.25C>T"), "ompK36_reduced_expression")
check("frameshift", porin_status("OmpK36:p.Ala183fs"), "ompK36_loss")
check("OmpK35 only", porin_status("OmpK35:p.Glu42fs"), "ompK35_loss")
check("none", porin_status("-"), "none")

print("\ncarbapenemase family assignment")
check("OXA-162 in OXA-48-like", carb_family("OXA-162"), "OXA-48-like")
check("OXA-181 in OXA-48-like", carb_family("OXA-181"), "OXA-48-like")
check("OXA-23 not folded in", carb_family("OXA-23"), "OXA-other")
check("combined label", carb_family("KPC-2;NDM-1"), "KPC+NDM")

print()
if fails:
    print(f"{len(fails)} TEST(S) FAILED: {fails}")
    sys.exit(1)
print("all tests passed")
