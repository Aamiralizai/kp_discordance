#!/usr/bin/env bash
# test_run_all_gate.sh — a run that skipped a required stage must not exit 0.
#
# WHY THIS EXISTS
# run_all.sh counted a missing external table in WARN and nowhere else, while
# the exit gate tested SKIPPED_STEPS. A rebuild that never ran the external
# validation therefore finished with status 0 and looked like a complete
# rebuild. The bug lived in the one branch an ordinary run never takes, so
# nothing caught it. This test takes that branch deliberately.
#
# It sources run_all.sh with RUN_ALL_LIB_ONLY=1, which defines the helpers and
# returns before any analysis runs, so the test needs no data and no network.
#
# USAGE
#     bash tests/test_run_all_gate.sh
set -uo pipefail
cd "$(dirname "$0")/.."
RUN_ALL=run_all.sh
fails=0
ok()  { printf '  ok   %s\n' "$1"; }
bad() { printf '  FAIL %s\n' "$1"; fails=$((fails+1)); }

# 1. the gate's behaviour, called directly
gate_status() {   # $1 = how many stages were skipped, $2 = ALLOW_SKIPPED
  (
    set -uo pipefail
    RUN_ALL_LIB_ONLY=1
    export RUN_ALL_LIB_ONLY
    # shellcheck disable=SC1090
    . "./$RUN_ALL" >/dev/null 2>&1
    SKIPPED_STEPS=()
    [ "$1" -gt 0 ] && SKIPPED_STEPS=("external validation (input absent)")
    ALLOW_SKIPPED=$2
    finish_gate >/dev/null 2>&1
  )
  echo $?
}
[ "$(gate_status 1 0)" = "3" ] && ok "a skipped required stage exits 3" \
                              || bad "a skipped required stage did not exit 3"
[ "$(gate_status 1 1)" = "0" ] && ok "ALLOW_SKIPPED=1 exits 0 for a partial run" \
                              || bad "ALLOW_SKIPPED=1 did not exit 0"
[ "$(gate_status 0 0)" = "0" ] && ok "a run with nothing skipped exits 0" \
                              || bad "a clean run did not exit 0"

# 2. the branch that caused the bug: reporting a skip without recording it
#    leaves the gate blind, so every "SKIPPED" message must record a step.
# A stage that skips itself prints a message in the house format
#     echo "   SKIPPED: ..."
# and must record the stage before its branch closes. finish_gate's own
# wording is not a stage skip, so the scan stops at its definition.
blind=$(awk '
  /^finish_gate\(\) \{/ { infn = 1; next }
  infn                   { if ($0 == "}") infn = 0; next }
  /SKIPPED_STEPS\+=/     { armed = 0; next }
  /^[[:space:]]*echo "[[:space:]]+SKIPPED/ { if (!armed) { armed = 1; at = NR } next }
  /^fi$/                 { if (armed) { print at; armed = 0 } }
' "$RUN_ALL")
if [ -z "$blind" ]; then
  ok "every reported skip is recorded for the exit gate"
else
  bad "a skip reported at line(s) $blind never reaches SKIPPED_STEPS"
fi

# 3. sourcing the script must not touch the filesystem. It used to stamp a
#    second, empty log and repoint latest.log at it, in the middle of the run
#    that was writing the real one, because this test sources the script.
sandbox=$(mktemp -d)
cp -r "$RUN_ALL" tests "$sandbox"/ 2>/dev/null
(
  cd "$sandbox"
  RUN_ALL_LIB_ONLY=1 . "./$RUN_ALL" >/dev/null 2>&1
) || true
if [ -z "$(find "$sandbox" -mindepth 1 -name results -o -mindepth 1 -name '*.log' 2>/dev/null)" ]; then
  ok "sourcing the script creates no files"
else
  bad "sourcing the script created $(find "$sandbox" -mindepth 1 -name results -o -name '*.log' | tr '\n' ' ')"
fi
rm -rf "$sandbox"

# 4. the regression tests must be able to stop the run. The script does not use
#    set -e, so a bare test call cannot.
if awk '/say "0\. regression tests/ { z = 1; next }
        z && /^(must |opt )/ { good++ }
        z && /^(python|bash) tests\// { bare++ }
        z && /^say / { exit }
        END { exit (bare > 0) }' "$RUN_ALL"; then
  ok "stage 0 runs its tests through a wrapper that can abort"
else
  bad "a stage 0 test is called bare, so a failure cannot stop the run"
fi

# 5. the external stage stays guarded by the presence of its input
grep -q 'if \[ -f data/external/ncbi_kp_ast.tsv \]' "$RUN_ALL" \
  && ok "external validation is guarded by its input file" \
  || bad "the external validation guard has changed or gone"

echo
if [ "$fails" -eq 0 ]; then
  echo "test_run_all_gate: all checks passed"
else
  echo "test_run_all_gate: $fails check(s) FAILED"
fi
exit $(( fails > 0 ))
