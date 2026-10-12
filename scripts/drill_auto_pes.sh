#!/usr/bin/env bash
# drill_auto_pes.sh - the Tier 1 auto-PES container drill (criterion family):
#   positive legs (first-make green, woven content, adiabatic default run with
#   the arm witness ABSENT, ELEC_METHOD = auto_tdhf run with the witness
#   PRESENT, R_TERM both branches, double-run determinism) + the runtime
#   coexistence abort (reference + auto = STOP 1) + the zero-touch close.
# Leaves the tree as found (trap restores); exits nonzero on any failed leg.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/src"; BLD="$SRC/build/gfortran"
FIX="$ROOT/systems/auto_sample"; SHELF="$ROOT/systems/inactive/auto_sample"
REF="$ROOT/systems/reference"; REFSHELF="$ROOT/systems/inactive/reference"
DRILL="$BLD/drill_case"
pass=0; fail=0
ok()  { echo "PASS: $1"; pass=$((pass+1)); }
bad() { echo "FAIL: $1"; fail=$((fail+1)); }

restore_tree() {
   [ -d "$FIX" ] && mv "$FIX" "$SHELF"
   [ -d "$REFSHELF" ] && mv "$REFSHELF" "$REF"
   rm -rf "$DRILL"
}
trap restore_tree EXIT

cd "$SRC"

# ---- leg 1: drag in, FIRST make green ----
mv "$SHELF" "$FIX"
if make >/tmp/drill_make.log 2>&1; then ok "first make green (drag-in, no second make)"
else bad "first make"; tail -20 /tmp/drill_make.log; exit 1; fi

# ---- leg 2: woven content ----
while IFS= read -r pat; do
   grep -qF "$pat" "$BLD/registry_gen.f90" && ok "woven: $pat" || bad "woven: $pat"
done <<'PATS'
use auto_sample_pot, only: pes, keys, term
use auto_sample_io, only: load, classify
use auto_sample_met, only: prop, match, sync, arm, method_name
call method_reg(method_name(), prop, match, sync)
if (trim(electronic%method) == method_name()) call arm()
call auto_slots()
use container_sched, only: container_sched_bind
call container_sched_bind(auto_keys, &
auto_force, &
PATS

# ---- leg 3: runtime coexistence abort (reference + auto = STOP 1 at assembly;
#      registry_gen_keys runs before the input read, so no input file needed) ----
"$BLD/venus.e" >/dev/null 2>/tmp/drill_coexist.err
if grep -q "slot is already filled" /tmp/drill_coexist.err; then ok "coexistence STOP 1"
else bad "coexistence STOP 1"; cat /tmp/drill_coexist.err; fi

# ---- leg 4: shelve reference, rebuild, scratch-case runs ----
mv "$REF" "$REFSHELF"
make >/dev/null 2>&1 || { bad "shelved-reference rebuild"; exit 1; }
mkdir -p "$DRILL/system"
cp "$FIX/system/frag_a.xyz" "$FIX/system/frag_b.xyz" "$DRILL/system/"
cat > "$DRILL/input_qct.txt" <<'EOF'
# auto_sample drill case - Tier 1 auto wiring
SYSTEM_DIR = system
TITLE = auto sample drill - origin-spring well
E_REL = 50.0
B_MAX = 1.0
R_SEP = 6.0
DT = 0.05
DIST_SCHEME = THERMALIZE,ROTATION
PROJECTILE = frag_a
TARGET = frag_b
T_EQ = 100.0
N_EQ = 5
DT_EQ = 0.05
J_ROT = 3
N_B = 1
N_TRAJ = 2
MAX_STEPS = 50
ISEED = 20261012
GWRITE_LEVEL = 2
# ---- auto_sample container keys ----
K_SPRING = 0.05
EOF
cd "$DRILL"
rm -f auto_arm_witness.txt
if "$BLD/venus.e" >stdout_a.txt 2>stderr_a.txt; then ok "adiabatic run (default ELEC_METHOD)"
else bad "adiabatic run"; cat stderr_a.txt; fi
[ ! -f auto_arm_witness.txt ] && ok "arm NOT called under adiabatic default" || bad "arm witness present under adiabatic"
grep -q "outcome counts: bound=2" stdout_a.txt && ok "classify stub ran (bound=2)" || bad "classify outcome line"
# determinism: run again, trajectories byte-identical
cp trajectories.txt /tmp/drill_traj_first.txt
"$BLD/venus.e" >/dev/null 2>&1
cmp -s /tmp/drill_traj_first.txt trajectories.txt && ok "double-run determinism" || bad "determinism"
# R_TERM branch: present key -> immediate termination, still classified
printf 'R_TERM = 2.0\n' >> input_qct.txt
if "$BLD/venus.e" >stdout_b.txt 2>stderr_b.txt; then ok "R_TERM run"
else bad "R_TERM run"; cat stderr_b.txt; fi
grep -q "outcome counts: bound=2" stdout_b.txt && ok "R_TERM branch classified" || bad "R_TERM outcome"
# method branch: select the package -> arm witness appears
sed -i '/^K_SPRING = 0.05/a ELEC_METHOD = auto_tdhf' input_qct.txt
rm -f auto_arm_witness.txt
if "$BLD/venus.e" >stdout_c.txt 2>stderr_c.txt; then ok "auto_tdhf run"
else bad "auto_tdhf run"; cat stderr_c.txt; fi
[ -f auto_arm_witness.txt ] && ok "arm witness written under ELEC_METHOD = auto_tdhf" || bad "arm witness missing"

# ---- leg 5: zero-touch close ----
cd "$SRC"; restore_tree; trap - EXIT
make >/dev/null 2>&1 || { bad "restored rebuild"; exit 1; }
dirty=$(git -C "$ROOT" status --porcelain | grep -v 'venus2026/' | grep -v '^??')
if [ -z "$dirty" ]; then ok "zero-touch (tracked tree clean, pre-existing deletions excluded)"
else bad "zero-touch"; echo "$dirty" | head; fi
echo "drill summary: pass=$pass fail=$fail"
[ "$fail" -eq 0 ]
