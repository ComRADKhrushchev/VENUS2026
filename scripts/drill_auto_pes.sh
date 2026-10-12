#!/usr/bin/env bash
# drill_auto_pes.sh - the Tier 1 auto-PES container drill (criterion family):
#   positive legs (clean-tree drag-in first make, woven content, pes-only
#   partial-shape weave, adiabatic default run with the arm witness ABSENT,
#   ELEC_METHOD = auto_tdhf run with the witness PRESENT, R_TERM both
#   branches, double-run determinism) + the runtime negatives (reference +
#   auto coexistence STOP 1; a package naming itself 'adiabatic' aborts at
#   method_reg assembly) + the prebuild negative family (eight named aborts -
#   the committed spec §5 negative coverage) + the zero-touch close.
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
   # the dragged fixture is a disposable working copy - restore the shelf
   # copy from git (any exit form leaves the tracked tree canonical)
   rm -rf "$FIX"
   git -C "$ROOT" checkout -q -- systems/inactive/auto_sample
   [ -d "$REFSHELF" ] && mv "$REFSHELF" "$REF"
   rm -rf "$DRILL"
}
trap restore_tree EXIT

drag() { [ -d "$FIX" ] || mv "$SHELF" "$FIX"; }   # idempotent drag-in

cd "$SRC"

# ---- leg 1: drag in, FIRST make green (clean tree - the bootstrap path) ----
make clean >/tmp/drill_clean.log 2>&1 || { bad "make clean (leg-1 precondition)"; exit 1; }
mv "$SHELF" "$FIX"
if make >/tmp/drill_make.log 2>&1; then ok "first make green (clean tree, drag-in, no second make)"
else bad "first make (precondition: clean tree — the drill now makes clean itself; a failure here is a real defect)"; tail -20 /tmp/drill_make.log; exit 1; fi

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
call container_sched_bind(keys = auto_keys, &
force = auto_force, &
PATS

# ---- leg 3: partial-shape weave (a pes-only container binds force to the
#      FORCE slot - regression for the keyword-form fix: with positional
#      actuals a non-prefix name set mis-bound one slot left, the force slot
#      stayed null) ----
mkdir -p "$DRILL"
mv "$FIX/io.f90" "$DRILL/park_io.f90"
mv "$FIX/tdhf_stub.f90" "$DRILL/park_tdhf.f90"
cp "$FIX/pot.f90" "$DRILL/park_pot.f90"
cat > "$FIX/pot.f90" <<'EOF'
module auto_sample_pot
   implicit none
   private
   public :: pes
contains
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      v = 0.0d0
      g = 0.0d0
   end subroutine pes
end module auto_sample_pot
EOF
(cd "$SRC" && bash ../scripts/prebuild.sh . build/gfortran) >"$DRILL/partial.log" 2>&1
nb=$(grep -c ' = auto_' "$BLD/registry_gen.f90" || true)
if grep -qF 'call container_sched_bind(force = auto_force)' "$BLD/registry_gen.f90" \
   && [ "$nb" -eq 1 ]; then
   ok "pes-only weave: force binds the force slot, sole bind actual (keyword form)"
else
   bad "pes-only weave (bind actuals=$nb, want 1 with force named)"
fi
mv "$DRILL/park_io.f90" "$FIX/io.f90"
mv "$DRILL/park_tdhf.f90" "$FIX/tdhf_stub.f90"
mv "$DRILL/park_pot.f90" "$FIX/pot.f90"
(cd "$SRC" && bash ../scripts/prebuild.sh . build/gfortran) >/dev/null 2>&1

# ---- leg 4: runtime coexistence abort (reference + auto = STOP 1 at assembly;
#      registry_gen_keys runs before the input read, so no input file needed) ----
"$BLD/venus.e" >/dev/null 2>/tmp/drill_coexist.err
if grep -q "slot is already filled" /tmp/drill_coexist.err; then ok "coexistence STOP 1"
else bad "coexistence STOP 1"; cat /tmp/drill_coexist.err; fi

# ---- leg 5: shelve reference, rebuild, scratch-case runs ----
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
# method-name collision: the package registering a name already in the table
# (the built-in 'adiabatic') aborts at assembly (method_reg duplicate)
sed -i "s/s = 'auto_tdhf'/s = 'adiabatic'/" "$FIX/tdhf_stub.f90"
if (cd "$SRC" && make) >"$DRILL/collide_make.log" 2>&1; then :; else bad "collision rebuild"; exit 1; fi
if "$BLD/venus.e" >stdout_d.txt 2>stderr_d.txt; then
   bad "adiabatic collision run (expected method_reg abort)"
else
   grep -q 'duplicate method name "adiabatic"' stderr_d.txt \
      && ok "collision abort: duplicate method name \"adiabatic\"" \
      || { bad "collision abort text"; cat stderr_d.txt; }
fi
rm -rf "$FIX"
git -C "$ROOT" checkout -q -- systems/inactive/auto_sample
drag
(cd "$SRC" && make) >/dev/null 2>&1 || { bad "post-collision canonical rebuild"; exit 1; }

# ---- leg 6: prebuild negative family (mutate -> named abort -> restore).
#      The dragged fixture is an untracked working copy; each leg's restore is
#      rm + git checkout of the shelf + re-drag (the group ends PARKED). ----
neg() { # $1 = id, $2 = expected abort fragment (fixed string); mutation applied by caller
   (cd "$SRC" && bash ../scripts/prebuild.sh . build/gfortran) >"$DRILL/neg$1.out" 2>&1
   rc=$?
   if [ "$rc" -ne 0 ] && grep -qF "$2" "$DRILL/neg$1.out"; then
      ok "prebuild negative $1: $2"
   else
      bad "prebuild negative $1 (rc=$rc, want '$2')"; sed -n '1,4p' "$DRILL/neg$1.out"
   fi
   rm -rf "$FIX"
   git -C "$ROOT" checkout -q -- systems/inactive/auto_sample
}

# 6.1: no pes anywhere (pot parked, io/met carry no force entry)
drag
mv "$FIX/pot.f90" "$DRILL/neg1_pot.f90"
neg 1 "wires nothing: expose exactly one 'pes'"

# 6.2: a second module defining pes (duplicate definition)
drag
cat >> "$FIX/pot.f90" <<'EOT'
module dup_mod
contains
subroutine pes(q, mass, g, v)
real(8), intent(in) :: q(:), mass(:)
real(8), intent(out) :: g(:), v
end subroutine
end module
EOT
neg 2 "'pes' is defined 2 times"

# 6.3: keys declared as a subroutine (wrong form)
drag
sed -i 's/^   function keys() result(kk)/   subroutine keys() result(kk)/' "$FIX/pot.f90"
neg 3 "must be a function, the other form"

# 6.4: a second candidate folder (two auto containers)
drag
cp -r "$FIX" "$ROOT/systems/auto_two"
neg 4 "two auto containers are installed"
rm -rf "$ROOT/systems/auto_two"

# 6.5: match renamed away (half a method package)
drag
sed -i 's/\bmatch\b/match2/g' "$FIX/tdhf_stub.f90"
neg 5 "prop and match must be defined together"

# 6.6: method_name renamed away (package without its registration name)
drag
sed -i 's/\bmethod_name\b/method_label/g' "$FIX/tdhf_stub.f90"
neg 6 "method_name is missing"

# 6.7: prop+match commented out, sync kept (lone sync member)
drag
sed -i -e 's/^   subroutine prop(sta, dt)$/! &/' -e 's/^   subroutine match(sta)$/! &/' "$FIX/tdhf_stub.f90"
neg 7 "sync is defined but prop+match are not"

# 6.8: only method_name kept (prop/match/sync/arm commented out)
drag
sed -i -e 's/^   subroutine prop(sta, dt)$/! &/' -e 's/^   subroutine match(sta)$/! &/' \
       -e 's/^   subroutine sync(s)$/! &/' -e 's/^   subroutine arm()$/! &/' "$FIX/tdhf_stub.f90"
neg 8 "method_name is defined but prop+match are not"

# fixture PARKED (neg 8 left it there); restore reference for the close
mv "$REFSHELF" "$REF"

# ---- leg 7: zero-touch close ----
cd "$SRC"; restore_tree; trap - EXIT
make >/dev/null 2>&1 || { bad "restored rebuild"; exit 1; }
dirty=$(git -C "$ROOT" status --porcelain | grep -v 'venus2026/' | grep -v '^??')
if [ -z "$dirty" ]; then ok "zero-touch (tracked tree clean, pre-existing deletions excluded)"
else bad "zero-touch"; echo "$dirty" | head; fi
echo "drill summary: pass=$pass fail=$fail"
[ "$fail" -eq 0 ]
