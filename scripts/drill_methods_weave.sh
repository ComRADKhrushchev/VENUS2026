#!/usr/bin/env bash
# drill_methods_weave.sh - the methods fixed-name weave drill (criterion family):
#   wv-1: a fixture member file (methods/samp_probe.f90 exposing member_register)
#         is woven on make - registry_gen carries the wrapper and the full build
#         is green; removing the fixture removes the wrapper (idempotence);
#   wv-2: a prefix file WITHOUT member_register (methods/samp_broken.f90) aborts
#         prebuild with the named contract message (the loud contract - added by
#         the contract task; this leg is RED until then, by TDD design).
# Leaves the tree as found (trap removes fixtures + rebuilds); nonzero on failure.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/src"; BLD="$SRC/build/gfortran"
pass=0; fail=0
ok()  { echo "PASS: $1"; pass=$((pass+1)); }
bad() { echo "FAIL: $1"; fail=$((fail+1)); }

restore_tree() {
   rm -f "$SRC/methods/samp_probe.f90" "$SRC/methods/samp_broken.f90"
   make -C "$SRC" >/dev/null 2>&1
}
trap restore_tree EXIT
cd "$SRC"

# ---- wv-1: fixture member woven ----
cat > methods/samp_probe.f90 <<'EOF'
!=====================================================================
! samp_probe.f90 - WEAVE DRILL FIXTURE (drill_methods_weave.sh only):
!   minimal member exposing the fixed-name entry; registered under a word
!   no input admits (the drill asserts the WEAVE, not the dispatch)
!=====================================================================
module samp_probe
   use sampler, only: samp_reg
   use state,    only: state_t
   implicit none
   private
   public :: member_register
contains
   subroutine probe_draw()
   end subroutine probe_draw

   subroutine probe_realize(sta)
      type(state_t), intent(inout) :: sta
   end subroutine probe_realize

   subroutine member_register()
      call samp_reg('probe', probe_draw, probe_realize)
   end subroutine member_register
end module samp_probe
EOF
if make >/tmp/wv_make1.log 2>&1; then ok "wv-1 make with fixture"; else bad "wv-1 make with fixture"; tail -20 /tmp/wv_make1.log; fi
for pat in 'subroutine auto_samp_probe()' 'use samp_probe, only: member_register' 'call auto_samp_probe()'; do
   grep -qF "$pat" "$BLD/registry_gen.f90" && ok "wv-1 woven: $pat" || bad "wv-1 woven: $pat"
done
rm -f methods/samp_probe.f90
if make >/tmp/wv_make1b.log 2>&1; then ok "wv-1 restore make"; else bad "wv-1 restore make"; fi
if grep -qF 'auto_samp_probe' "$BLD/registry_gen.f90"; then bad "wv-1 wrapper gone after fixture removal"; else ok "wv-1 wrapper gone after fixture removal"; fi

# ---- wv-2: prefix file without member_register aborts prebuild ----
printf 'module samp_broken\n   implicit none\n   private\nend module samp_broken\n' > methods/samp_broken.f90
if make >/tmp/wv_make2.log 2>&1; then
   bad "wv-2 contract fatal (make succeeded, expected the prebuild abort)"
else
   grep -qF 'carries no member_register' /tmp/wv_make2.log && ok "wv-2 contract fatal message" || { bad "wv-2 fatal message text"; tail -5 /tmp/wv_make2.log; }
fi
rm -f methods/samp_broken.f90

echo "drill_methods_weave: pass=$pass fail=$fail"
[ "$fail" -eq 0 ]
