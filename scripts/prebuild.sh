#!/usr/bin/env bash
# prebuild.sh <srcdir> <builddir> - the two build-generated artifacts (K3, E11 drag-in):
#   <builddir>/registry_gen.f90 : the aggregate registration module - every reg_*.f90
#       under methods/ and systems/ enters here. NAMING CONTRACT: a methods/ reg
#       module's entry point is <module_name>_all (the E11 reg_<system>_all rule;
#       reg_adiabatic -> reg_adiabatic_all) and it MAY also host <module_name>_init
#       (the member-parameter load - buffered input keys -> the member's own init).
#       A systems/ (container) reg module's entry point is <module_name>_slots - the
#       single handover that fills the schedule slots of interface/container_sched,
#       which then owns the run order and the initialization-assembly guard.
#       A fresh make picks a dragged-in folder up automatically.
#   <builddir>/depend.mk        : object-level use-dependencies for make - depending on
#       the DEFINING unit's .o (not the .mod) sidesteps make's two-product problem and
#       is parallel-safe (the .mod is written before the .o completes).
# Idempotent: reruns rewrite identical content.
set -euo pipefail
SRC="$1"; BLD="$2"
mkdir -p "$BLD"

# sync-protect: emit through <name>.tmp; overwrite the target ONLY when content differs
# (stable timestamps when unchanged - the Makefile runs this every make, phony)
sync_out() { # $1 target, $2 tmp
   if [ ! -f "$1" ] || ! cmp -s "$1" "$2"; then mv -f "$2" "$1"; else rm -f "$2"; fi
}

# ---------- gather sources ----------
tree_files=$(cd "$SRC" && ls -1 *.f90 */*.f90 2>/dev/null || true)
sys_files=$(cd "$SRC/.." && ls -1 systems/*/*.f90 2>/dev/null || true)

# ---------- module -> defining source map ----------
declare -A mod_owner
scan_modules() { # $1 = source path (relative for printing), reads file content from $2
   local rel="$1" f="$2" m
   for m in $(grep -E '^module[[:space:]]+[a-zA-Z_]' "$f" | tr -d '' | awk '{print $2}' || true); do
      mod_owner["$m"]="$rel"
   done
}
while IFS= read -r f; do
   [ -z "$f" ] && continue
   scan_modules "$f" "$SRC/$f"
done <<< "$tree_files"
while IFS= read -r f; do
   [ -z "$f" ] && continue
   scan_modules "../$f" "$SRC/../$f"
done <<< "$sys_files"

# ---------- registry_gen.f90 ----------
# methods/ reg modules: the <module_name>_all registrants (methods only - a container
#   registers nothing itself; its slots go through the schedule)
mods=""
for f in $tree_files; do
   case "$f" in methods/reg_*.f90) mods="$mods $(basename "${f%.f90}")";; esac
done

# a methods/ reg module MAY also host <module_name>_init - the member-parameter load
#   entry (buffered input keys -> the member's own init; the reg file that owns the
#   member's registration owns its initialization assembly). Probe each reg source for
#   the routine; the probe only builds the aggregate, never gates the drag-in.
probe_init() { # $1 = module name, $2 = source content path; prints the name when hosted
   if grep -Eq "^[[:space:]]*subroutine[[:space:]]+${1}_init\(" "$2"; then
      printf '%s\n' "$1"
   fi
   return 0
}
imods=""
for f in $tree_files; do
   case "$f" in methods/reg_*.f90) imods="$imods $(probe_init "$(basename "${f%.f90}")" "$SRC/$f")";; esac
done
has_init() { case " $imods " in *" $1 "*) return 0;; *) return 1;; esac }

# systems/ (container) reg modules: each hosts <module_name>_slots - the single
#   handover filling the schedule slots (interface/container_sched). The handover is
#   probed per source; the schedule then runs the slots at their lifecycle points.
probe_slots() { # $1 = module name, $2 = source content path; prints the name when hosted
   if grep -Eq "^[[:space:]]*subroutine[[:space:]]+${1}_slots\(" "$2"; then
      printf '%s\n' "$1"
   fi
   return 0
}
zmods=""
for f in $sys_files; do
   case "$f" in */reg_*.f90) zmods="$zmods $(probe_slots "$(basename "${f%.f90}")" "$SRC/../$f")";; esac
done

# ---------- auto-PES container scan (2026-10-12, docs/plans/2026-10-12-auto-pes-container.md) ----------
# A systems/* folder WITHOUT reg_*.f90 is a pure-PES container (Tier 1): the
# ten reserved names (pes/keys/load/term/classify + prop/match/sync/init/method_name)
# its modules expose get woven into registry_gen below. A folder WITH reg_*.f90
# is untouched (the explicit reg is authoritative - Tier 2).
auto_fatal() { # $1 = message (the design's named aborts, prebuild side)
   echo "prebuild.sh: auto-PES error: $1" >&2
   exit 1
}

# candidate folder: one level under systems/, holding .f90 (the reg_*.f90
# container form is retired 2026-10-14 - a leftover reg file is just another
# source file to the wildcard; folder identity is the reserved names alone)
auto_dir=""
for d in "$SRC/../systems"/*/; do
   [ -d "$d" ] || continue
   ls "$d"*.f90 >/dev/null 2>&1 || continue
   if [ -n "$auto_dir" ]; then
      auto_fatal "two auto containers are installed ('$auto_dir' and '$d') - one container per process"
   fi
   auto_dir="$d"
done

# reserved-name scan: emits "kind name module" per definition (kind s/f).
# Comment tails are stripped and the keyword may sit behind a type/purity
# prefix - "logical function term(sta)" is the idiomatic declaration form
# (reference_box writes it; a line-start anchor would silently MISS it).
# "end subroutine pes" never matches (the end-prefix exclusion). The name's
# argument paren is stripped, not demanded: a paren-less reserved-name definition
# still counts and fatals as a duplicate/wrong-form - loud, never a wrong weave.
# A wrong-kind hit (subroutine keys / function pes) emits its actual kind -
# the caller's form check fatals on it.
auto_scan_file() { # $1 = source path
   awk '
      BEGIN { split("pes load prop match sync init observables", sa); for (i in sa) sk[sa[i]] = 1
              split("keys term classify method_name", fa); for (i in fa) fk[fa[i]] = 1 }
      {
         line = $0
         sub(/!.*/, "", line)                       # strip the comment tail
         sub(/^[[:space:]]+/, "", line)             # normalize the lead
         lline = tolower(line)                      # the matching copy (keywords are case-insensitive)
         if (lline ~ /^module[[:space:]]+[a-z_]/) {
            split(line, mw, /[[:space:]]+/)         # original case for the emitted use line
            if (mw[2] != "procedure") mod = mw[2]
         }
         if (match(lline, /(subroutine|function)[[:space:]]+[a-z_]/)) {
            if (RSTART > 1 && substr(lline, 1, RSTART - 1) ~ /end[[:space:]]*$/) next
            split(substr(lline, RSTART), w, /[[:space:]]+/)
            kw = w[1]; nm = w[2]; sub(/\(.*/, "", nm)
            if (nm in sk || nm in fk) print (kw == "subroutine" ? "s" : "f"), nm, mod
         }
      }' "$1"
}

declare -A auto_mod                     # reserved name -> owning module
auto_pair=no
if [ -n "$auto_dir" ]; then
   defs=$(for f in "$auto_dir"*.f90; do auto_scan_file "$f"; done)
   for nm in pes keys load term classify prop match sync init method_name observables; do
      case "$nm" in pes|load|prop|match|sync|init|observables) want=s;; *) want=f;; esac
      got=$(printf '%s\n' "$defs" | awk -v n="$nm" '$2 == n {print}')
      cnt=$(printf '%s\n' "$got" | grep -c .) || true
      [ "$cnt" -eq 0 ] && continue
      [ "$cnt" -ge 2 ] && auto_fatal "reserved name '$nm' is defined $cnt times in $auto_dir (exactly one allowed)"
      read -r kind _ mod _ <<<"$(printf '%s\n' "$got")"
      [ "$kind" != "$want" ] && auto_fatal "reserved name '$nm' must be a $([ "$want" = s ] && echo subroutine || echo function), the other form was found in $auto_dir"
      [ -z "$mod" ] && auto_fatal "reserved name '$nm' sits outside any module in $auto_dir (the weave needs a module owner)"
      auto_mod[$nm]=$mod
   done
   n_pm=0
   [ -n "${auto_mod[prop]:-}" ] && n_pm=$((n_pm+1))
   [ -n "${auto_mod[match]:-}" ] && n_pm=$((n_pm+1))
   if [ "$n_pm" -eq 1 ]; then
      auto_fatal "prop and match must be defined together in $auto_dir (a half package is not a method)"
   fi
   [ "$n_pm" -eq 2 ] && auto_pair=yes
   if [ "$auto_pair" = yes ]; then
      [ -z "${auto_mod[method_name]:-}" ] && auto_fatal "prop+match are defined but method_name is missing in $auto_dir (the registration name is author data)"
   else
      [ -n "${auto_mod[sync]:-}" ] && auto_fatal "sync is defined but prop+match are not (a lone sync member is not a package) in $auto_dir"
      [ -n "${auto_mod[init]:-}" ] && auto_fatal "init is defined but prop+match are not (a lone init member is not a package) in $auto_dir"
      [ -n "${auto_mod[method_name]:-}" ] && auto_fatal "method_name is defined but prop+match are not in $auto_dir"
   fi
   [ -z "${auto_mod[pes]:-}" ] && auto_fatal "$auto_dir is installed (folder present) but wires nothing: expose exactly one 'pes' entry (the reserved names)"
fi

{
   echo "! registry_gen.f90 - BUILD-GENERATED by scripts/prebuild.sh (do not edit)"
   echo "! The aggregate registration module (E11/E14/E17-E19): every reg_*.f90 under"
   echo "! methods/ (method packages) and systems/ (the containers) enters here -"
   echo "! dragging a new folder in is the whole installation act."
   echo "! registry_gen_keys() runs BEFORE the parse: each container hands its slots"
   echo "! over (reg_<system>_slots fills the schedule slots), then the schedule runs"
   echo "! the key-vocabulary slot so the grammar admits the container key names."
   echo "! registry_gen_all() is the registration entry the driver calls at assembly:"
   echo "! every methods/ _all first, then the schedule's registration-time slots"
   echo "! (force slot, method rows, registration-time exports)."
   echo "! registry_gen_init() is the initialization-assembly aggregate: the"
   echo "! schedule's seam slots FIRST (the container loads its material - a"
   echo "! load-dependent PES must be computable before any member seam that"
   echo "! probes its force, e.g. the spectrum derivation; the guard lives"
   echo "! once in container_sched_init), then every methods/ _init entry"
   echo "! (unselected members skip inside their own guard)."
   if [ -n "$auto_dir" ]; then
      echo "! An auto container's weave (auto_slots) enters here when a systems/"
      echo "! folder is present - Tier 1, docs/plans/2026-10-12-auto-pes-container.md"
   fi
   echo "module registry_gen"
   echo "   use container_sched, only: container_sched_keys, container_sched_all, &"
   echo "                              container_sched_init"
   if [ -n "$auto_dir" ]; then
      echo "   use container_sched, only: container_sched_bind"
      [ -n "${auto_mod[term]:-}" ] && echo "   use force_interface, only: container_bind_pes, container_bind_term"
      [ -z "${auto_mod[term]:-}" ] && echo "   use force_interface, only: container_bind_pes"
      [ -n "${auto_mod[classify]:-}" ] && echo "   use final_state,     only: container_bind_classify"
      [ -n "${auto_mod[keys]:-}" ] && echo "   use input,           only: input_declare_keys"
      if [ "$auto_pair" = yes ]; then
         echo "   use elec_interface,  only: method_reg"
         [ -n "${auto_mod[init]:-}" ] && echo "   use config,          only: electronic"
      fi
      for m in $(printf '%s\n' "${auto_mod[@]}" | sort -u); do
         names=""; first=1
         for nm in pes keys load term classify prop match sync init method_name observables; do
            if [ "${auto_mod[$nm]:-}" = "$m" ]; then
               if [ "$first" -eq 1 ]; then names="$nm"; first=0; else names="$names, $nm"; fi
            fi
         done
         echo "   use $m, only: $names"
      done
   fi
   for m in $mods; do
      names="${m}_all"
      if has_init "$m"; then names="$names, ${m}_init"; fi
      echo "   use $m, only: $names"
   done
   for z in $zmods; do
      echo "   use $z, only: ${z}_slots"
   done
   echo "   implicit none"
   echo "   private"
   echo "   public :: registry_gen_all"
   echo "   public :: registry_gen_init"
   echo "   public :: registry_gen_keys"
   echo "contains"
   if [ -n "$auto_dir" ]; then
      args=""
      [ -n "${auto_mod[keys]:-}" ]     && args="$args keys"
      [ -n "${auto_mod[load]:-}" ]     && args="$args load"
      args="$args force"
      [ "$auto_pair" = yes ]         && args="$args method"
      [ -n "${auto_mod[term]:-}" ]     && args="$args term"
      [ -n "${auto_mod[init]:-}" ]     && args="$args init"
      [ -n "${auto_mod[classify]:-}" ] && args="$args classify"
      [ -n "${auto_mod[observables]:-}" ] && args="$args columns"
      echo "   subroutine auto_slots()"
      ia=0
      for a in $args; do
         ia=$((ia+1))
         if [ "$ia" -eq 1 ]; then
            printf '      call container_sched_bind(%s = auto_%s' "$a" "$a"
         else
            printf ', &\n                                %s = auto_%s' "$a" "$a"
         fi
      done
      echo ")"
      echo "   end subroutine auto_slots"
      [ -n "${auto_mod[keys]:-}" ] && { echo "   subroutine auto_keys()"
         echo "      call input_declare_keys(keys())"
         echo "   end subroutine auto_keys"; }
      [ -n "${auto_mod[load]:-}" ] && { echo "   subroutine auto_load()"
         echo "      call load()"
         echo "   end subroutine auto_load"; }
      echo "   subroutine auto_force()"
      echo "      call container_bind_pes(pes)"
      echo "   end subroutine auto_force"
      if [ "$auto_pair" = yes ]; then
         echo "   subroutine auto_method()"
         if [ -n "${auto_mod[sync]:-}" ]; then echo "      call method_reg(method_name(), prop, match, sync)"
         else echo "      call method_reg(method_name(), prop, match)"; fi
         echo "   end subroutine auto_method"
      fi
      [ -n "${auto_mod[term]:-}" ] && { echo "   subroutine auto_term()"
         echo "      call container_bind_term(term)"
         echo "   end subroutine auto_term"; }
      [ -n "${auto_mod[init]:-}" ] && { echo "   subroutine auto_init()"
         echo "      if (trim(electronic%method) == method_name()) call init()"
         echo "   end subroutine auto_init"; }
      [ -n "${auto_mod[classify]:-}" ] && { echo "   subroutine auto_classify()"
         echo "      call container_bind_classify(classify)"
         echo "   end subroutine auto_classify"; }
      [ -n "${auto_mod[observables]:-}" ] && { echo "   subroutine auto_columns()"
         echo "      call observables()"
         echo "   end subroutine auto_columns"; }
   fi
   echo "   subroutine registry_gen_keys()"
   [ -n "$auto_dir" ] && echo "      call auto_slots()"
   for z in $zmods; do echo "      call ${z}_slots()"; done
   echo "      call container_sched_keys()"
   echo "   end subroutine registry_gen_keys"
   echo "   subroutine registry_gen_all()"
   for m in $mods; do echo "      call ${m}_all()"; done
   echo "      call container_sched_all()"
   echo "   end subroutine registry_gen_all"
   echo "   subroutine registry_gen_init()"
   echo "      call container_sched_init()"
   for m in $imods; do echo "      call ${m}_init()"; done
   echo "   end subroutine registry_gen_init"
   echo "end module registry_gen"
} > "$BLD/registry_gen.f90.tmp"
sync_out "$BLD/registry_gen.f90" "$BLD/registry_gen.f90.tmp"

# the generated unit defines a module too - it must enter the module->owner map BEFORE
# dependency emission (a driver `use registry_gen` resolves through this entry)
scan_modules registry_gen "$BLD/registry_gen.f90"

# ---------- depend.mk ----------
objpath() { # tree: globals/state.f90 -> build/globals/state.o ; sys: systems/x/y.f90 -> build/sys/x/y.o ; generated: registry_gen
   case "$1" in
      registry_gen)      echo "$BLD/registry_gen.o" ;;
      ../*)              echo "$BLD/sys/$(echo "${1#../systems/}" | sed 's|\.f90$|.o|')" ;;
      *)                 echo "$BLD/$(echo "$1" | sed 's|\.f90$|.o|')" ;;
   esac
}
{
   echo "# depend.mk - BUILD-GENERATED by scripts/prebuild.sh (do not edit)"
   emit_deps() { # $1 = dep-owner source rel path, $2 = file content path
      local uses u d deps=""
      # || true: a use-free unit must not kill the scan
      uses=$(grep -E '^[[:space:]]*use[[:space:]]+[a-zA-Z_]' "$2" | tr -d '' | awk '{print $2}' | sed 's/,.*//' | sort -u || true)
      for u in $uses; do
         d="${mod_owner[$u]:-}"
         [ -n "$d" ] && deps="$deps $(objpath "$d")"
      done
      echo "$(objpath "$1"):$deps"
   }
   while IFS= read -r f; do
      [ -z "$f" ] && continue
      emit_deps "$f" "$SRC/$f"
   done <<< "$tree_files"
   while IFS= read -r f; do
      [ -z "$f" ] && continue
      emit_deps "../$f" "$SRC/../$f"
   done <<< "$sys_files"
   emit_deps registry_gen "$BLD/registry_gen.f90"
} > "$BLD/depend.mk.tmp"
sync_out "$BLD/depend.mk" "$BLD/depend.mk.tmp"
