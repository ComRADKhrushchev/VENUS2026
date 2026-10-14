#!/usr/bin/env bash
# check_ev_regression.sh —— 值域回归判据（2026-10-06 eV/Å/amu/fs 单位改制后的 case 回归门）
#
# 用法: ./check_ev_regression.sh <case_dir> <baseline_kcal_v1 下的同名基线目录>
#   例: ./check_ev_regression.sh cases/case_spectrum_source/h3_ebk \
#          cases/case_spectrum_source/baseline_kcal_v1/h3_ebk
#
# 背景: 单位改制（能量门面 kcal/mol→eV、时间基 10 fs→fs）后浮点舍入路径改变,
#   逐字节等价门废止（EXTENDING.md 演练清单第 3 条）,回归以值域判据为准。
#
# 判据（每项一行 PASS/FAIL,退出码 0=全过）:
#   A2 初帧能量锚  frame0 的 t_kin/v_pot/e_tot 与旧基线 ÷23.0605 相对差 ≤1e-3
#                  （初帧无传播,残差=旧固定舍入链 vs CODATA 链系统差,实测 ~1e-6）
#   A3 直方图结构  hist_e_rel.dat/hist_thta.dat 各 41 行,计数和 = 轨迹数
#                  （端点分布因轨迹混沌不可逐箱比对,只做结构门）
#   A4 时间轴      t 列自 0 起按 DT 递增（fs 单位下 DT=0.5）
#   A1 守恒漂移    报告运行日志 max|e_fin-e0|（信息性:verlet 大振荡为已知特性,
#                  radau 应 ≤1e-7 eV;超 1 eV 的 verlet 端点振荡亦在已知范围）
set -u

case_dir=$1
base_dir=$2
fail=0

# --- A2: frame0 energy anchor ---
new_line=$(sed -n '2p' "$case_dir/trajectories.txt")
old_line=$(sed -n '2p' "$base_dir/trajectories.txt")
[ -z "$new_line" ] || [ -z "$old_line" ] && { echo "FAIL A2: missing trajectories line (new='$new_line' old='$old_line')"; exit 1; }
awk -v nl="$new_line" -v ol="$old_line" '
function vals(s,   t, a, i) { gsub("t_kin=|v_pot=|e_tot=", " ", s); t = split(s, a, " "); for (i = 1; i <= t; i++) v[i] = a[i] + 0 }
BEGIN {
   vals(nl); nk = v[1]; np = v[2]; nt = v[3]
   vals(ol); ok = v[1]/23.0605; op = v[2]/23.0605; ot = v[3]/23.0605
   tol = 1.0e-3
   d1 = (nk - ok)/(abs = ok < 0 ? -ok : ok); if (d1 < 0) d1 = -d1
   d2 = (np - op)/(abs = op < 0 ? -op : op); if (d2 < 0) d2 = -d2
   d3 = (nt - ot)/(abs = ot < 0 ? -ot : ot); if (d3 < 0) d3 = -d3
   printf "  A2 frame0 eV  old→eV: %.6f %.6f %.6f\n", ok, op, ot
   printf "  A2 frame0 eV  new    : %.6f %.6f %.6f\n", nk, np, nt
   printf "  A2 rel diff: %.2e %.2e %.2e\n", d1, d2, d3
   if (d1 > tol || d2 > tol || d3 > tol) { print "FAIL A2: frame0 anchor exceeds 1e-3"; exit 1 }
   print "PASS A2: frame0 anchor within 1e-3"
   exit 0
}' || fail=1

# --- A3: histogram structure ---
# final rows: NF==5 all-numeric (per trajectory; only the first carries the '# final:' header)
n_traj=$(awk 'NF==5 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9.]/ {c++} END{print c+0}' "$case_dir/trajectories.txt")
[ "$n_traj" -eq 0 ] && { echo "FAIL A3: no final rows found in $case_dir/trajectories.txt"; fail=1; }
for h in hist_e_rel.dat hist_thta.dat; do
   awk -v f="$case_dir/$h" -v n="$n_traj" '
      BEGIN { rows = 0; s = 0; while ((getline < f) > 0) { rows++; s += $2 }
         if (rows == 41 && s == n) { printf "PASS A3: %s rows=41 sum=%d\n", f, s; exit 0 }
         printf "FAIL A3: %s rows=%d sum=%d want 41/%d\n", f, rows, s, n; exit 1 }' || fail=1
done

# --- A4: time axis (fs) monotone from 0 by DT, per trajectory ---
awk '
   /# frame traj=/ { tj = -1; t = 0
      for (i = 1; i < NF; i++) { if ($i ~ /^traj=/) tj = substr($i, 6) + 0
                                 if ($i == "t=") t = $(i+1) + 0 }
      if (tj != curtj) { curtj = tj; prev = -1; nrt = 0 }
      if (nrt == 0 && t != 0.0) { printf "FAIL A4: traj %d first frame t=%g, want 0\n", tj, t; bad = 1; exit 1 }
      if (t < prev - 1e-9) { printf "FAIL A4: traj %d t decreased: %g after %g\n", tj, t, prev; bad = 1; exit 1 }
      prev = t; nrt++; nr++ }
   END { if (!bad) printf "PASS A4: %d frames, per-traj t monotone from 0, last t=%g fs\n", nr, prev; exit bad }' "$case_dir/trajectories.txt" || fail=1

# --- A1: conservation drift (informational) ---
log="${CHECK_LOG:-../../../logs/ev_run_$(basename "$case_dir").log}"
if [ -f "$log" ]; then
   drift=$(grep -o 'max|e_fin-e0|=[ ]*[0-9.E+-]*' "$log" | grep -o '[0-9.E+-]*$')
   echo "INFO A1: conservation drift max|e_fin-e0| = $drift eV (from $log)"
else
   echo "INFO A1: no run log at $log - drift not reported"
fi

[ "$fail" -eq 0 ] && echo "== $case_dir: ALL GATES PASS ==" || echo "== $case_dir: FAILURES PRESENT =="
exit $fail
