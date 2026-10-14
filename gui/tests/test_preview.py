# -*- coding: utf-8 -*-
"""preview 测试:p_site/u_inc/速度劈分对 beam_laws+inc_* 公式的手算对照。"""
import math
from pathlib import Path

import pytest
from conftest import ROOT

from gui.core import preview, sysdir
from gui.core.inputwrite import parse_input_file


def test_pair_hand_calc():
    """case_harmonic_scatter 参数(法向入射):手算 p_site/v/劈分。"""
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")
    assert sd.ok
    values = dict(parse_input_file(ROOT / "cases" / "case_harmonic_scatter" / "input_qct.txt"))
    scene = preview.assemble_preview(sd, values, i_proj=0, i_targ=1)
    assert scene.beam is not None
    b = scene.beam
    # θ=THTA_MAX 未写 → 0(法向);χ=0;N_B=1 → 代表值 b=0.7·2.0、φ=45°
    assert b.theta == pytest.approx(0.0)
    assert b.b == pytest.approx(0.7 * 2.0)
    assert b.phi == pytest.approx(math.pi / 4)
    assert b.u_inc == pytest.approx([0.0, 0.0, -1.0])
    # 采样盘(N_B=1)显示恒对心(2026-10-14 细化裁决):p_site 取 b=0 值 = (0,0,R_SEP);
    # b 横向项只进 BeamView.b 信息面(固定 b 模式才会位移, 见 dedicated 测试)
    assert b.p_site == pytest.approx([0.0, 0.0, 8.0])
    # μ = 1.008(H2 折合质量 = 2.016²/4.032);E_REL 取自案例输入(8.67283884 eV
    # = 200 kcal/mol ÷ 23.0605, 2026-10-13 单位改基后的换算值);v = √(2·E_REL/(μ·dyn_ev))
    assert b.mu == pytest.approx(1.008)
    assert b.speed == pytest.approx(
        math.sqrt(2 * float(values["E_REL"]) / (1.008 * preview.DYN_EV)))
    # pair 劈分 w1=w2=0.5
    assert b.v_proj == pytest.approx([0.0, 0.0, -0.5 * b.speed])
    assert b.v_targ == pytest.approx([0.0, 0.0, 0.5 * b.speed])
    # 靶 COM 在原点,弹丸 COM 在 p_site
    assert scene.frag_coms[1] == pytest.approx([0.0, 0.0, 0.0])
    assert scene.frag_coms[0] == pytest.approx(b.p_site)
    # 内部几何保持(键长 0.7)
    p = scene.positions
    d = math.dist(p[0], p[1])
    assert d == pytest.approx(0.7)
    assert scene.b_max == pytest.approx(2.0)


def test_surface_fixed_site(tmp_path: Path):
    d = tmp_path / "system"
    d.mkdir()
    (d / "proj.poscar").write_text(
        "one H in a box\n1.0\n9.0 0.0 0.0\n0.0 9.0 0.0\n0.0 0.0 12.0\nH\n1\nCartesian\n4.5 4.5 1.0\n",
        encoding="utf-8",
    )
    sd = sysdir.scan_system_dir(d)
    assert sd.ok and sd.has_cell
    values = {
        "SCATTERING": "SURFACE", "CELL": "9.0 0.0 0.0, 0.0 9.0 0.0, 0.0 0.0 12.0",
        "E_REL": "100.0", "R_SEP": "6.0", "N_AIM": "2",
        "AIM_X": "1.5", "AIM_Y": "2.5", "N_B": "0",
    }
    scene = preview.assemble_preview(sd, values, i_proj=0)
    b = scene.beam
    assert b is not None
    # 定点落点 + b=B_MAX 未写 → 0(表面法定不写)→ p_site = (1.5, 2.5, 6.0)
    assert b.p_site == pytest.approx([1.5, 2.5, 6.0])
    assert b.mu == pytest.approx(1.008)
    assert b.speed == pytest.approx(math.sqrt(2 * 100.0 / (1.008 * preview.DYN_EV)))
    assert b.v_targ is None  # 表面:全额赋予弹丸
    assert b.v_proj == pytest.approx([0.0, 0.0, -b.speed])


def test_oblique_incidence(tmp_path: Path):
    """斜入射(θ=30°, χ=90°):u_inc 三分量手算。"""
    d = tmp_path / "system"
    d.mkdir()
    (d / "a.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    (d / "b.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    sd = sysdir.scan_system_dir(d)
    th = math.pi / 6
    values = {
        "E_REL": "10.0", "R_SEP": "5.0", "B_MAX": "1.0", "N_B": "0",
        "THTA_MAX": str(th), "CHI": str(math.pi / 2),
    }
    scene = preview.assemble_preview(sd, values, i_proj=0, i_targ=1)
    u = scene.beam.u_inc
    assert u == pytest.approx([-math.sin(th) * math.cos(math.pi / 2), -math.sin(th), -math.cos(th)])
    assert math.sqrt(sum(c * c for c in u)) == pytest.approx(1.0)


def test_single_fragment_no_beam():
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_b_layer" / "andersen" / "system")
    assert sd.ok
    values = dict(parse_input_file(ROOT / "cases" / "case_b_layer" / "andersen" / "input_qct.txt"))
    scene = preview.assemble_preview(sd, values)
    assert scene.beam is None
    assert any("单片段" in n for n in scene.notes)


def test_dyn_ev_mirror():
    """preview.DYN_EV 须等于 consts.f90 的 dyn_ev(单位制常数镜像,防漂移)。"""
    import re
    text = (ROOT / "src" / "globals" / "consts.f90").read_text(encoding="utf-8")
    m = re.search(r"dyn_ev\s*=\s*([0-9.]+)d0", text)
    assert m, "consts.f90 未找到 dyn_ev 定义"
    assert preview.DYN_EV == float(m.group(1))


# ---------- Task 2(2026-10-13 b 落点修正):角色字段 + 环心/取景中心派生 ----------
# b 定义绕靶 COM(inc_pair.f90:110-114:p_site = 落点 + b 偏移 + R_SEP 退距,
# 靶 COM 移到落点锚)——b 环与相机取景中心都应绑靶,视图层经下述纯派生消费。


def _pair_scene(i_proj: int, i_targ: int) -> preview.PreviewScene:
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")
    values = dict(parse_input_file(ROOT / "cases" / "case_harmonic_scatter" / "input_qct.txt"))
    return preview.assemble_preview(sd, values, i_proj=i_proj, i_targ=i_targ)


def test_scene_role_fields_pair_both_orders():
    """pair 两序:i_proj/i_targ 如实回填;靶 COM=落点锚(原点),弹丸在 p_site。"""
    for ip, it in ((0, 1), (1, 0)):
        sc = _pair_scene(ip, it)
        assert (sc.i_proj, sc.i_targ) == (ip, it)
        assert sc.beam is not None
        assert sc.frag_coms[it] == pytest.approx([0.0, 0.0, 0.0])
        assert sc.frag_coms[ip] == pytest.approx(sc.beam.p_site)


def test_scene_role_fields_surface_unresolved_single(tmp_path: Path):
    """表面(SCATTERING=SURFACE):弹丸在场但无靶片段角色(i_targ=-1);未解析/单片段:双 -1。"""
    d = tmp_path / "system"
    d.mkdir()
    (d / "proj.poscar").write_text(
        "one H in a box\n1.0\n9.0 0.0 0.0\n0.0 9.0 0.0\n0.0 0.0 12.0\nH\n1\nCartesian\n4.5 4.5 1.0\n",
        encoding="utf-8",
    )
    sd = sysdir.scan_system_dir(d)
    values = {"SCATTERING": "SURFACE", "CELL": "9 0 0, 0 9 0, 0 0 12",
              "E_REL": "100.0", "R_SEP": "6.0", "N_B": "0"}
    sc = preview.assemble_preview(sd, values, i_proj=0)
    assert (sc.i_proj, sc.i_targ) == (0, -1)

    sc2 = _pair_scene(-1, -1)
    assert (sc2.i_proj, sc2.i_targ) == (-1, -1) and sc2.beam is None

    sd3 = sysdir.scan_system_dir(ROOT / "cases" / "case_b_layer" / "andersen" / "system")
    vals3 = dict(parse_input_file(ROOT / "cases" / "case_b_layer" / "andersen" / "input_qct.txt"))
    sc3 = preview.assemble_preview(sd3, vals3, i_proj=0)
    assert (sc3.i_proj, sc3.i_targ) == (-1, -1) and sc3.beam is None


def test_ring_center_binds_target_com():
    """b 环心 = 靶 COM(两序同律);角色未解析/表面/单片段 → None(不画环)。"""
    for ip, it in ((0, 1), (1, 0)):
        sc = _pair_scene(ip, it)
        rc = preview.ring_center(sc)
        assert rc is not None, (ip, it)
        assert math.dist(rc, sc.frag_coms[it]) == pytest.approx(0.0, abs=1e-12)
        if it == 0:  # 旧缺陷:硬编码 frag_coms[1] 会把环套在弹丸上
            assert math.dist(rc, sc.frag_coms[1]) > 0.5
    assert preview.ring_center(_pair_scene(-1, -1)) is None


def test_view_center_policy_pair_and_fallback():
    """取景中心:pair=靶 COM(两序);未解析/空场景回退全原子均值。"""
    for ip, it in ((0, 1), (1, 0)):
        sc = _pair_scene(ip, it)
        assert preview.view_center(sc) == pytest.approx(sc.frag_coms[it])
    sc = _pair_scene(-1, -1)
    mean = [sum(p[c] for p in sc.positions) / len(sc.positions) for c in range(3)]
    assert preview.view_center(sc) == pytest.approx(mean)
    assert preview.view_center(preview.PreviewScene()) == [0.0, 0.0, 0.0]


def test_view_center_policy_surface_drop(tmp_path: Path):
    """表面(SCATTERING=SURFACE):横向 = 落点(N_AIM=2 定点),高度 = 衬底(非弹丸)原子均值。"""
    d = tmp_path / "system"
    d.mkdir()
    (d / "proj.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    (d / "slab.poscar").write_text(
        "slab\n1.0\n9.0 0.0 0.0\n0.0 9.0 0.0\n0.0 0.0 12.0\nAg\n2\nCartesian\n"
        "0.5 0.5 0.5\n1.5 1.5 1.5\n",
        encoding="utf-8",
    )
    sd = sysdir.scan_system_dir(d)
    assert sd.ok and sd.has_cell  # POSCAR 带胞为扫描侧事实(信息性;范式由 SCATTERING 键定)
    values = {"SCATTERING": "SURFACE", "CELL": "9 0 0, 0 9 0, 0 0 12",
              "E_REL": "100.0", "R_SEP": "6.0", "N_AIM": "2",
              "AIM_X": "1.5", "AIM_Y": "2.5", "N_B": "0"}
    sc = preview.assemble_preview(sd, values, i_proj=0)  # 名序:proj(0) < slab(1)
    assert (sc.i_proj, sc.i_targ) == (0, -1) and sc.beam is not None
    vc = preview.view_center(sc)
    assert vc[:2] == pytest.approx([1.5, 2.5])  # 横向 = 落点
    a, b = sc.frag_slices[sc.i_proj]  # 高度 = 衬底原子均值(独立复推)
    zs = [p[2] for k, p in enumerate(sc.positions) if not a <= k < b]
    assert vc[2] == pytest.approx(sum(zs) / len(zs))
    assert vc[2] == pytest.approx(1.0)


# ---------- b 显示项(2026-10-13 恒对心批 → 2026-10-14 细化裁决) ----------
# 用户细化: N_B=0(固定 b)→ 弹丸位置随 B_MAX 偏移(物理正确);
# N_B=1(碰撞盘采样)→ 恒对心不位移(b 逐轨迹抽样), 范围由 b 环可视化。


def _pair_scene_bmax(n_b: int, b_max: str) -> preview.PreviewScene:
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")
    values = dict(parse_input_file(ROOT / "cases" / "case_harmonic_scatter" / "input_qct.txt"))
    values.update(B_MAX=b_max, N_B=str(n_b))
    return preview.assemble_preview(sd, values, i_proj=0, i_targ=1)


def test_sampling_placement_invariant_under_bmax():
    """采样盘(N_B=1)恒对心:仅改 B_MAX → 原子摆位/瞄准十字逐点不变。"""
    ref = _pair_scene_bmax(1, "1.0")
    for b_max in ("2.0", "3.5"):
        sc = _pair_scene_bmax(1, b_max)
        assert sc.positions == ref.positions, b_max  # 位级恒等(强于 approx)
        assert sc.aim_marker == ref.aim_marker
        # 信息面/可视化面确实随 B_MAX 变:代表 b 与环半径
        assert sc.beam.b == pytest.approx(0.7 * float(b_max))
        assert sc.b_max == pytest.approx(float(b_max))


def test_fixed_b_offsets_projectile_under_bmax():
    """固定 b(N_B=0)按真值偏移:弹丸横向位置=b, 随 B_MAX 变;靶恒在原点。"""
    for b_max in ("1.0", "2.0", "3.5"):
        sc = _pair_scene_bmax(0, b_max)
        r_sep = sc.beam.p_site[2]
        assert sc.frag_coms[1] == pytest.approx([0.0, 0.0, 0.0])  # 靶不动
        # χ 未写 → 0(法向入射): 偏移沿 +x, p_site = (b, 0, R_SEP)
        assert sc.frag_coms[0] == pytest.approx([float(b_max), 0.0, r_sep])
        assert sc.aim_marker == pytest.approx([float(b_max), 0.0, 0.0])
        assert sc.beam.b == pytest.approx(float(b_max))  # 信息面与显示一致


def test_headon_notes_semantics():
    """注记:采样=恒对心语义;固定 b=偏移语义;旧代表值措辞不在场。"""
    sc = _pair_scene_bmax(1, "2.0")
    assert any("恒对心" in n and "范围见环" in n for n in sc.notes)
    assert not any("0.7·B_MAX" in n for n in sc.notes)
    sc0 = _pair_scene_bmax(0, "2.0")
    assert any("N_B=0" in n and "偏移" in n for n in sc0.notes)
    # 表面(B_MAX 法定不写 → b=0):除注记措辞外零行为变化,b_max 不设环
    sc_s = preview.assemble_preview(
        sysdir.scan_system_dir(ROOT / "cases" / "case_surface_window" / "system"),
        dict(parse_input_file(ROOT / "cases" / "case_surface_window" / "input_qct.txt")),
        i_proj=0,
    )
    assert sc_s.beam is not None
    assert sc_s.b_max is None
    assert sc_s.aim_marker[0] == pytest.approx(sc_s.beam.drop[0])
    assert sc_s.aim_marker[1] == pytest.approx(sc_s.beam.drop[1])
    assert sc_s.beam.b == pytest.approx(0.0)
