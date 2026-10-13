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
    s2 = math.sqrt(0.5)
    assert b.u_inc == pytest.approx([0.0, 0.0, -1.0])
    assert b.p_site == pytest.approx([1.4 * s2, 1.4 * s2, 8.0])
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
