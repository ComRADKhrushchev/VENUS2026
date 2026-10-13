# -*- coding: utf-8 -*-
"""sysdir 测试:元素质量表对源表、案例 system/ 解析、范式与 box 派生、负面规则。"""
import math
import re
from pathlib import Path

import pytest
from conftest import ROOT

from gui.core import sysdir


def test_elem_mass_mirror():
    """质量表数值逐项等于 specio.f90 DATA 表(程序常数,防转录漂移)。"""
    text = (ROOT / "src" / "utils" / "specio.f90").read_text(encoding="utf-8")
    m = re.search(r"data tbl_symb / (.*?) /", text, re.S)
    syms = re.findall(r"'(\w+)'", m.group(1))
    m2 = re.search(r"data tbl_mass / (.*?) /", text, re.S)
    vals = [float(x.replace("d", "e")) for x in re.findall(r"-?\d+\.\d+(?:d[-+]?\d+)?", m2.group(1))]
    assert syms == sysdir.ELEM_SYMBOLS
    assert vals == sysdir.ELEM_MASSES
    assert sysdir.ELEM_MASS["D"] == 2.014
    assert sysdir.elem_mass("AU") == pytest.approx(196.966569)
    assert sysdir.elem_mass("H") == pytest.approx(1.008)


def test_scan_pair_case():
    """case_harmonic_scatter:双 xyz 片段、名序。"""
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")
    assert sd.ok, sd.errors
    assert sd.frag_names == ["frag_a", "frag_b"]
    assert sd.frags[0].nat == 2 and sd.frags[1].nat == 2
    assert sd.frags[0].coords[1][0] == pytest.approx(0.7)
    assert sd.frags[1].coords[1][0] == pytest.approx(0.8)
    assert sd.frags[0].masses == pytest.approx([1.008, 1.008])


def test_scan_surface_case():
    """case_spectrum_source/h2ag_ebk:H2.poscar 带胞,box 余弦定理(范式自 SCATTERING 键)。"""
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_spectrum_source" / "h2ag_ebk" / "system")
    assert sd.ok, sd.errors
    assert sd.frag_names == ["H2"]
    assert sd.frags[0].nat == 2
    c = sd.cell
    assert c is not None
    # 六方 Ag(111) 原胞 a=2.9437:|a1| = |a2| = 2.9437
    assert sd.box_lx == pytest.approx(2.9437, abs=1e-6)
    assert sd.box_ly == pytest.approx(2.9437, abs=1e-6)
    assert sd.skew == pytest.approx(60.0 * math.pi / 180.0, abs=1e-4)


def test_parse_poscar_cartesian_and_velocity_skip(tmp_path: Path):
    """Cartesian 模式 + 缩放因子 + Selective dynamics + 速度块读跳。"""
    poscar = "\n".join([
        "test slab",
        "0.5",
        "6.0 0.0 0.0",
        "0.0 6.0 0.0",
        "0.0 0.0 9.0",
        "H Au",
        "1 1",
        "Selective dynamics",
        "Cartesian",
        "0.0 0.0 0.0 T T F",
        "1.0 2.0 3.0",
        "0.1 0.0 0.0",
        "0.0 0.1 0.0",
    ])
    frag, cell = sysdir.parse_poscar(poscar, "Slab.poscar")
    assert frag.symbols == ["H", "Au"]
    assert frag.coords[0] == pytest.approx([0.0, 0.0, 0.0])
    assert frag.coords[1] == pytest.approx([0.5, 1.0, 1.5])
    assert cell[0] == pytest.approx([3.0, 0.0, 0.0])


def test_parse_poscar_direct(tmp_path: Path):
    poscar = "\n".join([
        "x", "1.0",
        "10.0 0.0 0.0", "0.0 8.0 0.0", "0.0 0.0 12.0",
        "H", "1",
        "Direct",
        "0.5 0.25 0.1",
    ])
    frag, cell = sysdir.parse_poscar(poscar, "H.poscar")
    assert frag.coords[0] == pytest.approx([5.0, 2.0, 1.2])


def test_negative_rules(tmp_path: Path):
    # 未知文件类型
    d = tmp_path / "sys1"
    d.mkdir()
    (d / "notes.txt").write_text("x", encoding="utf-8")
    sd = sysdir.scan_system_dir(d)
    assert not sd.ok and any("未知体系文件类型" in e for e in sd.errors)
    # 名折叠重复
    d2 = tmp_path / "sys2"
    d2.mkdir()
    (d2 / "a.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")
    (d2 / "A.poscar").write_text("x\n1.0\n2 0 0\n0 2 0\n0 0 2\nH\n1\nCartesian\n0 0 0\n", encoding="utf-8")
    sd2 = sysdir.scan_system_dir(d2)
    assert any("重复" in e for e in sd2.errors)
    # 双胞
    d3 = tmp_path / "sys3"
    d3.mkdir()
    pos = "x\n1.0\n2 0 0\n0 2 0\n0 0 2\nH\n1\nCartesian\n0 0 0\n"
    (d3 / "a.poscar").write_text(pos, encoding="utf-8")
    (d3 / "b.poscar").write_text(pos, encoding="utf-8")
    sd3 = sysdir.scan_system_dir(d3)
    assert any("带胞" in e for e in sd3.errors)
    # 空文件夹
    d4 = tmp_path / "sys4"
    d4.mkdir()
    sd4 = sysdir.scan_system_dir(d4)
    assert not sd4.ok


def test_classify_entry():
    assert sysdir.classify_entry("H2.xyz") == ("xyz", "H2")
    assert sysdir.classify_entry("SLAB.POSCAR") == ("poscar", "SLAB")
    assert sysdir.classify_entry("POSCAR.au111") == ("poscar", "POSCAR")  # 扩展名剥离,基名即反应物名
    with pytest.raises(ValueError):
        sysdir.classify_entry("data.txt")
    with pytest.raises(ValueError):
        sysdir.classify_entry(".xyz")
