# -*- coding: utf-8 -*-
"""条件显隐(隐藏即清空/不写出)、模式下拉化、物理量标签——三诉求的行为契约。

规约(2026-10-13 用户三诉求 + 实现裁决):
  - 被其余设置判"不需要"的行整行隐藏;隐藏即清空该键值(与范式切换清空同规约),
    且 collect 跳过隐藏键——隐藏的键绝不写出(折叠组仍收集, 两机制互不干扰);
  - 成员参数的词基显隐(未选词)保持旧语义(不清空、可收集), 与条件显隐分开记账;
  - ORDER 为随积分器变的下拉(VERLET 1/2·SYMPLECTIC 4/6/8·RADAU 1-15,
    权威 src/dynamics/propagator.f90:126-165/186/223-224);
  - 行标签显示物理量, KEY 名入 tooltip("写出键: X")。
"""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")
pytest.importorskip("pyqtgraph")

from conftest import ROOT

_APP = None


def _page():
    global _APP
    from PyQt6.QtWidgets import QApplication
    from gui.pages.input_page import InputPage

    if QApplication.instance() is None:
        _APP = QApplication([])
    return InputPage()


def _gas(page):
    page._sd_line.setText(str(ROOT / "cases" / "case_harmonic_scatter" / "system"))
    page._rescan()


def _surf(page):
    page._sd_line.setText(str(ROOT / "cases" / "case_spectrum_source" / "h2ag_ebk" / "system"))
    page._rescan()


def test_paradigm_hides_paradigm_keys():
    """范式决定: 表面藏 B_MAX/靶片段; 气相藏晶格常数/落点模式/定点坐标; 单片段再藏双角色。"""
    page = _page()
    _gas(page)
    for k in ("A_LAT", "N_AIM", "AIM_X", "AIM_Y"):
        assert page.row_hidden(k), f"气相应隐藏 {k}"
    for k in ("PROJECTILE", "TARGET", "B_MAX"):
        assert not page.row_hidden(k), f"气相双片段应显示 {k}"
    _surf(page)
    assert page.row_hidden("B_MAX") and page.row_hidden("TARGET"), "表面应隐藏 B_MAX 与靶片段"
    for k in ("PROJECTILE", "A_LAT", "N_AIM"):
        assert not page.row_hidden(k), f"表面应显示 {k}"


def test_hidden_keys_cleared_and_never_collected():
    """隐藏即清空; 对隐藏行程序化塞值也不进 collect(隐藏的键绝不写出)。"""
    page = _page()
    _gas(page)
    page._w["AIM_X"].setText("3.0")  # 气相下 AIM_X 整行隐藏
    assert "AIM_X" not in page.collect()[0], "隐藏键不得进入 collect"
    _surf(page)
    page._w["N_AIM"].setCurrentText("2")
    assert not page.row_hidden("AIM_X"), "N_AIM=2 应显示定点坐标"
    page._w["AIM_X"].setText("3.0")
    page._w["N_AIM"].setCurrentText("0")
    assert page.row_hidden("AIM_X") and page._w["AIM_X"].text() == "", "隐藏时刻清空"
    assert "AIM_X" not in page.collect()[0]


def test_n_aim_governs_aim_rows():
    page = _page()
    _surf(page)
    assert page.row_hidden("AIM_X") and page.row_hidden("AIM_Y"), "落点模式默认(0)不显示定点"
    page._w["N_AIM"].setCurrentText("2")
    assert not page.row_hidden("AIM_X") and not page.row_hidden("AIM_Y")
    page._w["N_AIM"].setCurrentText("1")
    assert page.row_hidden("AIM_X") and page.row_hidden("AIM_Y")


def test_integrator_governs_order_and_adaptive():
    """积分阶恒显但选项随积分器变; 自适应步长界仅 RADAU(dt_min 只被 radau_step 消费)。"""
    from PyQt6.QtWidgets import QComboBox

    page = _page()
    assert isinstance(page._w["ORDER"], QComboBox), "ORDER 应为下拉"
    assert not page.row_hidden("ORDER")
    assert page.row_hidden("DT_MIN") and page.row_hidden("DT_MAX"), "VERLET 默认无自适应界"

    def items(key):
        c = page._w[key]
        return [c.itemText(i) for i in range(c.count())]

    page._w["INTEGRATOR"].setCurrentText("SYMPLECTIC")
    assert items("ORDER") == ["", "4", "6", "8"]
    page._w["INTEGRATOR"].setCurrentText("RADAU")
    assert items("ORDER") == [""] + [str(i) for i in range(1, 16)]
    assert not page.row_hidden("DT_MIN") and not page.row_hidden("DT_MAX")
    page._w["INTEGRATOR"].setCurrentText("VERLET")
    assert items("ORDER") == ["", "1", "2"]
    assert page.row_hidden("DT_MIN") and page.row_hidden("DT_MAX")


def test_bath_and_spectrum_governs_rows():
    page = _page()
    assert all(page.row_hidden(k) for k in ("N_BATH", "DT_BATH", "T_BATH", "NU_COLL")), \
        "热浴 NONE 时四键全藏"
    page._w["BATH"].setCurrentText("ANDERSEN")
    assert not any(page.row_hidden(k) for k in ("N_BATH", "DT_BATH", "T_BATH", "NU_COLL"))
    assert page.row_hidden("SPECTRUM_FILE"), "谱来源默认 COMPUTE 不显示谱表文件"
    page._w["SPECTRUM_SOURCE"].setCurrentText("MANUAL")
    assert not page.row_hidden("SPECTRUM_FILE")


def test_beam_modes_governs_rows():
    page = _page()
    assert page.row_hidden("T_TRANS") and page.row_hidden("V_WIDTH"), \
        "碰撞能默认定值: 平动温度与速度展宽都不需要"
    page._w["N_E_REL"].setCurrentText("0")
    assert not page.row_hidden("T_TRANS") and page.row_hidden("V_WIDTH")
    page._w["N_E_REL"].setCurrentText("2")
    assert page.row_hidden("T_TRANS") and not page.row_hidden("V_WIDTH")
    assert not page.row_hidden("CHI"), "方位角默认固定: 显示固定方位角"
    page._w["N_CHI"].setCurrentText("1")
    assert page.row_hidden("CHI"), "方位角均匀采样: 不需要固定值"
    assert not page.row_hidden("THTA_MAX"), "极角两种模式共用, 恒显"


def test_member_modes_governs_rows():
    """成员内模式: 振子默认热→温度; 1 定态→能级; 反应坐标 1 热→温度; 未选词整行藏(旧语义)。"""
    page = _page()
    _gas(page)
    for c in page._scheme_combos:
        c.setCurrentText("SURFACE_OSCILLATOR")
    assert not page.row_hidden("T_OSC") and page.row_hidden("N_LEVEL"), \
        "振子模式默认 0(热): 显示温度, 藏能级"
    assert not page.row_hidden("K_OSC"), "力常数恒显"
    page._w["N_OSC"].setCurrentText("1")
    assert page.row_hidden("T_OSC") and not page.row_hidden("N_LEVEL")
    for c in page._scheme_combos:
        c.setCurrentText("BARRIER_EXCITATION")
    assert page.row_hidden("T_OSC") and page.row_hidden("K_OSC"), "换词后旧词成员行隐藏"
    page._w["N_E_BAR"].setCurrentText("1")
    assert not page.row_hidden("T_BAR")
    page._w["N_E_BAR"].setCurrentText("0")
    assert page.row_hidden("T_BAR")
    assert not page.row_hidden("E_STAB") and not page.row_hidden("E_BAR")


def test_physical_labels_and_tooltips():
    """行标签显示物理量不含 KEY 名; KEY 名入 tooltip。"""
    page = _page()
    for key, text in (
        ("A_LAT", "晶格常数 [Å]"),
        ("N_AIM", "落点模式"),
        ("T_TRANS", "平动温度 [K]"),
        ("ORDER", "积分阶"),
        ("CHI", "固定方位角 [rad]"),
        ("BATH", "热浴"),
    ):
        lab = page._labels[key]
        assert lab.text() == text, f"{key} 标签应为物理量 {text!r}, 实为 {lab.text()!r}"
        assert key in lab.toolTip(), f"{key} 的 KEY 名应入 tooltip"


def test_member_labels_physical():
    """成员参数标签从 hint 派生物理量名(剥去尾部括注), KEY 名入 tooltip。"""
    page = _page()
    _gas(page)
    page._scheme_combos[0].setCurrentText("THERMALIZE")
    assert page._labels["T_EQ"].text() == "平衡温度 [K]"
    assert "T_EQ" in page._labels["T_EQ"].toolTip()
    page._scheme_combos[0].setCurrentText("EBK")
    assert page._labels["N_VIB"].text() == "EBK 振动量子数"


def test_mark_stars_physical_names():
    """星号标注加粗物理量名(单位后缀不进粗体)。"""
    page = _page()
    assert page._labels["E_REL"].text() == "<b>碰撞能 *</b> [eV]"
    _surf(page)
    assert page._labels["A_LAT"].text() == "<b>晶格常数 *</b> [Å]"
    assert page._labels["N_AIM"].text() == "落点模式"
