# -*- coding: utf-8 -*-
"""行布局批(2026-10-13 用户裁决: 同类设置一行多列, 自然成组打包):

  - 行 = QHBox 容器, 内含 1-3 个"单元对"(标签+控件); 隐藏某键的单元对时同行
    其余键左移补位(QBoxLayout 不为隐藏件留位); 整行键全隐藏时行容器塌缩;
  - 模式类下拉选项显示"值 · 含义"(如 0 · 热分布), 写出取原值(itemData);
    set_choice(key, value) 为程序化选择入口;
  - 行内提示文字并入悬停提示(标签与控件同载"写出键: X"+提示), 释放横向空间。
语义不变: 条件显隐/折叠/校验联动与 row_hidden API 照旧。
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
    page._w["SCATTERING"].setCurrentText("SURFACE")


def test_natural_packing():
    """自然成组键同行: 坐标对/角色对/三计数键/积分三连/模式+值; 异组不同行。"""
    page = _page()
    same = page._same_row
    assert same("E_REL", "R_SEP"), "碰撞能与初始间距应同行"
    assert same("PROJECTILE", "TARGET"), "弹丸与靶应同行"
    assert same("AIM_X", "AIM_Y"), "定点落点 x/y 应同行"
    assert same("N_TRAJ", "ISEED") and same("ISEED", "MAX_STEPS"), "三计数键应同行"
    assert same("DT", "INTEGRATOR") and same("INTEGRATOR", "ORDER"), "步长/积分器/阶应同行"
    assert same("DT_MIN", "DT_MAX")
    assert same("N_THTA", "THTA_MAX"), "极角模式与极角应同行"
    assert same("N_CHI", "CHI"), "方位角模式与固定方位角应同行"
    assert same("BATH", "SPECTRUM_SOURCE")
    assert same("N_BATH", "DT_BATH")
    assert same("T_BATH", "NU_COLL")
    assert same("N_E_REL", "T_TRANS") and same("T_TRANS", "V_WIDTH")
    assert not same("E_REL", "B_MAX"), "不同性质键不该被拉进同一行"


def test_member_rows_packed_by_word():
    """逐片段参数随方案行内联(2026-10-13 括号裁定): 单词全参数一次铺开,
    无参数词不产生字段;占位符带物理量说明。"""
    page = _page()
    _gas(page)
    page._scheme_combos[0].setCurrentText("BARRIER_EXCITATION")
    f = page._scheme_params[0]
    assert set(f) == {"E_STAB", "N_E_BAR", "E_BAR", "T_BAR"}
    assert f["T_BAR"].placeholderText(), "占位符载有物理量说明"
    page._scheme_combos[0].setCurrentText("ROTATION")
    assert set(page._scheme_params[0]) == {"J_ROT"}


def test_mode_combo_value_display():
    """模式下拉: 选项显示'值 · 含义', 写出取原值; set_choice 为程序化入口。"""
    page = _page()
    c = page._w["N_E_REL"]
    texts = [c.itemText(i) for i in range(c.count())]
    assert any("热" in t and t.startswith("0") for t in texts), f"选项应含 '0 · 热分布', 实为 {texts}"
    page.set_choice("N_E_REL", "0")
    assert page._w["N_E_REL"].currentData() == "0"
    vals, _ = page.collect()
    assert vals.get("N_E_REL") == "0", "写出值须为原值 0 而非显示文本"
    page.set_choice("N_AIM", "1")
    aim = page._w["N_AIM"]
    assert "胞内" in aim.itemText(aim.currentIndex())
    assert page.row_hidden("N_AIM"), "未扫描: N_AIM 条件隐藏"
    _surf(page)  # 表面散射下 N_AIM 在场
    page.set_choice("N_AIM", "1")
    assert page.collect()[0].get("N_AIM") == "1"
    # 数据面不变的普通枚举(INTEGRATOR)仍可 setCurrentText
    page._w["INTEGRATOR"].setCurrentText("RADAU")
    assert page.collect()[0].get("INTEGRATOR") == "RADAU"


def test_hidden_cell_shifts_and_row_collapses():
    """同行键隐藏: 单元对消失(剩者补位); 整行键全隐藏时行容器塌缩。"""
    page = _page()
    page.set_choice("N_E_REL", "0")
    assert not page.row_hidden("T_TRANS") and page.row_hidden("V_WIDTH")
    page.set_choice("N_E_REL", "2")
    assert page.row_hidden("T_TRANS") and not page.row_hidden("V_WIDTH")
    page.set_choice("N_E_REL", "1")
    # 束流行只剩碰撞能模式, 模式键恒显 → 行容器在场
    assert not page.row_hidden("N_E_REL")
    # AIM 行两键全隐藏(表面+落点模式≠2) → 行容器塌缩
    _surf(page)
    assert page.row_hidden("AIM_X") and page.row_hidden("AIM_Y")
    page.set_choice("N_AIM", "2")
    assert not page.row_hidden("AIM_X") and not page.row_hidden("AIM_Y")


def test_hints_folded_into_tooltips():
    """旁注并入悬停: 标签 tooltip 载 '写出键: X' 与提示; 不再有第三列提示标签。"""
    page = _page()
    tip = page._labels["THTA_MAX"].toolTip()
    assert "写出键: THTA_MAX" in tip, f"tooltip 应含键名, 实为 {tip!r}"
    assert "采样界" in tip or "极角" in tip, "tooltip 应并入原旁注约束"
    assert "写出键: DT" in page._labels["DT"].toolTip()
