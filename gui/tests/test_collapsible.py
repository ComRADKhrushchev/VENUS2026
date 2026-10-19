# -*- coding: utf-8 -*-
"""可折叠参数分组:默认态 / 标题切换 / 折叠不改变键值语义 / 校验出错自动展开。"""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")
pytest.importorskip("pyqtgraph")

from conftest import ROOT


_APP = None  # 无父窗口的 QApplication 须挂引用, 否则被回收时连带删除全部控件


def _page():
    global _APP
    from PyQt6.QtWidgets import QApplication
    from gui.pages.input_page import InputPage

    if QApplication.instance() is None:
        _APP = QApplication([])
    return InputPage()


def test_section_defaults():
    """体系定义/采样方案/输出恒展开不可折叠;运行控制可折叠默认展开;四组默认折叠。"""
    page = _page()
    secs = getattr(page, "_sections", None)
    assert secs is not None, "缺少 _sections(可折叠分组未实现)"
    for sid in ("sysdef", "scheme", "output"):
        assert secs[sid].expanded, f"{sid} 应恒展开"
        assert not secs[sid].collapsible, f"{sid} 应不可折叠"
    assert secs["run"].collapsible and secs["run"].expanded, "运行控制可折叠且默认展开"
    for sid in ("bath", "beam", "container"):
        assert secs[sid].collapsible, f"{sid} 应可折叠"
        assert not secs[sid].expanded, f"{sid} 应默认折叠"


def test_section_toggle():
    """toggle/点击标题切换内容显隐;不可折叠组忽略切换。"""
    page = _page()
    secs = page._sections
    secs["run"].toggle()
    assert not secs["run"].expanded and secs["run"].content.isHidden()
    secs["run"].toggle()
    assert secs["run"].expanded and not secs["run"].content.isHidden()
    # 不可折叠组:切换无效,内容恒可见
    secs["sysdef"].toggle()
    assert secs["sysdef"].expanded


def test_collapsed_values_still_collected():
    """折叠只藏可见性不禁用:折叠组内已填值仍进入 collect(条件隐藏键则绝不写出)。"""
    page = _page()
    assert not page._sections["beam"].expanded
    page._w["THTA_MAX"].setText("1.2")  # 极角: 恒显键, 折叠组语义的有效探针
    vals, _ = page.collect()
    assert vals.get("THTA_MAX") == "1.2"


def test_scheme_word_shows_inline_params():
    """选中带参数的方案词 → 该片段行内联参数字段在场(2026-10-14: 说明组已删,
    取代旧"成员组自动展开"语义)。"""
    page = _page()
    page._sd_line.setText(str(ROOT / "cases" / "case_harmonic_scatter" / "system"))
    page._rescan()
    pc, _pl = page._scheme_param_rows[0]
    assert pc.isHidden()
    page._scheme_combos[0].setCurrentText("WIGNER")
    assert not pc.isHidden()
    assert page._scheme_params[0], "WIGNER 应带内联参数字段"


def test_validation_error_expands_section():
    """校验出错的键所在折叠组自动展开(N_E_REL=0 缺 T_TRANS → 束流组展开)。"""
    page = _page()
    page._sd_line.setText(str(ROOT / "cases" / "case_harmonic_scatter" / "system"))
    page._rescan()
    page._w["PROJECTILE"].setCurrentText("frag_a")
    page._w["TARGET"].setCurrentText("frag_b")
    page._w["E_REL"].setText("200.0")
    page._w["R_SEP"].setText("8.0")
    page._w["DT"].setText("0.05")
    page._w["B_MAX"].setText("2.0")
    page._w["N_B"].setCurrentText("1")
    page.set_choice("N_E_REL", "0")  # 热分布须 T_TRANS>0,此处留空 → 报错
    for c in page._scheme_combos:
        c.setCurrentText("WIGNER")
    assert not page._sections["beam"].expanded
    page.run_validation()
    assert "T_TRANS" in page._verdict.toPlainText()
    assert page._sections["beam"].expanded
