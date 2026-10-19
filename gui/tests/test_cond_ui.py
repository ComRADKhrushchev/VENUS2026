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
    page._w["SCATTERING"].setCurrentText("SURFACE")


def test_paradigm_hides_paradigm_keys():
    """范式决定: 表面藏 B_MAX/靶片段; 气相藏晶格常数/落点模式/定点坐标; 单片段再藏双角色。"""
    page = _page()
    _gas(page)
    for k in ("CELL", "N_AIM", "AIM_X", "AIM_Y"):
        assert page.row_hidden(k), f"气相应隐藏 {k}"
    for k in ("PROJECTILE", "TARGET", "B_MAX", "N_B"):
        assert not page.row_hidden(k), f"气相双片段应显示 {k}"
    _surf(page)
    assert page.row_hidden("B_MAX") and page.row_hidden("TARGET"), "表面应隐藏 B_MAX 与靶片段"
    assert page.row_hidden("N_B"), "表面应隐藏 N_B(b 合法值仅 0, 瞄准模式无意义)"
    for k in ("PROJECTILE", "CELL", "N_AIM"):
        assert not page.row_hidden(k), f"表面应显示 {k}"


def test_hidden_keys_cleared_and_never_collected():
    """隐藏即清空; 对隐藏行程序化塞值也不进 collect(隐藏的键绝不写出)。"""
    page = _page()
    _gas(page)
    page._w["AIM_X"].setText("3.0")  # 气相下 AIM_X 整行隐藏
    assert "AIM_X" not in page.collect()[0], "隐藏键不得进入 collect"
    _surf(page)
    page.set_choice("N_AIM", "2")
    assert not page.row_hidden("AIM_X"), "N_AIM=2 应显示定点坐标"
    page._w["AIM_X"].setText("3.0")
    page.set_choice("N_AIM", "0")
    assert page.row_hidden("AIM_X") and page._w["AIM_X"].text() == "", "隐藏时刻清空"
    assert "AIM_X" not in page.collect()[0]


def test_n_aim_governs_aim_rows():
    page = _page()
    _surf(page)
    assert page.row_hidden("AIM_X") and page.row_hidden("AIM_Y"), "落点模式默认(0)不显示定点"
    page.set_choice("N_AIM", "2")
    assert not page.row_hidden("AIM_X") and not page.row_hidden("AIM_Y")
    page.set_choice("N_AIM", "1")
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
    page.set_choice("N_E_REL", "0")
    assert not page.row_hidden("T_TRANS") and page.row_hidden("V_WIDTH")
    page.set_choice("N_E_REL", "2")
    assert page.row_hidden("T_TRANS") and not page.row_hidden("V_WIDTH")
    assert not page.row_hidden("CHI"), "方位角默认固定: 显示固定方位角"
    page.set_choice("N_CHI", "1")
    assert page.row_hidden("CHI"), "方位角均匀采样: 不需要固定值"
    assert not page.row_hidden("THTA_MAX"), "极角两种模式共用, 恒显"


def test_member_params_inline_per_fragment():
    """逐片段内联参数(2026-10-13 括号裁定): 选词后该词参数字段内联出现,
    换词即重建;无参数字段不属于 _w(参数不是顶层键)。"""
    page = _page()
    _gas(page)
    c0 = page._scheme_combos[0]
    c0.setCurrentText("BARRIER_EXCITATION")
    fields0 = page._scheme_params[0]
    assert set(fields0) == {"E_STAB", "N_E_BAR", "E_BAR", "T_BAR"}, \
        "鞍点词四参数内联"
    assert "E_BAR" not in page._w, "参数不是顶层键"
    c0.setCurrentText("WIGNER")
    assert set(page._scheme_params[0]) == {"E_SCALE"}, "换词即重建"
    # 片段 2 独立选词,互不影响
    page._scheme_combos[1].setCurrentText("EBK")
    assert set(page._scheme_params[1]) == {"N_VIB", "N_ROT"}
    assert set(page._scheme_params[0]) == {"E_SCALE"}


def test_member_labels_physical():
    """内联参数标签 = 参数名, 物理量说明入 tooltip;collect 构造括号值。"""
    page = _page()
    _gas(page)
    page._scheme_combos[0].setCurrentText("WIGNER")
    page._scheme_combos[1].setCurrentText("ROTATION")   # 每片段都须选词才写出
    page._scheme_params[0]["E_SCALE"].setText("1.0")
    vals, _c = page.collect()
    assert "WIGNER[E_SCALE=1.0]" in vals["DIST_SCHEME"]


def test_mark_stars_physical_names():
    """星号标注加粗物理量名(单位后缀不进粗体)。"""
    page = _page()
    assert page._labels["E_REL"].text() == "<b>碰撞能 *</b> [eV]"
    _surf(page)
    assert page._labels["CELL"].text() == "<b>晶胞矢量 *</b> [Å]"
    assert page._labels["N_AIM"].text() == "落点模式"


def test_retired_keys_absent():
    """2026-10-13 用户裁决: GUI 表单移除主势能面(PES_MAIN)/势能面数(N_SURF)/
    动能·势能·总能输出键(E_KIN/E_POT/E_TOT)——键仍在 input.f90 键表与 keyspec
    镜像中(Fortran 侧不变), GUI 不再提供入口, 运行走键表默认。"""
    page = _page()
    for k in ("N_SURF", "PES_MAIN", "E_KIN", "E_POT", "E_TOT"):
        assert k not in page._w, f"{k} 不应再有输入控件"
        assert k not in page._labels, f"{k} 不应再有行标签"
        assert k not in page._cell_of_key, f"{k} 不应再占表单行"
        assert k not in page._checks, f"{k} 不应再有复选框"
    vals, _ = page.collect()
    for k in ("N_SURF", "PES_MAIN", "E_KIN", "E_POT", "E_TOT"):
        assert k not in vals


def test_observables_default_prefill():
    """记录列默认预填 动能/势能/总能(2026-10-14 裁决: t_kin/v_pot/e_tot,
    recorder.f90 注册表实名); 清空即回到不写出(程序默认)。"""
    page = _page()
    assert page._w["OBSERVABLES"].text() == "t_kin,v_pot,e_tot"
    vals, _ = page.collect()
    assert vals.get("OBSERVABLES") == "t_kin,v_pot,e_tot"
    page._w["OBSERVABLES"].clear()
    vals2, _ = page.collect()
    assert "OBSERVABLES" not in vals2, "清空后不应写出该键"


def test_member_section_removed():
    """方案词参数说明组已删(2026-10-14: 参数内联后仅剩一行说明的残留组);
    说明并入采样方案组悬停。"""
    page = _page()
    assert "member" not in page._sections
    tip = page._sections["scheme"]._btn.toolTip()
    assert "内联" in tip and "裸词" in tip


def test_cell_autofill_from_poscar(tmp_path):
    """表面散射 CELL 预填(2026-10-14): 扫描到带胞 POSCAR 且 CELL 为空→自动预填
    (键仍权威); 手动清空不被刷新抢填; 范式切换隐藏清空后切回可复填。"""
    d = tmp_path / "system"
    d.mkdir()
    (d / "slab.poscar").write_text('slab\n1.0\n9.0 0.0 0.0\n0.0 9.0 0.0\n0.0 0.0 12.0\nH\n2\nCartesian\n0.5 0.5 0.1\n5.0 5.0 0.2', encoding="utf-8")
    page = _page()
    _gas(page)
    for k in ("CELL", "N_AIM", "AIM_X", "AIM_Y"):
        assert page.row_hidden(k), f"气相应隐藏 {k}"
    for k in ("PROJECTILE", "TARGET", "B_MAX", "N_B"):
        assert not page.row_hidden(k), f"气相双片段应显示 {k}"
    _surf(page)
    assert page.row_hidden("B_MAX") and page.row_hidden("TARGET"), "表面应隐藏 B_MAX 与靶片段"
    assert page.row_hidden("N_B"), "表面应隐藏 N_B(b 合法值仅 0, 瞄准模式无意义)"
    for k in ("PROJECTILE", "CELL", "N_AIM"):
        assert not page.row_hidden(k), f"表面应显示 {k}"


def test_hidden_keys_cleared_and_never_collected():
    """隐藏即清空; 对隐藏行程序化塞值也不进 collect(隐藏的键绝不写出)。"""
    page = _page()
    _gas(page)
    page._w["AIM_X"].setText("3.0")  # 气相下 AIM_X 整行隐藏
    assert "AIM_X" not in page.collect()[0], "隐藏键不得进入 collect"
    _surf(page)
    page.set_choice("N_AIM", "2")
    assert not page.row_hidden("AIM_X"), "N_AIM=2 应显示定点坐标"
    page._w["AIM_X"].setText("3.0")
    page.set_choice("N_AIM", "0")
    assert page.row_hidden("AIM_X") and page._w["AIM_X"].text() == "", "隐藏时刻清空"
    assert "AIM_X" not in page.collect()[0]


def test_n_aim_governs_aim_rows():
    page = _page()
    _surf(page)
    assert page.row_hidden("AIM_X") and page.row_hidden("AIM_Y"), "落点模式默认(0)不显示定点"
    page.set_choice("N_AIM", "2")
    assert not page.row_hidden("AIM_X") and not page.row_hidden("AIM_Y")
    page.set_choice("N_AIM", "1")
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
    page.set_choice("N_E_REL", "0")
    assert not page.row_hidden("T_TRANS") and page.row_hidden("V_WIDTH")
    page.set_choice("N_E_REL", "2")
    assert page.row_hidden("T_TRANS") and not page.row_hidden("V_WIDTH")
    assert not page.row_hidden("CHI"), "方位角默认固定: 显示固定方位角"
    page.set_choice("N_CHI", "1")
    assert page.row_hidden("CHI"), "方位角均匀采样: 不需要固定值"
    assert not page.row_hidden("THTA_MAX"), "极角两种模式共用, 恒显"


def test_member_params_inline_per_fragment():
    """逐片段内联参数(2026-10-13 括号裁定): 选词后该词参数字段内联出现,
    换词即重建;无参数字段不属于 _w(参数不是顶层键)。"""
    page = _page()
    _gas(page)
    c0 = page._scheme_combos[0]
    c0.setCurrentText("BARRIER_EXCITATION")
    fields0 = page._scheme_params[0]
    assert set(fields0) == {"E_STAB", "N_E_BAR", "E_BAR", "T_BAR"}, \
        "鞍点词四参数内联"
    assert "E_BAR" not in page._w, "参数不是顶层键"
    c0.setCurrentText("WIGNER")
    assert set(page._scheme_params[0]) == {"E_SCALE"}, "换词即重建"
    # 片段 2 独立选词,互不影响
    page._scheme_combos[1].setCurrentText("EBK")
    assert set(page._scheme_params[1]) == {"N_VIB", "N_ROT"}
    assert set(page._scheme_params[0]) == {"E_SCALE"}


def test_member_labels_physical():
    """内联参数标签 = 参数名, 物理量说明入 tooltip;collect 构造括号值。"""
    page = _page()
    _gas(page)
    page._scheme_combos[0].setCurrentText("WIGNER")
    page._scheme_combos[1].setCurrentText("ROTATION")   # 每片段都须选词才写出
    page._scheme_params[0]["E_SCALE"].setText("1.0")
    vals, _c = page.collect()
    assert "WIGNER[E_SCALE=1.0]" in vals["DIST_SCHEME"]


def test_mark_stars_physical_names():
    """星号标注加粗物理量名(单位后缀不进粗体)。"""
    page = _page()
    assert page._labels["E_REL"].text() == "<b>碰撞能 *</b> [eV]"
    _surf(page)
    assert page._labels["CELL"].text() == "<b>晶胞矢量 *</b> [Å]"
    assert page._labels["N_AIM"].text() == "落点模式"


def test_retired_keys_absent():
    """2026-10-13 用户裁决: GUI 表单移除主势能面(PES_MAIN)/势能面数(N_SURF)/
    动能·势能·总能输出键(E_KIN/E_POT/E_TOT)——键仍在 input.f90 键表与 keyspec
    镜像中(Fortran 侧不变), GUI 不再提供入口, 运行走键表默认。"""
    page = _page()
    for k in ("N_SURF", "PES_MAIN", "E_KIN", "E_POT", "E_TOT"):
        assert k not in page._w, f"{k} 不应再有输入控件"
        assert k not in page._labels, f"{k} 不应再有行标签"
        assert k not in page._cell_of_key, f"{k} 不应再占表单行"
        assert k not in page._checks, f"{k} 不应再有复选框"
    vals, _ = page.collect()
    for k in ("N_SURF", "PES_MAIN", "E_KIN", "E_POT", "E_TOT"):
        assert k not in vals


def test_observables_default_prefill():
    """记录列默认预填 动能/势能/总能(2026-10-14 裁决: t_kin/v_pot/e_tot,
    recorder.f90 注册表实名); 清空即回到不写出(程序默认)。"""
    page = _page()
    assert page._w["OBSERVABLES"].text() == "t_kin,v_pot,e_tot"
    vals, _ = page.collect()
    assert vals.get("OBSERVABLES") == "t_kin,v_pot,e_tot"
    page._w["OBSERVABLES"].clear()
    vals2, _ = page.collect()
    assert "OBSERVABLES" not in vals2, "清空后不应写出该键"


def test_member_section_removed():
    """方案词参数说明组已删(2026-10-14: 参数内联后仅剩一行说明的残留组);
    说明并入采样方案组悬停。"""
    page = _page()
    assert "member" not in page._sections
    tip = page._sections["scheme"]._btn.toolTip()
    assert "内联" in tip and "裸词" in tip
