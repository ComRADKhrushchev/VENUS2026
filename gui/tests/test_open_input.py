# -*- coding: utf-8 -*-
"""打开输入文件批(2026-10-14 用户诉求): 现有 input_qct.txt 一键回填全表单;
"体系文件夹"更名"反应物构型"(标签/占位符/悬停/判词; 写出键 SYSTEM_DIR 锚不动)。"""
import os
import shutil

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


def test_reactant_label_renamed():
    """标签更名且写出键锚保留; 占位符与信息行初始措辞随改。"""
    page = _page()
    assert page._labels["SYSTEM_DIR"].text() == "<b>反应物构型 *</b>"
    tip = page._labels["SYSTEM_DIR"].toolTip()
    assert "写出键: SYSTEM_DIR" in tip and "反应物构型" in tip
    assert page._sd_line.placeholderText() == "反应物构型路径"
    assert page._info.text() == "未扫描反应物构型"


def test_open_input_file_populates(tmp_path, monkeypatch):
    """端到端回填: 金标案例复制为相对 SYSTEM_DIR 形("system")——构型行回位/
    方案词按序/全部写出键落于对应控件/容器键入表/装载折叠组自动展开。"""
    from gui.core import inputwrite

    case = ROOT / "cases" / "case_harmonic_scatter"
    dst = tmp_path / "input_qct.txt"
    lines = []
    for ln in (case / "input_qct.txt").read_text(encoding="utf-8").splitlines():
        if ln.strip().upper().startswith("SYSTEM_DIR="):
            ln = "SYSTEM_DIR=system"
        lines.append(ln)
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")
    shutil.copytree(case / "system", tmp_path / "system")
    parsed = inputwrite.parse_input_file(dst)

    monkeypatch.setattr(
        "gui.pages.input_page.QFileDialog.getOpenFileName",
        staticmethod(lambda *a, **k: (str(dst), "")))

    page = _page()
    page._open_input_file()

    # 相对形按输入文件位置回位; 扫描成功
    assert page._sd_line.text() == str(tmp_path / "system")
    assert page._info.text().startswith("✓"), page._info.text()
    # 方案词按片段序回填(括号参数形: 词入下拉, 参数入内联字段)
    from gui.core.keyspec import parse_dist_scheme

    items = parse_dist_scheme(parsed.get("DIST_SCHEME", "")) \
        if parsed.get("DIST_SCHEME") else []
    assert [c.currentText() for c in page._scheme_combos] == [w for w, _ in items]
    for i, (w, params) in enumerate(items):
        if i < len(page._scheme_params):
            for name, e in page._scheme_params[i].items():
                if name in params:
                    assert e.text() == params[name], f"片段{i} 参数 {name}"
    # 全键巡回: 每个非空写出键落于对应控件(复选/角色键文本/行编辑文本/下拉数据值/容器表)
    for k, v in parsed.items():
        if k in ("SYSTEM_DIR", "DIST_SCHEME") or not v:
            continue
        w = page._w.get(k)
        if k in page._checks:
            assert page._checks[k].isChecked() == v.strip().upper().startswith("T"), k
        elif w is None:  # 容器键 → 表内行值
            t = page._container_table
            hit = [t.item(r, 1).text() if t.item(r, 1) else ""
                   for r in range(t.rowCount())
                   if t.item(r, 0) and t.item(r, 0).text().strip().upper() == k]
            assert hit == [v], f"容器键 {k}={v} 未入表: {hit}"
        elif k in ("PROJECTILE", "TARGET"):
            assert w.currentText() == v, k
        elif hasattr(w, "text"):
            assert w.text() == v, k
        else:
            assert w.currentData() == v, f"{k}: {w.currentData()!r} != {v!r}"
    # 装载到折叠组的键自动展开(该组有键入档时)
    for sid, sec in page._sections.items():
        keys = {k for k, s in page._section_of_key.items() if s is sec}
        if sid not in ("sysdef", "scheme", "output") and keys & set(parsed):
            assert sec.expanded, f"{sid} 组装载后应自动展开"
