# -*- coding: utf-8 -*-
"""整程序 offscreen 冒烟:两页构建、扫描/预览联动、截图落盘人工复核。"""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")
pytest.importorskip("pyqtgraph")

from conftest import ROOT

from gui.main import MainWindow
from gui.view3d.preview_view import PreviewView


def _shot(win, path):
    win.repaint()
    pix = win.grab()
    pix.save(str(path))
    assert path.is_file() and path.stat().st_size > 10_000, f"截图过小: {path}"


def test_main_window_smoke(tmp_path):
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    win = MainWindow()
    win.show()
    _shot(win, tmp_path / "smoke_main.png")

    # 页1:扫描 harmonic_scatter → 填必填 → 校验/预览(2D 投影后备可离屏渲染)
    page1 = win.centralWidget().widget(0)
    page1._view.set_mode(PreviewView.MODE_2D)
    # 星号标注约定:无条件必填静态加粗带 *;条件必填随范式增减
    assert page1._labels["E_REL"].text() == "<b>E_REL *</b> [kcal/mol]"
    assert page1._labels["A_LAT"].text() == "A_LAT [Å]"      # 表面专用,未扫描不带星
    page1._sd_line.setText(str(ROOT / "cases" / "case_harmonic_scatter" / "system"))
    page1._rescan()
    assert "气相散射" in page1._info.text()
    assert page1._labels["PROJECTILE"].text() == "<b>PROJECTILE *</b>"
    assert page1._labels["TARGET"].text() == "<b>TARGET *</b>"
    page1._w["PROJECTILE"].setCurrentText("frag_a")
    page1._w["TARGET"].setCurrentText("frag_b")
    page1._w["E_REL"].setText("200.0")
    page1._w["R_SEP"].setText("8.0")
    page1._w["DT"].setText("0.05")
    page1._w["B_MAX"].setText("2.0")
    page1._w["N_B"].setCurrentText("1")
    for c in page1._scheme_combos:
        c.setCurrentText("THERMALIZE" if not c.currentIndex() else "ROTATION")
    page1._scheme_combos[1].setCurrentText("ROTATION")
    page1.refresh_preview()
    assert "示意" in page1._notes.text()
    page1.run_validation()
    txt = page1._verdict.toPlainText()
    assert txt.startswith("✓"), txt
    _shot(win, tmp_path / "smoke_input_pair.png")

    # 页1:表面散射切换(h2ag_ebk)
    page1._sd_line.setText(str(ROOT / "cases" / "case_spectrum_source" / "h2ag_ebk" / "system"))
    page1._rescan()
    assert "表面散射" in page1._info.text()
    assert not page1._w["B_MAX"].isEnabled()
    assert page1._w["B_MAX"].placeholderText() == "表面散射不支持"
    assert page1._labels["A_LAT"].text() == "<b>A_LAT *</b> [Å]"
    assert page1._labels["TARGET"].text() == "TARGET"        # 表面散射:无 TARGET 星
    page1._w["PROJECTILE"].setCurrentText("H2")
    page1._w["A_LAT"].setText("2.9437")
    page1._w["N_AIM"].setCurrentText("1")
    for c in page1._scheme_combos:
        c.setCurrentText("EBK")
    page1.refresh_preview()
    page1.run_validation()
    txt = page1._verdict.toPlainText()
    assert txt.startswith("✓"), txt
    _shot(win, tmp_path / "smoke_input_surface.png")

    # 页2:扫描 auto_sample → 十名全在场
    page2 = win.centralWidget().widget(1)
    page2._dir.setText(str(ROOT / "systems" / "inactive" / "auto_sample"))
    page2.scan()
    assert page2._scan is not None and page2._scan.ok, page2._scan.errors if page2._scan else None
    assert "契约满足" in page2._scan_info.text()
    _shot(win, tmp_path / "smoke_pes.png")

    # 页2:负面对照(pot-only 缺 pes 的空壳)
    empty = tmp_path / "empty_shell"
    empty.mkdir()
    (empty / "x.f90").write_text("module x\nend module\n", encoding="utf-8")
    page2._dir.setText(str(empty))
    page2.scan()
    assert not page2._scan.ok
    win.close()


def test_projection2d_pixel_assertions(tmp_path):
    """2D 投影后备的像素级断言:瞄准点橙/箭头蓝/键灰/原子白 + 双面板各占半幅。"""
    from PyQt6.QtCore import QSize
    from PyQt6.QtWidgets import QApplication
    from gui.core import preview, sysdir
    from gui.core.inputwrite import parse_input_file
    from gui.view3d.preview_view import _Projection2D

    app = QApplication.instance() or QApplication([])
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")
    vals = dict(parse_input_file(ROOT / "cases" / "case_harmonic_scatter" / "input_qct.txt"))
    scene = preview.assemble_preview(sd, vals, i_proj=0, i_targ=1)
    w = _Projection2D()
    w.resize(QSize(900, 500))
    w.set_scene(scene)
    img = w.grab().toImage()

    def near(c, t, tol=40):
        return all(abs(c[i] - t[i]) <= tol for i in range(3))

    counts = {"orange_aim": 0, "blue_arrow": 0, "gray_bond": 0}
    left_ink = right_ink = 0
    for y in range(img.height()):
        for x in range(img.width()):
            px = img.pixel(x, y)
            c = ((px >> 16) & 255, (px >> 8) & 255, px & 255)
            if near(c, (255, 140, 26)):
                counts["orange_aim"] += 1
            elif near(c, (26, 128, 255)):
                counts["blue_arrow"] += 1
            elif near(c, (90, 90, 90)):
                counts["gray_bond"] += 1
            if sum(c) < 640:
                if x < img.width() // 2:
                    left_ink += 1
                else:
                    right_ink += 1
    assert counts["orange_aim"] > 30, "瞄准点/盘未渲染"
    assert counts["blue_arrow"] > 30, "速度箭头未渲染"
    assert counts["gray_bond"] > 10, "键线未渲染"
    assert left_ink > 100 and right_ink > 100, "双面板未各占半幅"
