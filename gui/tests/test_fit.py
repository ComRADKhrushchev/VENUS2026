# -*- coding: utf-8 -*-
"""适配批(2026-10-13 用户报障): 窗口下沿不被任务栏挡住 + 左栏免横向滚动。

  - 主窗口最小高度 ≤ 700px(页2 包滚动区塌缩最小尺寸), 初始尺寸钳制在
    屏幕可用几何内(availableGeometry, 任务栏除外), 打开即完整可见;
  - 左栏内容最小宽 ≤ 560px(默认分栏份额), 平时无需横向滚动;
  - 状态栏单位文案与 eV/Å/amu/fs 新制一致(2026-10-13 改基)。
"""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")
pytest.importorskip("pyqtgraph")

_APP = None


def _app():
    global _APP
    from PyQt6.QtWidgets import QApplication

    if QApplication.instance() is None:
        _APP = QApplication([])
    return QApplication.instance()


def test_window_min_height_small():
    """窗口最小高度须 ≤700px(1080p 去任务栏 ~1030px 可完整容纳)。"""
    _app()
    from gui.main import MainWindow

    win = MainWindow()
    assert win.minimumSizeHint().height() <= 700, (
        f"窗口最小高 {win.minimumSizeHint().height()}px 仍会顶出屏外")
    win.close()


def test_initial_size_within_available_geometry():
    """初始尺寸钳制在屏幕可用几何内(下沿不压任务栏)。"""
    app = _app()
    from gui.main import MainWindow

    win = MainWindow()
    avail = app.primaryScreen().availableGeometry()
    assert win.width() <= avail.width(), f"宽 {win.width()} > 可用 {avail.width()}"
    assert win.height() <= avail.height(), f"高 {win.height()} > 可用 {avail.height()}"
    win.close()


def test_left_content_min_width_fits_default_share():
    """左栏内容最小宽 ≤560px(默认分栏份额), 平时不出横向滚动条。"""
    _app()
    from PyQt6.QtWidgets import QScrollArea
    from gui.pages.input_page import InputPage

    page = InputPage()
    scr = page.findChildren(QScrollArea)[0]
    w = scr.widget().minimumSizeHint().width()
    assert w <= 560, f"左栏内容最小宽 {w}px > 默认份额 560px, 必出横向滚动条"


def test_pes_page_min_height_small():
    """页2 包滚动区后最小高度塌缩(不再是窗口最小高的元凶)。"""
    _app()
    from gui.pages.pes_page import PesPage

    p = PesPage()
    assert p.minimumSizeHint().height() <= 400, (
        f"页2 最小高 {p.minimumSizeHint().height()}px 应由滚动区塌缩")


def test_statusbar_units_text_current():
    """状态栏单位文案随 2026-10-13 eV/Å/amu/fs 改基(不得残留 10 fs 旧制)。"""
    _app()
    from gui.main import MainWindow

    win = MainWindow()
    msg = win.statusBar().currentMessage()
    assert "eV" in msg and "10 fs" not in msg, f"状态栏单位文案陈旧: {msg!r}"
    win.close()
