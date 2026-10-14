# -*- coding: utf-8 -*-
"""主题批(2026-10-13 用户裁决: 方案 A 令牌化浅色主题 + 图标系统 + 密度重排):

  - 唯一样式入口 = gui/theme.py(令牌 + build_qss + apply_theme), 页面代码
    零 setStyleSheet(局部差异走 objectName/dynamic property);
  - 令牌为纯量常量(无 Qt 可读), 间距标尺 SP_S/M/L/XL = 4/8/12/16;
  - 密度: 页边距/组内边距取自 theme 常量; 左栏横向滚动条关闭且内容适配;
  - 主次按钮: 写出/开始编译 class=primary(实心强调);
  - 校验语义色: _verdict 与摘要行 state=ok/error(文本不变, 仅样式);
  - 分组头 QToolButton objectName=SectionHeader + chevron 图标(箭头 API 退役)。
语义不变: 标签文本/tooltip/显隐/收集/折叠契约由既有测试看守, 本文件只锁样式面。
"""
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")

from conftest import ROOT
from gui import theme

_APP = None


def _app():
    global _APP
    from PyQt6.QtWidgets import QApplication
    if QApplication.instance() is None:
        _APP = QApplication([])
    return QApplication.instance()


def _page():
    from gui.pages.input_page import InputPage
    return InputPage()


# ---------- 令牌与 QSS(无 Qt 依赖面) ----------

def test_tokens_complete():
    need = ("accent", "accent_hover", "accent_pressed", "accent_soft", "surface",
            "card", "band", "border", "text_primary", "text_secondary", "text_muted",
            "success", "error", "tooltip_bg", "tooltip_fg")
    for k in need:
        assert theme.COLORS[k].startswith("#") and len(theme.COLORS[k]) == 7, k


def test_spacing_scale():
    assert (theme.SP_S, theme.SP_M, theme.SP_L, theme.SP_XL) == (4, 8, 12, 16)


def test_qss_key_selectors():
    qss = theme.build_qss()
    for sel in ("QToolButton#SectionHeader", 'QPushButton[class="primary"]',
                'QPlainTextEdit[state="ok"]', 'QPlainTextEdit[state="error"]',
                'QLabel[class="caption"]', "QTabBar::tab:selected",
                "QLineEdit:focus", 'QPlainTextEdit[class="log"]'):
        assert sel in qss, sel
    assert theme.COLORS["accent"] in qss, "QSS 应由令牌装配"


def test_icon_assets_present():
    for name in ("app", "tab-input", "tab-pes", "check", "save", "search",
                 "wand", "hammer", "folder", "chevron-down", "chevron-right"):
        assert (theme.ASSETS / f"{name}.svg").exists(), name


# ---------- 样式入口唯一性(防回漂: 页面代码不得再写 setStyleSheet) ----------

def test_no_local_stylesheets_in_pages():
    for rel in ("gui/pages/input_page.py", "gui/pages/pes_page.py",
                "gui/view3d/preview_view.py"):
        src = (ROOT / rel).read_text(encoding="utf-8")
        assert "setStyleSheet" not in src, f"{rel} 应走全局 QSS(theme.py)"


# ---------- 应用面(Qt offscreen) ----------

def test_apply_theme_offscreen():
    import inspect
    from PyQt6.QtWidgets import QApplication
    app = _app()
    theme.apply_theme(app)
    assert app.styleSheet(), "应用级 QSS 已装配"
    # QSS 装配后 app.style() 是 QStyleSheetStyle 代理(遮住底座), Fusion 底座
    # 以源镜像断言看守(与 test_dyn_ev_mirror 同法)
    assert app.style().metaObject().className() == "QStyleSheetStyle"
    assert 'setStyle("Fusion")' in inspect.getsource(theme.apply_theme)


def test_input_page_uses_theme_metrics():
    page = _page()
    m = page.layout().contentsMargins()
    assert (m.left(), m.top(), m.right(), m.bottom()) == (
        theme.SP_L, theme.SP_L, theme.SP_L, theme.SP_L)
    lm = page._scroll.widget().layout().contentsMargins()
    assert lm.left() == theme.SP_M


def test_section_header_styled_object_and_chevron():
    page = _page()
    btn = page._sections["run"]._btn
    assert btn.objectName() == "SectionHeader"
    btn.toggle()  # 折叠
    assert not btn.icon().isNull(), "折叠态 chevron 图标在场"
    btn.toggle()  # 展开
    assert not btn.icon().isNull(), "展开态 chevron 图标在场"


def test_primary_buttons():
    from PyQt6.QtWidgets import QPushButton
    page = _page()
    btns = [b for b in page.findChildren(QPushButton) if b.text().startswith("写出")]
    assert btns and btns[0].property("class") == "primary"


def test_verdict_state_hook():
    page = _page()
    page.run_validation()  # 未扫描+必填空 → issues 非空
    assert page._verdict.property("state") == "error"
    assert page._verdict_summary.property("state") == "error"
    assert "处问题" in page._verdict_summary.text()
    assert page._verdict.toPlainText().startswith("✗")  # 结果框文本照旧


def test_no_horizontal_scrollbar():
    """默认分栏份额下左栏不出横向滚动条(拉窄时的后备按窗口适配裁决保留)。"""
    from gui.main import MainWindow
    app = _app()
    theme.apply_theme(app)
    win = MainWindow()
    win.resize(1500, 950)
    win.show()
    for _ in range(4):
        app.processEvents()
    page = win.centralWidget().widget(0)
    need = page._scroll.widget().minimumSizeHint().width()
    have = page._scroll.viewport().width()
    assert need <= have, f"左栏内容最小宽 {need}px 超出可视宽 {have}px(应可折行/收缩)"
    assert not page._scroll.horizontalScrollBar().isVisible()
