# -*- coding: utf-8 -*-
"""theme.py - GUI 主题中枢: 色彩/间距/排版令牌 + 全局 QSS(Fusion 底座)。

设计:
  唯一样式入口(2026-10-13 用户裁决: 方案 A 令牌化浅色主题)——页面代码不写
  setStyleSheet, 局部差异走 objectName(如 #SectionHeader)或 dynamic property
  (class="primary|caption|card|log", state="ok|error"), 全部在 build_qss()
  收口; 令牌为纯量常量(本模块顶层不 import Qt, 无 Qt 环境可读可测);
  apply_theme() 一次装配 Fusion + 字体 + 调色板 + QSS。
"""

from __future__ import annotations

from pathlib import Path

ASSETS = Path(__file__).resolve().parent / "assets" / "icons"

# ---------- 色彩令牌(角色化命名; 调主题只动这里) ----------
COLORS: dict[str, str] = {
    "accent": "#2563EB",          # 交互强调: 主按钮/聚焦/活动页签
    "accent_hover": "#1D4ED8",
    "accent_pressed": "#1E40AF",
    "accent_soft": "#EFF4FE",     # 强调浅底: 悬停/选中/表头
    "surface": "#F3F4F6",         # 窗底
    "card": "#FFFFFF",            # 卡面/输入面
    "band": "#E9EDF4",            # 分组标题色带
    "band_hover": "#DFE6F2",
    "border": "#D1D5DB",
    "gridline": "#E5E7EB",
    "head": "#F8FAFC",            # 表头/次级浅底
    "text_primary": "#1F2937",
    "text_secondary": "#4B5563",
    "text_muted": "#6B7280",
    "success": "#16A34A",
    "success_soft": "#F0FDF4",
    "error": "#DC2626",
    "error_soft": "#FEF2F2",
    "tooltip_bg": "#1F2937",
    "tooltip_fg": "#F9FAFB",
    "scroll_handle": "#C6CCD8",
    "scroll_handle_hover": "#AAB2C2",
    "gl_bg": "#D7E1EF",           # 3D 预览视口底: 浅而不近白(白系原子可辨)
    "gl_axis": "#44546A",          # 3D 预览 xyz 轴/刻度线: 浅底可辨的深蓝灰
    "gl_cell": "#5A6B85",          # 3D 预览周期晶格网格: 区别于轴线, 浅底可辨
    "gl_guide": "#8B5E8C",         # 3D 预览入射角辅助线(θ/χ 弧): 梅紫, 与三锚色区分
}

# 间距标尺(px): 所有布局间距取自此处, 页面不写裸字面量
SP_S, SP_M, SP_L, SP_XL = 4, 8, 12, 16

# 排版(pt, 随系统 DPI 缩放): 正文 9 / 组头 10 / 说明 8, 三档为限
FONT_FAMILY = "Microsoft YaHei UI"
FONT_BODY_PT = 9
FONT_SECTION_PT = 10
FONT_CAPTION_PT = 8


def build_qss() -> str:
    """全局样式表: 纯字符串装配, 不依赖 Qt。"""
    c = COLORS
    return f"""
/* ---- 基面 ---- */
QWidget {{ color: {c['text_primary']}; font-size: {FONT_BODY_PT}pt; }}
QScrollArea, QScrollArea > QWidget {{ border: none; background: transparent; }}

/* ---- 分组标题条(_Section 的 QToolButton, objectName=SectionHeader) ---- */
QToolButton#SectionHeader {{
    background: {c['band']};
    border: none;
    border-left: 3px solid {c['accent']};
    border-radius: 4px;
    padding: 5px {SP_M}px;
    font-size: {FONT_SECTION_PT}pt;
    font-weight: 600;
    color: {c['text_primary']};
    text-align: left;
}}
QToolButton#SectionHeader:hover {{ background: {c['band_hover']}; }}
QToolButton#SectionHeader:checked, QToolButton#SectionHeader:unchecked {{ background: {c['band']}; }}
QToolButton#SectionHeader:pressed {{ background: {c['band_hover']}; }}

/* ---- 输入控件 ---- */
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
    background: {c['card']};
    border: 1px solid {c['border']};
    border-radius: 4px;
    padding: 2px 6px;
    selection-background-color: {c['accent']};
    selection-color: #FFFFFF;
}}
QLineEdit:hover, QComboBox:hover, QSpinBox:hover, QDoubleSpinBox:hover {{ border-color: {c['accent']}; }}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 1px solid {c['accent']};
    background: {c['accent_soft']};
}}
QLineEdit:disabled, QComboBox:disabled {{ color: {c['text_muted']}; border-color: {c['gridline']}; }}
QComboBox::drop-down {{ border: none; background: transparent; width: 18px; }}
QComboBox QAbstractItemView {{
    background: {c['card']};
    border: 1px solid {c['border']};
    selection-background-color: {c['accent_soft']};
    selection-color: {c['text_primary']};
}}

/* ---- 按钮: 默认次级描边, class=primary 实心强调 ---- */
QPushButton {{
    background: {c['card']};
    border: 1px solid {c['border']};
    border-radius: 4px;
    padding: 4px 14px;
    color: {c['text_primary']};
}}
QPushButton:hover {{ background: {c['accent_soft']}; border-color: {c['accent']}; color: {c['accent']}; }}
QPushButton:pressed {{ background: {c['band_hover']}; }}
QPushButton:focus {{ border: 1px solid {c['accent']}; }}
QPushButton:disabled {{ color: {c['text_muted']}; border-color: {c['gridline']}; background: {c['head']}; }}
QPushButton[class="primary"] {{
    background: {c['accent']};
    border: 1px solid {c['accent']};
    color: #FFFFFF;
    font-weight: 600;
}}
QPushButton[class="primary"]:hover {{ background: {c['accent_hover']}; border-color: {c['accent_hover']}; color: #FFFFFF; }}
QPushButton[class="primary"]:pressed {{ background: {c['accent_pressed']}; border-color: {c['accent_pressed']}; }}
QPushButton[class="primary"]:focus {{ border: 1px solid {c['accent_pressed']}; }}
QPushButton[class="primary"]:disabled {{ background: {c['band']}; border-color: {c['band']}; color: {c['text_muted']}; }}

/* ---- 页签: 活动页下缘强调线 ---- */
QTabWidget::pane {{ border: none; }}
QTabBar {{ background: transparent; }}
QTabBar::tab {{
    background: transparent;
    color: {c['text_muted']};
    padding: 6px 16px;
    border: none;
    border-bottom: 2px solid transparent;
    margin-right: {SP_S}px;
}}
QTabBar::tab:hover {{ color: {c['text_primary']}; }}
QTabBar::tab:selected {{ color: {c['accent']}; border-bottom: 2px solid {c['accent']}; font-weight: 600; }}

/* ---- 卡片容器 ---- */
QGroupBox {{
    background: {c['card']};
    border: 1px solid {c['border']};
    border-radius: 6px;
    margin-top: {SP_M}px;
    font-weight: 600;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: {SP_M}px;
    padding: 0 {SP_S}px;
    color: {c['text_secondary']};
    background: {c['card']};
}}
QFrame[class="card"] {{ background: {c['card']}; border: 1px solid {c['border']}; border-radius: 6px; }}

/* ---- 表格 ---- */
QTableView {{
    background: {c['card']};
    border: 1px solid {c['border']};
    border-radius: 4px;
    gridline-color: {c['gridline']};
    selection-background-color: {c['accent_soft']};
    selection-color: {c['text_primary']};
}}
QHeaderView::section {{
    background: {c['head']};
    border: none;
    border-bottom: 1px solid {c['border']};
    border-right: 1px solid {c['gridline']};
    padding: 4px {SP_M}px;
    color: {c['text_secondary']};
    font-weight: 600;
}}
QTableCornerButton::section {{ background: {c['head']}; border: none; }}

/* ---- 文本域: 日志等宽(class=log); 校验结果按 state 语义色 ---- */
QPlainTextEdit {{
    background: {c['card']};
    border: 1px solid {c['border']};
    border-radius: 4px;
    padding: {SP_S}px;
    selection-background-color: {c['accent']};
    selection-color: #FFFFFF;
}}
QPlainTextEdit:focus {{ border: 1px solid {c['accent']}; }}
QPlainTextEdit[class="log"] {{ font-family: "Consolas", monospace; font-size: {FONT_BODY_PT}pt; }}
QPlainTextEdit[state="ok"] {{ border-color: {c['success']}; background: {c['success_soft']}; }}
QPlainTextEdit[state="error"] {{ border-color: {c['error']}; background: {c['error_soft']}; }}

/* ---- 标签: 说明文字(class=caption); 校验摘要行按 state 语义色 ---- */
QLabel[class="caption"] {{ color: {c['text_muted']}; font-size: {FONT_CAPTION_PT}pt; }}
QLabel[class="verdict-summary"] {{ font-weight: 600; }}
QLabel[class="verdict-summary"][state="ok"] {{ color: {c['success']}; }}
QLabel[class="verdict-summary"][state="error"] {{ color: {c['error']}; }}

/* ---- 滚动条(细杆) / 分隔手柄 / 状态栏 / 气泡 ---- */
QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px 0; }}
QScrollBar::handle:vertical {{ background: {c['scroll_handle']}; border-radius: 3px; min-height: 28px; }}
QScrollBar::handle:vertical:hover {{ background: {c['scroll_handle_hover']}; }}
QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 0 2px; }}
QScrollBar::handle:horizontal {{ background: {c['scroll_handle']}; border-radius: 3px; min-width: 28px; }}
QScrollBar::handle:horizontal:hover {{ background: {c['scroll_handle_hover']}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}

QSplitter::handle:horizontal {{ width: 6px; background: transparent; }}
QSplitter::handle:horizontal:hover {{ background: {c['accent_soft']}; border-radius: 3px; }}

QStatusBar {{ background: transparent; border-top: 1px solid {c['border']}; color: {c['text_muted']}; }}

QToolTip {{
    background: {c['tooltip_bg']};
    color: {c['tooltip_fg']};
    border: none;
    padding: 4px {SP_M}px;
    border-radius: 3px;
}}

/* ---- 3D 预览: 刻度角标(GL 视口右上, Qt 控件不依赖 GL 文本) ---- */
QLabel#GlScaleLabel {{
    background: {c['card']};
    color: {c['text_secondary']};
    border: 1px solid {c['border']};
    border-radius: 3px;
    padding: 1px 6px;
    font-size: {FONT_CAPTION_PT}pt;
}}
"""


def icon(name: str):
    """按名取 assets/icons/ 下的 SVG 图标; 文件缺失返回空 QIcon(不崩)。"""
    from PyQt6.QtGui import QIcon
    p = ASSETS / f"{name}.svg"
    return QIcon(str(p)) if p.exists() else QIcon()


def apply_theme(app) -> None:
    """一次装配: Fusion 底座 + 调色板 + 应用字体 + 全局 QSS。"""
    from PyQt6.QtGui import QColor, QFont, QPalette

    app.setStyle("Fusion")
    c = COLORS
    pal = QPalette()
    pal.setColor(QPalette.ColorRole.Window, QColor(c["surface"]))
    pal.setColor(QPalette.ColorRole.Base, QColor(c["card"]))
    pal.setColor(QPalette.ColorRole.AlternateBase, QColor(c["head"]))
    pal.setColor(QPalette.ColorRole.Text, QColor(c["text_primary"]))
    pal.setColor(QPalette.ColorRole.WindowText, QColor(c["text_primary"]))
    pal.setColor(QPalette.ColorRole.Button, QColor(c["card"]))
    pal.setColor(QPalette.ColorRole.ButtonText, QColor(c["text_primary"]))
    pal.setColor(QPalette.ColorRole.Highlight, QColor(c["accent"]))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor("#FFFFFF"))
    pal.setColor(QPalette.ColorRole.ToolTipBase, QColor(c["tooltip_bg"]))
    pal.setColor(QPalette.ColorRole.ToolTipText, QColor(c["tooltip_fg"]))
    pal.setColor(QPalette.ColorRole.PlaceholderText, QColor(c["text_muted"]))
    pal.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, QColor(c["text_muted"]))
    pal.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, QColor(c["text_muted"]))
    app.setPalette(pal)
    app.setFont(QFont(FONT_FAMILY, FONT_BODY_PT))
    app.setStyleSheet(build_qss())
