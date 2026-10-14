# -*- coding: utf-8 -*-
"""main.py - GUI 入口:python -m gui.main(或 python gui/main.py)。

两页:①输入文件生成(左表单/右 3D 起始散射预览)②自定义 PES 容器
检查(扫描→程序补全→隔离树编译)。样式走 gui.theme(令牌化浅色主题,
Fusion 底座 + 全局 QSS——唯一样式入口)。
测试版发布开关 SHOW_PES_PAGE(2026-10-14 用户裁决: 先发布仅含输入
生成/校验的测试版):False 隐藏页2(代码保留, 单页隐藏页签条), True 恢复。
"""

from __future__ import annotations

import sys
from pathlib import Path

# python gui/main.py 直跑时把仓库根放入 sys.path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from PyQt6.QtWidgets import QApplication, QMainWindow, QStatusBar, QTabWidget

from gui import theme
from gui.pages.input_page import InputPage
from gui.pages.pes_page import PesPage

# 测试版发布开关:False = 仅页1(输入生成/校验);恢复两页改 True
SHOW_PES_PAGE = False


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VENUS2026 GUI — 输入生成"
                            + (" / 自定义 PES 容器检查" if SHOW_PES_PAGE else "(测试版)"))
        self.setWindowIcon(theme.icon("app"))
        # 初始尺寸钳制在屏幕可用几何内(任务栏除外——下沿不得被任务栏挡住)
        screen = QApplication.primaryScreen()
        avail = screen.availableGeometry() if screen else None
        w, h = 1500, 950
        if avail is not None:
            w, h = min(w, avail.width()), min(h, avail.height())
        self.resize(w, h)
        self.setMinimumSize(min(900, w), min(700, h))
        if avail is not None:
            self.move(avail.topLeft())
        tabs = QTabWidget()
        tabs.addTab(InputPage(), theme.icon("tab-input"), "输入文件生成")
        if SHOW_PES_PAGE:
            tabs.addTab(PesPage(), theme.icon("tab-pes"), "PES 容器检查")
        else:
            tabs.tabBar().setVisible(False)  # 单页发布:隐藏页签条
        self.setCentralWidget(tabs)
        self.setStatusBar(QStatusBar())
        self.statusBar().showMessage(
            "内部单位 eV/Å/amu/fs;预览为示意"
            + (";编译验证走隔离副本" if SHOW_PES_PAGE else ""))


def main() -> int:
    app = QApplication(sys.argv)
    theme.apply_theme(app)  # Fusion + 调色板 + 字体 + 全局 QSS(先于建窗)
    app.setWindowIcon(theme.icon("app"))
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
