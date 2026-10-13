# -*- coding: utf-8 -*-
"""main.py - GUI 入口:python -m gui.main(或 python gui/main.py)。

两页:①输入文件生成(左表单/右 3D 起始散射预览)②自定义 PES 容器
检查(扫描→槽位生成→隔离树编译)。
"""

from __future__ import annotations

import sys
from pathlib import Path

# python gui/main.py 直跑时把仓库根放入 sys.path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from PyQt6.QtWidgets import QApplication, QMainWindow, QStatusBar, QTabWidget

from gui.pages.input_page import InputPage
from gui.pages.pes_page import PesPage


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VENUS2026 GUI — 输入生成 / 自定义 PES 容器检查")
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
        tabs.addTab(InputPage(), "输入文件生成")
        tabs.addTab(PesPage(), "PES 容器检查")
        self.setCentralWidget(tabs)
        self.setStatusBar(QStatusBar())
        self.statusBar().showMessage(
            "内部单位 eV/Å/amu/fs;预览为示意(平衡几何+代表值);编译验证走隔离副本,不动工作树"
        )


def main() -> int:
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
