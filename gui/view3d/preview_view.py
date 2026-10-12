# -*- coding: utf-8 -*-
"""preview_view.py - 起始散射预览视图(3D OpenGL 主 + 2D 投影后备)

设计:
  set_scene(PreviewScene) 同时喂给两后端:GLViewWidget 原子球/键/晶胞/箭头
  (世界坐标尺寸);QPainter 双面板投影(俯视 xy / 侧视 xz)——显卡/离屏无 GL
  时的后备,也是冒烟测试的可验证渲染路径。坐标笛卡尔 Å,z=表面法向。
"""

from __future__ import annotations

import numpy as np

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QVBoxLayout, QWidget
from pyqtgraph.opengl import GLLinePlotItem, GLScatterPlotItem, GLViewWidget

from ..core.preview import PreviewScene

# CPK 色 + 共价半径(成键判据;显示半径 0.45×共价+0.28)
CPK = {
    "H": (255, 255, 255), "D": (255, 255, 255),
    "C": (77, 77, 77), "N": (33, 77, 255), "O": (255, 33, 33),
    "F": (77, 230, 77), "Cl": (51, 204, 51), "Br": (153, 26, 26),
    "S": (230, 204, 26), "P": (255, 128, 0), "Au": (255, 209, 26),
    "Ag": (204, 204, 217), "Cu": (204, 115, 51), "Pt": (204, 204, 217),
    "Ni": (128, 128, 140), "Na": (102, 102, 255),
}
COVR = {
    "H": 0.31, "D": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Cl": 1.02, "Br": 1.20, "S": 1.05, "P": 1.07, "Au": 1.36,
    "Ag": 1.45, "Cu": 1.32, "Pt": 1.36, "Ni": 1.24, "Na": 1.66,
}


def _display_radius(sym: str) -> float:
    return 0.45 * COVR.get(sym, 0.8) + 0.28


def _bonds(pos: np.ndarray, scene: PreviewScene) -> list[tuple[int, int]]:
    out = []
    for a, b in scene.frag_slices:
        for i in range(a, b):
            for j in range(i + 1, b):
                d = float(np.linalg.norm(pos[i] - pos[j]))
                lim = COVR.get(scene.symbols[i], 0.8) + COVR.get(scene.symbols[j], 0.8) + 0.45
                if d < lim:
                    out.append((i, j))
    return out


def _cell_corners(cell) -> np.ndarray:
    c = np.asarray(cell, dtype=float)
    o = np.zeros(3)
    return np.array([
        o, c[0], c[0] + c[1], c[1],
        c[2], c[0] + c[2], c[0] + c[1] + c[2], c[1] + c[2],
    ])


_CELL_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0),
               (4, 5), (5, 6), (6, 7), (7, 4),
               (0, 4), (1, 5), (2, 6), (3, 7)]


class _Gl3D(QWidget):
    """GL 后端(薄封装)。"""

    def __init__(self):
        super().__init__()
        self.view = GLViewWidget()
        self.view.setCameraParams(distance=30.0, elevation=25, azimuth=-60)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.view)
        self._items = []

    def clear(self):
        for it in self._items:
            self.view.removeItem(it)
        self._items = []

    def _add(self, item):
        self.view.addItem(item)
        self._items.append(item)

    def _seg(self, p0, p1, color, width):
        self._add(GLLinePlotItem(pos=np.array([p0, p1], dtype=float),
                                 color=color, width=width, antialias=True))

    def set_scene(self, scene: PreviewScene):
        self.clear()
        pos = np.asarray(scene.positions, dtype=float)
        if pos.size == 0:
            return
        colors = [list(CPK.get(s, (190, 190, 195))) + [1.0] for s in scene.symbols]
        sizes = [_display_radius(s) for s in scene.symbols]
        self._add(GLScatterPlotItem(pos=pos, size=np.asarray(sizes),
                                    color=np.asarray(colors), pxMode=False))
        for i, j in _bonds(pos, scene):
            self._seg(pos[i], pos[j], (0.35, 0.35, 0.35, 1.0), 2.0)
        if scene.cell is not None:
            cs = _cell_corners(scene.cell)
            for i, j in _CELL_EDGES:
                self._seg(cs[i], cs[j], (0.55, 0.65, 0.75, 0.9), 1.2)
        if scene.aim_marker is not None:
            m = np.asarray(scene.aim_marker, dtype=float)
            d = 0.6
            self._seg(m + [-d, 0, 0], m + [d, 0, 0], (1.0, 0.55, 0.1, 1.0), 2.5)
            self._seg(m + [0, -d, 0], m + [0, d, 0], (1.0, 0.55, 0.1, 1.0), 2.5)
        if scene.b_max and scene.b_max > 0:
            r = float(scene.b_max)
            ctr = scene.frag_coms[1] if len(scene.frag_coms) > 1 else np.zeros(3)
            th = np.linspace(0, 2 * np.pi, 72)
            ring = np.stack([ctr[0] + r * np.cos(th), ctr[1] + r * np.sin(th),
                             np.full_like(th, ctr[2])], axis=1)
            self._add(GLLinePlotItem(pos=ring, color=(1.0, 0.55, 0.1, 0.8), width=2.0))
        for com, vel in zip(scene.frag_coms, scene.frag_vels):
            v = np.asarray(vel, dtype=float)
            n = float(np.linalg.norm(v))
            if n < 1e-9:
                continue
            u = v / n
            p0 = np.asarray(com, dtype=float)
            p1 = p0 + min(0.8 * n, 6.0) * u
            self._seg(p0, p1, (0.1, 0.5, 1.0, 1.0), 3.0)
            w = np.array([-u[1], u[0], 0.0])
            if np.linalg.norm(w) < 1e-9:
                w = np.array([1.0, 0.0, 0.0])
            w /= np.linalg.norm(w)
            for s in (+1, -1):
                self._seg(p1, p1 - 0.35 * u + s * 0.18 * w, (0.1, 0.5, 1.0, 1.0), 3.0)
        ctr = pos.mean(axis=0)
        span = float(np.linalg.norm(pos.max(axis=0) - pos.min(axis=0))) or 10.0
        self.view.opts["center"] = ctr
        self.view.opts["distance"] = max(2.2 * span, 12.0)


class _Projection2D(QWidget):
    """QPainter 双面板(俯视 xy / 侧视 xz)后备渲染。"""

    def __init__(self):
        super().__init__()
        self._scene: PreviewScene | None = None
        self.setMinimumHeight(320)

    def set_scene(self, scene: PreviewScene):
        self._scene = scene
        self.update()

    def clear(self):
        self._scene = None
        self.update()

    def paintEvent(self, ev):
        p = QPainter(self)
        p.fillRect(self.rect(), QColor(252, 252, 252))
        if self._scene is None or not self._scene.positions:
            p.setPen(QColor(120, 120, 120))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "(无预览数据)")
            p.end()
            return
        w = self.width() / 2 - 10
        h = self.height() - 30
        self._panel(p, QRectF(8, 22, w, h), "俯视 (x–y)", (0, 1), (0, 1))
        self._panel(p, QRectF(w + 18, 22, w, h), "侧视 (x–z)", (0, 2), (2, 1))
        p.end()

    def _panel(self, p, rect: QRectF, title: str, ax, flip):
        scene = self._scene
        pos = np.asarray(scene.positions, dtype=float)
        p.setPen(QColor(90, 90, 90))
        p.drawText(QPointF(rect.left(), rect.top() - 6), title)
        pts = pos[:, list(ax)]
        extra = []
        if scene.cell is not None:
            extra.append(_cell_corners(scene.cell)[:, list(ax)])
        if scene.aim_marker is not None:
            extra.append(np.asarray([scene.aim_marker])[..., list(ax)])
        all_pts = np.vstack([pts] + extra) if extra else pts
        lo = all_pts.min(axis=0) - 1.0
        hi = all_pts.max(axis=0) + 1.0
        span = np.maximum(hi - lo, 1e-6)
        scale = min(rect.width() / span[0], rect.height() / span[1])

        def X(q):
            u = (q - lo) / span
            return QPointF(rect.left() + u[0] * span[0] * scale,
                           rect.bottom() - u[1] * span[1] * scale)

        # 晶胞投影
        if scene.cell is not None:
            pen = QPen(QColor(140, 160, 180))
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            cs = _cell_corners(scene.cell)[:, list(ax)]
            for i, j in _CELL_EDGES:
                p.drawLine(X(cs[i]), X(cs[j]))
        # 瞄准点 + B_MAX 盘
        if scene.aim_marker is not None:
            m = np.asarray(scene.aim_marker, dtype=float)[list(ax)]
            pen = QPen(QColor(255, 140, 26), 2)
            p.setPen(pen)
            c = X(m)
            p.drawLine(c - QPointF(6, 0), c + QPointF(6, 0))
            p.drawLine(c - QPointF(0, 6), c + QPointF(0, 6))
        if scene.b_max and scene.b_max > 0:
            ctr = scene.frag_coms[1] if len(scene.frag_coms) > 1 else np.zeros(3)
            rpx = scene.b_max * scale
            p.setPen(QPen(QColor(255, 140, 26), 2))
            p.drawEllipse(X(np.asarray(ctr)[list(ax)]), rpx, rpx)
        # 键
        p.setPen(QPen(QColor(90, 90, 90), 2))
        for i, j in _bonds(pos, scene):
            p.drawLine(X(pts[i]), X(pts[j]))
        # 原子
        f = QFont("Arial", 7)
        for k, s in enumerate(scene.symbols):
            r = max(_display_radius(s) * scale, 3.5)
            col = CPK.get(s, (190, 190, 195))
            p.setPen(QPen(QColor(40, 40, 40)))
            p.setBrush(QBrush(QColor(*col)))
            c = X(pts[k])
            p.drawEllipse(c, r, r)
            p.setPen(QPen(QColor(20, 20, 20)))
            p.setFont(f)
            p.drawText(QPointF(c.x() - 5, c.y() - r - 2), s)
        # COM 速度箭头
        for com, vel in zip(scene.frag_coms, scene.frag_vels):
            v = np.asarray(vel, dtype=float)
            n = float(np.linalg.norm(v))
            if n < 1e-9:
                continue
            u = v / n
            p0 = np.asarray(com, dtype=float)
            p1 = p0 + min(0.8 * n, 6.0) * u
            a, b = X(p0[list(ax)]), X(p1[list(ax)])
            p.setPen(QPen(QColor(26, 128, 255), 3))
            p.drawLine(a, b)
            ang = np.arctan2(b.y() - a.y(), b.x() - a.x())
            for s in (0.45, -0.45):
                p.drawLine(b, QPointF(b.x() - 10 * np.cos(ang + s),
                                      b.y() - 10 * np.sin(ang + s)))


class PreviewView(QWidget):
    """预览面板:模式切换(3D OpenGL / 2D 投影)。"""

    MODE_3D = 0
    MODE_2D = 1

    def __init__(self, parent=None):
        super().__init__(parent)
        self._gl = _Gl3D()
        self._2d = _Projection2D()
        self._mode = QComboBox()
        self._mode.addItems(["3D(OpenGL,可拖拽旋转)", "2D 投影(无 GL 后备)"])
        self._mode.currentIndexChanged.connect(self._switch)
        self._stack_handle = [self._gl, self._2d]
        self._current = self._gl
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        top = QHBoxLayout()
        top.addWidget(QLabel("预览模式:"))
        top.addWidget(self._mode)
        top.addStretch(1)
        lay.addLayout(top)
        lay.addWidget(self._gl, 1)
        lay.addWidget(self._2d, 1)
        self._2d.hide()
        self._scene: PreviewScene | None = None

    def _switch(self, idx: int):
        self._current = self._stack_handle[idx]
        if idx == 0:
            self._2d.hide()
            self._gl.show()
        else:
            self._gl.hide()
            self._2d.show()

    def set_mode(self, idx: int):
        self._mode.setCurrentIndex(idx)

    def clear(self):
        self._scene = None
        self._gl.clear()
        self._2d.clear()

    def set_scene(self, scene: PreviewScene):
        self._scene = scene
        self._gl.set_scene(scene)
        self._2d.set_scene(scene)
