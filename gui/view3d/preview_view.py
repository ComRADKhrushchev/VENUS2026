# -*- coding: utf-8 -*-
"""preview_view.py - 起始散射预览视图(3D OpenGL 主 + 2D 投影后备)

设计:
  set_scene(PreviewScene) 同时喂给两后端:GLViewWidget 原子球/键/晶胞/箭头
  (世界坐标尺寸);QPainter 双面板投影(俯视 xy / 侧视 xz)——显卡/离屏无 GL
  时的后备,也是冒烟测试的可验证渲染路径。坐标笛卡尔 Å,z=表面法向。
  3D 侧浅底渲染(2026-10-13):背景/轴色走 theme 令牌(gl_bg/gl_axis),原子
  球色走 SCHEMES 方案表(工具条"配色"下拉),带尺度 xyz 轴 + 刻度角标。
  取景恒以靶为中心(2026-10-13 修):b 环心/相机中心/2D 取景窗对称化都绑
  preview.ring_center/view_center(靶 COM/落点;角色未解析回退全原子均值)。
"""

from __future__ import annotations

import math

import numpy as np

from PyQt6.QtCore import QEvent, QPointF, QRectF, Qt
from PyQt6.QtGui import QBrush, QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QVBoxLayout, QWidget
from pyqtgraph.opengl import GLLinePlotItem, GLScatterPlotItem, GLViewWidget

from .. import theme
from ..core.preview import PreviewScene, ring_center, view_center

# 经典 CPK 色(2D 投影后备查表,画布冻结不随方案变)+ 共价半径
# (成键判据;显示半径 0.45×共价+0.28)
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

# ---------- 原子配色方案表(3D 球色查表集中一处;两方案同键域 = CPK 键集 + 未知回退) ----------
# cpk = CPK 经典浅底适配版,对经典 CPK 的已记录偏离(浅底可辨):
#   H/D 白 (255,255,255)→浅灰;Ag/Pt 银白 (204,204,217)→加深;Au 略压暗;
#   未知回退 (190,190,195)→中灰。2D 投影后备不走本表(画布冻结,仍用上方经典 CPK)。
# okabe = 元素高区分度方案(Okabe-Ito 色盲安全系衍生):8 基色 + 同族明暗两档;
#   黄 #F0E442 压暗为 (184,146,0)、黑软化 (32,32,32) 以适配浅底;D 与 H 同色
#   (同位素,沿 CPK 惯例);P/Au/Cu/Ni/Na 为各基色族的暗档衍生。
SCHEMES: dict[str, dict] = {
    "cpk": {
        "label": "CPK(经典·浅底适配)",
        "colors": {
            **CPK,
            "H": (196, 196, 196), "D": (196, 196, 196),
            "Ag": (172, 175, 190), "Pt": (172, 175, 190),
            "Au": (240, 190, 20),
        },
        "fallback": (140, 140, 150),
    },
    "okabe": {
        "label": "高区分(色盲安全)",
        "colors": {
            "H": (86, 180, 233), "D": (86, 180, 233),   # OI 天蓝
            "C": (32, 32, 32),                           # OI 黑(软化)
            "N": (0, 114, 178),                          # OI 蓝
            "O": (213, 94, 0),                           # OI 朱红
            "F": (230, 159, 0),                          # OI 橙
            "Cl": (0, 158, 115),                         # OI 蓝绿
            "Br": (204, 121, 167),                       # OI 红紫
            "S": (184, 146, 0),                          # OI 黄(压暗)
            "P": (94, 60, 153),                          # 红紫族暗档
            "Au": (146, 101, 32),                        # 黄族暗档(金褐)
            "Ag": (148, 153, 161),                       # 中性钢灰
            "Cu": (140, 62, 22),                         # 朱红族暗档(赭)
            "Pt": (96, 96, 106),                         # 深灰
            "Ni": (0, 97, 108),                          # 蓝绿族暗档(深青)
            "Na": (0, 48, 85),                           # 蓝族暗档(藏青)
        },
        "fallback": (105, 105, 115),
    },
}
DEFAULT_SCHEME = "cpk"


def nice_step(span: float, target: int = 5) -> float:
    """刻度整洁步长:1/2/5×10^n 中取不小于 span/target 的最小者。

    刻度段数由此落在约 (target/2, target];跨度非正退化为 1.0(纯函数,不依赖 Qt)。
    """
    if span <= 0:
        return 1.0
    raw = span / max(target, 1)
    mag = 10.0 ** math.floor(math.log10(raw))
    for m in (1.0, 2.0, 5.0):
        if m * mag >= raw * (1.0 - 1e-12):
            return m * mag
    return 10.0 * mag


def symmetric_window(lo, hi, center):
    """2D 取景窗对称化(纯函数):各轴半宽 = max(center−lo, hi−center)。

    取景中心(view_center:靶 COM/落点)由此恒落面板几何中心,原包围盒仍含于窗内。
    """
    c = np.asarray(center, dtype=float)
    half = np.maximum(c - lo, hi - c)
    return c - half, c + half


def _rgb01(hex_color: str) -> tuple[float, float, float, float]:
    """主题色令牌(HEX)→ GL 顶点色 (r,g,b,1)。"""
    c = QColor(hex_color)
    return (c.red() / 255.0, c.green() / 255.0, c.blue() / 255.0, 1.0)


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
    """GL 后端(薄封装):浅底令牌背景 + 方案球色 + 带尺度 xyz 轴。"""

    def __init__(self):
        super().__init__()
        self._scheme = DEFAULT_SCHEME
        self.view = GLViewWidget()
        self.view.setCameraParams(distance=30.0, elevation=25, azimuth=-60)
        self.view.setBackgroundColor(QColor(theme.COLORS["gl_bg"]))
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.view)
        self._items = []
        self._axis_items = []
        # 刻度角标:Qt 控件浮于视口右上(QSS 经 objectName=GlScaleLabel)
        self._scale_label = QLabel(self.view)
        self._scale_label.setObjectName("GlScaleLabel")
        self._scale_label.hide()
        self.view.installEventFilter(self)

    def eventFilter(self, obj, ev):
        if obj is self.view and ev.type() == QEvent.Type.Resize:
            self._reposition_scale()
        return super().eventFilter(obj, ev)

    def _reposition_scale(self):
        self._scale_label.adjustSize()
        self._scale_label.move(self.view.width() - self._scale_label.width() - 6, 6)

    def set_scheme(self, key: str):
        """记住原子配色方案(重绘由上层重喂场景触发)。"""
        if key in SCHEMES:
            self._scheme = key

    def clear(self):
        for it in self._items:
            self.view.removeItem(it)
        self._items = []
        self._axis_items = []
        self._scale_label.hide()

    def _add(self, item):
        self.view.addItem(item)
        self._items.append(item)

    def _seg(self, p0, p1, color, width):
        self._add(GLLinePlotItem(pos=np.array([p0, p1], dtype=float),
                                 color=color, width=width, antialias=True))

    def _axis_add(self, item):
        self._add(item)
        self._axis_items.append(item)

    def _draw_axes(self, lo, hi, step: float):
        """带尺度 xyz 轴。
        ! 1. 过场景原点三轴,各向覆盖 [min(包围盒, 0), max(包围盒, 0)]
        ! 2. 每轴按整洁步长画刻度短线(mode='lines' 批段);原点三轴交汇处不画
        ! 3. 轴/刻度色 = theme 令牌 gl_axis(浅底可辨)
        """
        col = _rgb01(theme.COLORS["gl_axis"])
        tick = 0.22 * step
        for k in range(3):
            a_lo, a_hi = min(float(lo[k]), 0.0), max(float(hi[k]), 0.0)
            if a_hi - a_lo < 1e-9:
                continue  # 该向场景零跨度,无轴可标
            p0, p1 = [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]
            p0[k], p1[k] = a_lo, a_hi
            self._axis_add(GLLinePlotItem(pos=np.array([p0, p1], dtype=float),
                                          color=col, width=1.6, antialias=True))
            perp = [0.0, 0.0, 0.0]
            perp[(k + 1) % 3] = tick
            segs = []
            for n in range(math.ceil(a_lo / step - 1e-9),
                           math.floor(a_hi / step + 1e-9) + 1):
                if n == 0:
                    continue
                q = [0.0, 0.0, 0.0]
                q[k] = n * step
                segs.append((q, (q[0] + perp[0], q[1] + perp[1], q[2] + perp[2])))
            if segs:
                self._axis_add(GLLinePlotItem(
                    pos=np.asarray(segs, dtype=float).reshape(-1, 3),
                    color=col, width=1.2, antialias=True, mode="lines"))

    def set_scene(self, scene: PreviewScene):
        self.clear()
        pos = np.asarray(scene.positions, dtype=float)
        if pos.size == 0:
            return
        sch = SCHEMES[self._scheme]
        # GL 顶点色须 0.0-1.0(表存 0-255;旧码直传 0-255 会被 GL 钳到全白)
        colors = [[c / 255.0 for c in sch["colors"].get(s, sch["fallback"])] + [1.0]
                  for s in scene.symbols]
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
            rc = ring_center(scene)  # b 环心绑靶 COM(角色缺失不画环)
            if rc is not None:
                r = float(scene.b_max)
                ctr = np.asarray(rc, dtype=float)
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
        # 带尺度 xyz 轴 + 刻度角标(步长 = 场景跨度的整洁值;lo/hi 复用到相机取景)
        lo, hi = pos.min(axis=0), pos.max(axis=0)
        span = float(np.linalg.norm(hi - lo)) or 10.0
        step = nice_step(span)
        self._draw_axes(lo, hi, step)
        self._scale_label.setText(f"刻度 {step:g} Å")
        self._scale_label.show()
        self._reposition_scale()
        # 相机中心绑取景中心(靶 COM/落点;未解析回退全原子均值)——靶恒居中,
        # 不再被 R_SEP 远处的弹丸拉偏(2026-10-13 修)
        ctr = np.asarray(view_center(scene), dtype=float)
        # opts['center'] 须为 QVector3D——pyqtgraph viewMatrix 对其调 .x()/.y()/.z(),
        # ndarray 会让每帧 paintGL 抛 AttributeError(空白视图)
        from pyqtgraph.Qt import QtGui
        self.view.opts["center"] = QtGui.QVector3D(float(ctr[0]), float(ctr[1]), float(ctr[2]))
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
        # 取景窗对称化:取景中心(靶 COM/落点)恒落面板几何中心(2026-10-13 修)
        vc = np.asarray(view_center(scene), dtype=float)[list(ax)]
        lo, hi = symmetric_window(all_pts.min(axis=0) - 1.0, all_pts.max(axis=0) + 1.0, vc)
        span = np.maximum(hi - lo, 1e-6)
        scale = min(rect.width() / span[0], rect.height() / span[1])
        # 内容盒居中落面板(松弛轴两侧均分留白,不再贴左下锚)
        off = (rect.width() - span[0] * scale) / 2, (rect.height() - span[1] * scale) / 2

        def X(q):
            u = (q - lo) / span
            return QPointF(rect.left() + off[0] + u[0] * span[0] * scale,
                           rect.bottom() - off[1] - u[1] * span[1] * scale)

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
            rc = ring_center(scene)  # b 环心绑靶 COM(角色缺失不画环)
            if rc is not None:
                rpx = scene.b_max * scale
                p.setPen(QPen(QColor(255, 140, 26), 2))
                p.drawEllipse(X(np.asarray(rc, dtype=float)[list(ax)]), rpx, rpx)
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
    """预览面板:模式切换(3D OpenGL / 2D 投影)+ 3D 原子配色方案。"""

    MODE_3D = 0
    MODE_2D = 1

    def __init__(self, parent=None):
        super().__init__(parent)
        self._gl = _Gl3D()
        self._2d = _Projection2D()
        self._mode = QComboBox()
        self._mode.addItems(["3D(OpenGL,可拖拽旋转)", "2D 投影(无 GL 后备)"])
        self._mode.currentIndexChanged.connect(self._switch)
        self._scheme = QComboBox()
        for key, sch in SCHEMES.items():
            self._scheme.addItem(sch["label"], key)  # itemData 承载方案键
        self._scheme.currentIndexChanged.connect(self._apply_scheme)
        self._stack_handle = [self._gl, self._2d]
        self._current = self._gl
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        top = QHBoxLayout()
        mode_lab = QLabel("预览模式:")
        mode_lab.setProperty("class", "caption")
        top.addWidget(mode_lab)
        top.addWidget(self._mode)
        top.addSpacing(theme.SP_L)
        scheme_lab = QLabel("配色:")
        scheme_lab.setProperty("class", "caption")
        top.addWidget(scheme_lab)
        top.addWidget(self._scheme)
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

    def _apply_scheme(self, idx: int):
        key = self._scheme.itemData(idx)
        if key is None:
            return
        self._gl.set_scheme(key)
        if self._scene is not None:
            self._gl.set_scene(self._scene)

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
