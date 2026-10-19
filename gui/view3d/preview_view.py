# -*- coding: utf-8 -*-
"""preview_view.py - 起始散射预览视图(3D OpenGL 主 + 2D 投影后备)

设计:
  set_scene(PreviewScene) 同时喂给两后端:GLViewWidget 原子球/键/晶胞/箭头
  (世界坐标尺寸);QPainter 双面板投影(俯视 xy / 侧视 xz)——显卡/离屏无 GL
  时的后备,也是冒烟测试的可验证渲染路径。坐标笛卡尔 Å,z=表面法向。
  3D 侧浅底渲染(2026-10-13):背景/轴色走 theme 令牌(gl_bg/gl_axis),原子
  球色走 SCHEMES 方案表(工具条"配色"下拉),带尺度 xyz 轴 + 刻度角标。
  原子逐符号单色件(2026-10-13 真机修复):pyqtgraph 0.14 pxMode=False 的
  逐点色(VBO)路径真机不生效(原子全白);且散点默认 additive 混合在浅底
  gl_bg 上把球色加底饱和成白——每符号一件均匀色元组 + 标量尺寸 + opaque
  混合(常量顶点属性路径),渲出的即方案色本身。
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


def _slerp_arc(center, d_from, d_to, radius, n=24):
    """center 处从 d_from 方向到 d_to 方向的弧点列 (n,3)(最短角路径; 纯函数)。"""
    a = np.asarray(d_from, dtype=float)
    b = np.asarray(d_to, dtype=float)
    a /= np.linalg.norm(a)
    b /= np.linalg.norm(b)
    ang = float(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)))
    if ang < 1e-9:
        return np.zeros((0, 3))
    axis = np.cross(a, b)
    if np.linalg.norm(axis) < 1e-9:            # 反向: 任取正交轴
        axis = np.cross(a, [0.0, 0.0, 1.0])
        if np.linalg.norm(axis) < 1e-9:
            axis = np.cross(a, [0.0, 1.0, 0.0])
    axis /= np.linalg.norm(axis)
    out = []
    for t in np.linspace(0.0, ang, n):
        rot = (a * np.cos(t) + np.cross(axis, a) * np.sin(t))
        out.append(center + radius * rot)
    return np.array(out, dtype=float)


def guide_elements(scene):
    """入射角辅助线元素(2026-10-14): 锚=落点; 表面范式抬至衬底顶面+间隙
    (防被不透明表面原子遮挡), 气相(pair)保持 z=0 层。无束流/角色缺失返 None。

    元素: normal_seg(法线短枝)/dir_seg(弹丸方向短枝)/theta_arc(法线→弹丸
    方向弧)/chi_arc(+x→方位向水平弧)/chi_ray(方位射线), 以及 theta/chi 原值。
    """
    if scene.beam is None or not (0 <= scene.i_proj < len(scene.frag_coms)):
        return None
    d = scene.beam.drop
    z_lift = 0.0
    if scene.i_targ < 0:                        # 表面范式: 衬底=非弹丸原子
        a, b = scene.frag_slices[scene.i_proj]
        zs = [p[2] for k, p in enumerate(scene.positions) if not a <= k < b]
        if zs:
            z_lift = max(zs) + 0.4
    anchor = np.array([float(d[0]), float(d[1]), z_lift])
    v = np.asarray(scene.frag_coms[scene.i_proj], dtype=float) - anchor
    L = float(np.linalg.norm(v))
    if L < 1e-9:
        return None
    vhat = v / L
    r = min(max(0.35 * L, 1.2), 4.0)
    stub = 1.6 * r
    normal_seg = np.array([anchor, anchor + (0.0, 0.0, stub)])
    dir_seg = np.array([anchor, anchor + stub * vhat])
    theta_arc = _slerp_arc(anchor, (0.0, 0.0, 1.0), vhat, r)
    e_chi = np.array([np.cos(scene.beam.chi), np.sin(scene.beam.chi), 0.0])
    chi_arc = _slerp_arc(anchor, (1.0, 0.0, 0.0), e_chi, 0.8 * r)
    chi_ray = np.array([anchor, anchor + stub * e_chi])
    return {"anchor": anchor, "normal_seg": normal_seg, "dir_seg": dir_seg,
            "theta_arc": theta_arc, "chi_arc": chi_arc, "chi_ray": chi_ray,
            "theta": scene.beam.theta, "chi": scene.beam.chi}


def lattice_segments(cell, pos, cap: int = 16):
    """覆盖原子横向足迹的周期晶格网格(2026-10-14, 取代原点单胞框)。

    沿 a/b 晶向整数格线铺满足迹(floor/ceil), 底/顶两层 + 四角立柱;
    每方向跨度封顶 cap(防病态大 slab)。返回 (段数组 (m,2,3), z_lo, z_hi);
    横向晶矢奇异(零面积投影)返回空段。纯函数(numpy), 可离屏单测。
    """
    c = np.asarray(cell, dtype=float)
    M = c[:2, :2]                      # 行 = a1/a2 横向晶矢(c 为表面法向)
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    P = np.asarray(pos, dtype=float)
    if abs(det) < 1e-9 or P.size == 0:
        return np.zeros((0, 2, 3)), 0.0, 0.0
    fr = P[:, :2] @ np.linalg.inv(M)   # 分数坐标(行向量式: p = f·M)
    i0 = int(np.floor(fr[:, 0].min()))
    j0 = int(np.floor(fr[:, 1].min()))
    i1 = i0 + min(int(np.ceil(fr[:, 0].max())) - i0, cap)
    j1 = j0 + min(int(np.ceil(fr[:, 1].max())) - j0, cap)
    z_lo, z_hi = float(P[:, 2].min()), float(P[:, 2].max())
    flat = z_hi - z_lo < 1e-9

    def corner(i, j, z):
        xy = i * M[0] + j * M[1]
        return (float(xy[0]), float(xy[1]), float(z))

    zs = (z_lo,) if flat else (z_lo, z_hi)
    segs: list[tuple[tuple, tuple]] = []
    for z in zs:
        for j in range(j0, j1 + 1):
            segs.append((corner(i0, j, z), corner(i1, j, z)))
        for i in range(i0, i1 + 1):
            segs.append((corner(i, j0, z), corner(i, j1, z)))
    if not flat:
        for i in (i0, i1):
            for j in (j0, j1):
                segs.append((corner(i, j, z_lo), corner(i, j, z_hi)))
    return np.array(segs, dtype=float), z_lo, z_hi

# b 环圆周采样点数(构造处与单测的环件识别共用, 勿在两处各写字面量)
RING_SEGMENTS = 72


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
        # 原子逐符号单色件(2026-10-13 真机修复, 原子全白双因):(1) pyqtgraph
        # 0.14 pxMode=False 的逐点色(VBO)路径真机不生效——均匀色元组走
        # glVertexAttrib4f 常量顶点属性, 支持良好; (2) 散点默认 additive 混合
        # 在浅底 gl_bg 上把球色加底饱和成白(src·1+dst)——原子件改 opaque
        # (不混合+深度测试), 渲出的即方案色本身。同符号原子共用显示半径
        # → 每符号一件, size 为标量
        idx_by_sym: dict[str, list[int]] = {}
        for k, s in enumerate(scene.symbols):
            idx_by_sym.setdefault(s, []).append(k)
        for sym, idx in idx_by_sym.items():
            c = sch["colors"].get(sym, sch["fallback"])
            self._add(GLScatterPlotItem(
                pos=pos[idx], size=_display_radius(sym), pxMode=False,
                color=(c[0] / 255.0, c[1] / 255.0, c[2] / 255.0, 1.0),
                glOptions="opaque"))
        for i, j in _bonds(pos, scene):
            self._seg(pos[i], pos[j], (0.35, 0.35, 0.35, 1.0), 2.0)
        if scene.cell is not None:
            # 周期晶格网格覆盖表面片段足迹(2026-10-14; 单件批段, 色走 theme 令牌)
            # _rgb01 已含 alpha=1.0——不得再追加(5 元色组会让真 GL 的
            # glVertexAttrib4f 抛 TypeError, 逐项被吞→不可见, 2026-10-14 实机报障)
            segs, _, _ = lattice_segments(scene.cell, pos)
            if len(segs):
                rgba = _rgb01(theme.COLORS["gl_cell"])
                self._add(GLLinePlotItem(pos=segs.reshape(-1, 3), color=rgba,
                                         width=1.4, antialias=True, mode="lines"))
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
                th = np.linspace(0, 2 * np.pi, RING_SEGMENTS)
                ring = np.stack([ctr[0] + r * np.cos(th), ctr[1] + r * np.sin(th),
                                 np.full_like(th, ctr[2])], axis=1)
                self._add(GLLinePlotItem(pos=ring, color=(1.0, 0.55, 0.1, 0.8), width=2.0))
        # 入射角辅助线(2026-10-14): 法线/弹丸方向短枝 + θ 弧 + χ 弧与方位射线
        gd = guide_elements(scene)
        if gd is not None:
            rgba = _rgb01(theme.COLORS["gl_guide"])  # 4 元组(含 alpha=1), 勿再追加
            for seg in (gd["normal_seg"], gd["dir_seg"], gd["chi_ray"]):
                self._add(GLLinePlotItem(pos=np.asarray(seg, dtype=float),
                                         color=rgba, width=1.6, antialias=True,
                                         mode="lines"))
            for arc in (gd["theta_arc"], gd["chi_arc"]):
                if len(arc) >= 2:
                    self._add(GLLinePlotItem(pos=arc, color=rgba, width=1.6,
                                             antialias=True))
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
        txt = f"刻度 {step:g} Å"
        if scene.beam is not None:  # 入射角读数(度数直观; 弧度在表单/键表)
            txt += (f" · θ={np.degrees(scene.beam.theta):.1f}°"
                    f" · χ={np.degrees(scene.beam.chi):.1f}°")
        self._scale_label.setText(txt)
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
        # 晶格: 俯视画周期网格(覆盖足迹, 取景含网格端点), 侧视画胞框 z 范围
        grid2d = None
        if scene.cell is not None and tuple(ax) == (0, 1):
            grid2d, _, _ = lattice_segments(scene.cell, pos)
            if len(grid2d):
                extra.append(grid2d.reshape(-1, 3)[:, list(ax)])
        if scene.cell is not None and grid2d is None:
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

        # 晶格投影: 俯视=周期网格, 侧视=胞框
        if scene.cell is not None:
            pen = QPen(QColor(105, 120, 148))
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            if tuple(ax) == (0, 1):
                for s in (grid2d if grid2d is not None and len(grid2d)
                          else lattice_segments(scene.cell, pos)[0]):
                    p.drawLine(X(s[0][:2]), X(s[1][:2]))
            else:
                cs = _cell_corners(scene.cell)[:, list(ax)]
                for i, j in _CELL_EDGES:
                    p.drawLine(X(cs[i]), X(cs[j]))
        # 入射角辅助线(2026-10-14): 俯视=χ 弧+方位线, 侧视=法线/方向短枝+θ 弧
        gd = guide_elements(scene)
        if gd is not None:
            pen = QPen(QColor(theme.COLORS["gl_guide"]))
            p.setPen(pen)

            def draw_seg(s):
                p.drawLine(X((s[0][ax[0]], s[0][ax[1]])),
                           X((s[1][ax[0]], s[1][ax[1]])))

            def draw_poly(pts):
                for u, v in zip(pts[:-1], pts[1:]):
                    p.drawLine(X((u[ax[0]], u[ax[1]])), X((v[ax[0]], v[ax[1]])))

            if tuple(ax) == (0, 1):
                draw_seg(gd["chi_ray"])
                if len(gd["chi_arc"]) >= 2:
                    draw_poly(gd["chi_arc"])
            else:
                draw_seg(gd["normal_seg"])
                draw_seg(gd["dir_seg"])
                if len(gd["theta_arc"]) >= 2:
                    draw_poly(gd["theta_arc"])
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
        self._mode.addItem("3D 视图")
        self._mode.setItemData(0, "OpenGL 实时渲染,可拖拽旋转/缩放", Qt.ItemDataRole.ToolTipRole)
        self._mode.addItem("2D 投影")
        self._mode.setItemData(1, "无 GL 后备(QPainter 双面板投影)", Qt.ItemDataRole.ToolTipRole)
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
