# -*- coding: utf-8 -*-
"""Task 1(2026-10-13 3D 预览渲染增强)单测:

  - nice_step: 1/2/5×10^n 整洁刻度步长(纯函数, 不需 QApplication);
  - SCHEMES: 两方案同键域(现 CPK 键集 + 未知元素回退), CPK 浅底适配偏离已锁,
    gl_bg/gl_axis 令牌浅而不近白/浅底可辨;
  - _Gl3D: offscreen 构造 + set_scene 不崩, opts 契约(center Vector /
    distance 数值)不变, 背景落 gl_bg 令牌, 轴件在场, 刻度角标文本
    = nice_step 派生值;
  - PreviewView: 配色下拉默认 CPK, 切换后 GL 球色重绘。
"""
import math
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6")
pytest.importorskip("pyqtgraph")

from gui import theme
from gui.view3d.preview_view import (
    DEFAULT_SCHEME,
    SCHEMES,
    PreviewView,
    _Gl3D,
    nice_step,
)

CPK_KEYS = {"H", "D", "C", "N", "O", "F", "Cl", "Br", "S", "P",
            "Au", "Ag", "Cu", "Pt", "Ni", "Na"}


def _rel_lum(rgb) -> float:
    return (0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]) / 255.0


def _hex_rgb(h: str):
    return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))


# ---------- 纯函数: 刻度整洁值(1/2/5×10^n) ----------

@pytest.mark.parametrize("span,expect", [
    (10.0, 2.0), (20.0, 5.0), (7.0, 2.0), (3.0, 1.0),
    (100.0, 20.0), (250.0, 50.0),
    (1.0, 0.2), (0.5, 0.1), (0.23, 0.05),
    (0.0, 1.0), (-4.0, 1.0),  # 退化跨度: 不崩, 给 1.0
])
def test_nice_step_values(span, expect):
    assert nice_step(span) == pytest.approx(expect)


@pytest.mark.parametrize("span", [0.07, 0.9, 2.71, 13.7, 99.9, 1234.5])
def test_nice_step_is_one_two_five(span):
    s = nice_step(span)
    assert s > 0
    mant = s / 10.0 ** math.floor(math.log10(s))
    assert any(mant == pytest.approx(m0) for m0 in (1.0, 2.0, 5.0)), (span, s)


@pytest.mark.parametrize("span", [0.07, 0.9, 2.71, 13.7, 99.9, 1234.5])
def test_nice_step_tick_density(span):
    """步长不细于 2.5 段/跨度, 不粗于 5 段/跨度(刻度可读且不过密)。"""
    d = span / nice_step(span)
    assert 2.4 <= d <= 5.0, (span, d)


# ---------- 配色方案表: 键域完备 / 令牌可辨性 ----------

def test_schemes_key_domain():
    """两方案同键域 = 现 CPK 键集; 每键 3-tuple int 0-255; 各带未知回退。"""
    assert len(SCHEMES) >= 2
    assert DEFAULT_SCHEME == "cpk" and DEFAULT_SCHEME in SCHEMES
    domains = []
    for key, sch in SCHEMES.items():
        assert set(sch["colors"]) == CPK_KEYS, key
        domains.append(set(sch["colors"]))
        fb = sch["fallback"]
        assert len(fb) == 3 and all(isinstance(v, int) and 0 <= v <= 255 for v in fb), key
        for sym, col in sch["colors"].items():
            assert len(col) == 3, (key, sym)
            assert all(isinstance(v, int) and 0 <= v <= 255 for v in col), (key, sym)
        assert sch["label"], key
    assert all(d == domains[0] for d in domains), "两方案须同键域"


def test_cpk_scheme_light_adapted():
    """CPK 浅底适配: H/D 白→浅灰(经典 CPK 的已记录偏离), 亮度不贴底色。"""
    col = SCHEMES["cpk"]["colors"]
    assert col["H"] != (255, 255, 255), "白 H 在浅底不可辨, 须适配"
    assert col["H"] == col["D"]
    assert _rel_lum(col["H"]) <= 0.85, col["H"]


def test_gl_tokens_light_but_not_white():
    """gl_bg 浅而不近白(白系原子可辨); gl_axis 足够深(浅底可辨)。"""
    bg = _hex_rgb(theme.COLORS["gl_bg"])
    assert 0.5 <= _rel_lum(bg) < 0.90, bg
    ax = _hex_rgb(theme.COLORS["gl_axis"])
    assert _rel_lum(ax) <= 0.5, ax


# ---------- GL 构造冒烟(offscreen 不崩, 契约不变) ----------

def _scene():
    from gui.core.preview import PreviewScene
    return PreviewScene(
        symbols=["H", "O", "H"],
        masses=[1.0, 16.0, 1.0],
        positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 1.0], [0.0, 3.5, 4.2]],
        frag_slices=[(0, 3)], frag_names=["tri"],
        frag_coms=[[0.833, 1.167, 1.733]], frag_vels=[[0.1, 0.0, 0.0]],
    )


def test_gl_offscreen_smoke():
    """offscreen 构造+set_scene 不崩; opts 类型契约不变; 背景/轴/角标落位。"""
    from PyQt6.QtGui import QColor
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    gl = _Gl3D()
    # 背景: 主题令牌 gl_bg → GLViewWidget opts['bgcolor'](RGBF 元组)
    assert tuple(gl.view.opts["bgcolor"]) == QColor(theme.COLORS["gl_bg"]).getRgbF()
    for key in SCHEMES:
        gl.set_scheme(key)
        gl.set_scene(_scene())
        center = gl.view.opts["center"]
        assert all(callable(getattr(center, a, None)) for a in ("x", "y", "z")), \
            f"opts['center'] 为 {type(center)}(方案 {key})"
        assert isinstance(gl.view.opts["distance"], (int, float))
    # 带尺度轴: 三主线 + 至少一组刻度短线
    assert len(gl._axis_items) >= 4
    # 刻度角标: 跨度 |(2.5,3.5,4.2)|≈6.01 → 整洁步长 2
    assert gl._scale_label.text() == "刻度 2 Å"
    assert not gl._scale_label.isHidden()
    gl.clear()
    assert gl._scale_label.isHidden(), "清场后角标应隐藏"
    assert gl._axis_items == []


def test_scheme_combo_wiring():
    """PreviewView 工具条配色下拉: 默认 CPK; 切换后 GL 球色重绘。"""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    pv = PreviewView()
    assert pv._scheme.count() == len(SCHEMES)
    assert pv._scheme.currentIndex() == 0
    assert pv._scheme.itemData(0) == DEFAULT_SCHEME
    pv.set_scene(_scene())
    first_cpk = np.asarray(pv._gl._items[0].color)[0]  # GLScatterPlotItem 球色
    for i in range(1, pv._scheme.count()):
        pv._scheme.setCurrentIndex(i)
        other = np.asarray(pv._gl._items[0].color)[0]
        assert not np.allclose(first_cpk, other), \
            f"方案 {pv._scheme.itemData(i)} 切换后球色未变"
    # 场景在空态: 切换不崩(无场景可重绘)
    pv.clear()
    pv._scheme.setCurrentIndex(0)


# ---------- Task 2(2026-10-13 b 落点修正): b 环/相机中心绑靶 + 2D 取景对称化 ----------

def _pair_scene_roles(i_proj: int, i_targ: int):
    """手搭双片段场景(带角色): 片段 0 = 靶@近原点, 片段 1 = 弹丸@R_SEP 远处。"""
    from gui.core.preview import PreviewScene
    return PreviewScene(
        symbols=["H", "H", "H", "H"],
        masses=[1.0, 1.0, 1.0, 1.0],
        positions=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0],
                   [2.5, 2.5, 8.0], [3.2, 2.5, 8.0]],
        frag_slices=[(0, 2), (2, 4)], frag_names=["t", "p"],
        frag_coms=[[0.35, 0.0, 0.0], [2.85, 2.5, 8.0]],
        frag_vels=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        b_max=2.0, i_proj=i_proj, i_targ=i_targ,
    )


def test_gl_ring_center_binds_target():
    """GL b 环(72 点圆)质心 = 靶 COM(两序); 角色未解析 → 无环。"""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    gl = _Gl3D()
    for ip, it in ((1, 0), (0, 1)):
        gl.set_scene(_pair_scene_roles(ip, it))
        rings = [np.asarray(it_.pos) for it_ in gl._items
                 if hasattr(it_, "pos") and it_.pos is not None and len(it_.pos) == 72]
        assert len(rings) == 1, (ip, it)
        ctr_fit = (rings[0].min(axis=0) + rings[0].max(axis=0)) / 2  # 采样圆圆心
        assert ctr_fit == pytest.approx(_pair_scene_roles(ip, it).frag_coms[it], abs=0.01)
    gl.set_scene(_pair_scene_roles(-1, -1))  # 角色缺失: b_max 在场也不画环
    rings = [it_ for it_ in gl._items
             if hasattr(it_, "pos") and it_.pos is not None and len(it_.pos) == 72]
    assert rings == []


def test_gl_camera_center_binds_view_center():
    """相机中心(opts['center']) = view_center(靶 COM); 类型契约(QVector3D)不变。"""
    from PyQt6.QtWidgets import QApplication

    from gui.core.preview import view_center

    app = QApplication.instance() or QApplication([])
    gl = _Gl3D()
    for ip, it in ((1, 0), (0, 1)):
        sc = _pair_scene_roles(ip, it)
        gl.set_scene(sc)
        c = gl.view.opts["center"]
        assert all(callable(getattr(c, a, None)) for a in ("x", "y", "z"))
        assert [c.x(), c.y(), c.z()] == pytest.approx(view_center(sc), abs=1e-6)  # QVector3D 单精度
    sc = _pair_scene_roles(-1, -1)  # 回退: 全原子均值
    gl.set_scene(sc)
    c = gl.view.opts["center"]
    assert [c.x(), c.y(), c.z()] == pytest.approx(
        np.asarray(sc.positions, dtype=float).mean(axis=0), abs=1e-6)


def test_symmetric_window():
    """2D 取景对称化: 视图中心恒落窗几何中心, 原包围盒仍含于窗内。"""
    from gui.view3d.preview_view import symmetric_window

    lo = np.array([-1.0, -1.0])
    hi = np.array([3.0, 2.0])
    for c in (np.array([0.0, 0.0]), np.array([2.9, -0.5])):  # 含中心贴边的情形
        nlo, nhi = symmetric_window(lo, hi, c)
        assert np.allclose((nlo + nhi) / 2, c)
        assert np.all(nlo <= lo) and np.all(nhi >= hi)
        half = np.maximum(c - lo, hi - c)
        assert np.allclose(nhi - nlo, 2 * half)
