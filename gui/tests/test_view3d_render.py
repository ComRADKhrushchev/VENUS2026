# -*- coding: utf-8 -*-
"""Task 1(2026-10-13 3D 预览渲染增强)单测:

  - nice_step: 1/2/5×10^n 整洁刻度步长(纯函数, 不需 QApplication);
  - SCHEMES: 两方案同键域(现 CPK 键集 + 未知元素回退), CPK 浅底适配偏离已锁,
    gl_bg/gl_axis 令牌浅而不近白/浅底可辨;
  - _Gl3D: offscreen 构造 + set_scene 不崩, opts 契约(center Vector /
    distance 数值)不变, 背景落 gl_bg 令牌, 轴件在场, 刻度角标文本
    = nice_step 派生值;
  - PreviewView: 配色下拉默认 CPK, 切换后 GL 球色重绘。
  - Task 5(2026-10-13 真机配色修复): 原子为逐符号单色 GLScatterPlotItem
    (均匀色元组+标量尺寸+opaque 混合), 每符号一件, 色=方案查表值
    (未知符号走回退), 两方案都断言。
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
    RING_SEGMENTS,
    SCHEMES,
    PreviewView,
    _Gl3D,
    _display_radius,
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
    """PreviewView 工具条配色下拉: 默认 CPK; 切换后逐符号球色随方案重绘。"""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    pv = PreviewView()
    assert pv._scheme.count() == len(SCHEMES)
    assert pv._scheme.currentIndex() == 0
    assert pv._scheme.itemData(0) == DEFAULT_SCHEME
    sc = _scene()  # 符号 H/O/H → 逐符号两件
    pv.set_scene(sc)

    def sym_color(sym: str):
        for it in _scatter_items(pv._gl):
            if _item_symbol(it, sc) == sym:
                return np.asarray(it.color, dtype=float)[:3]
        raise AssertionError(f"未找到符号 {sym} 的散点件")

    for i in range(pv._scheme.count()):
        pv._scheme.setCurrentIndex(i)
        key = pv._scheme.itemData(i)
        for sym in ("H", "O"):
            assert sym_color(sym) == pytest.approx(
                [c / 255.0 for c in SCHEMES[key]["colors"][sym]]), (key, sym)
    # 场景在空态: 切换不崩(无场景可重绘)
    pv.clear()
    pv._scheme.setCurrentIndex(0)


# ---------- Task 5(2026-10-13 真机配色修复): 逐符号单色件 ----------

def _scatter_items(gl):
    """set_scene 后的原子散点件(键/箭头/环是 GLLinePlotItem, 排除)。"""
    from pyqtgraph.opengl import GLScatterPlotItem

    return [it for it in gl._items if isinstance(it, GLScatterPlotItem)]


def _item_symbol(it, scene):
    """把逐符号散点件映回符号(按其原子位置集合匹配场景符号分组)。"""
    pts = {tuple(np.round(p, 6)) for p in np.asarray(it.pos, dtype=float)}
    groups: dict[str, set] = {}
    for k, s in enumerate(scene.symbols):
        groups.setdefault(s, set()).add(tuple(np.round(scene.positions[k], 6)))
    for sym, want in groups.items():
        if pts == want:
            return sym
    return None


def test_gl_atoms_per_symbol_items_colors():
    """原子逐符号单色件: 每符号一件; 色=方案查表值(0-1 均匀元组), 未知符号
    走回退; 标量尺寸+pxMode=False; 两方案都断言; clear 清场。"""
    from PyQt6.QtWidgets import QApplication

    from gui.core.preview import PreviewScene

    app = QApplication.instance() or QApplication([])
    sc = PreviewScene(
        symbols=["H", "O", "H", "Zz"],  # Zz = 未知符号 → 回退色
        masses=[1.0, 16.0, 1.0, 100.0],
        positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 1.0], [0.0, 3.5, 4.2], [4.0, 4.0, 0.0]],
        frag_slices=[(0, 4)], frag_names=["mix"],
    )
    gl = _Gl3D()
    for key, sch in SCHEMES.items():
        gl.set_scheme(key)
        gl.set_scene(sc)
        scats = _scatter_items(gl)
        assert len(scats) == 3, (key, len(scats))  # H / O / Zz 各一件
        seen = set()
        for it in scats:
            sym = _item_symbol(it, sc)
            assert sym is not None, (key, np.asarray(it.pos).tolist())
            seen.add(sym)
            want = sch["colors"].get(sym, sch["fallback"])
            assert np.asarray(it.color) == pytest.approx(
                [c / 255.0 for c in want] + [1.0]), (key, sym)
            assert it.size == pytest.approx(_display_radius(sym)), (key, sym)
            assert it.pxMode is False, (key, sym)
        assert seen == {"H", "O", "Zz"}, (key, seen)
    gl.clear()
    assert _scatter_items(gl) == []


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
    """GL b 环(RING_SEGMENTS 点圆)质心 = 靶 COM(两序); 角色未解析 → 无环。"""
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    gl = _Gl3D()
    for ip, it in ((1, 0), (0, 1)):
        gl.set_scene(_pair_scene_roles(ip, it))
        rings = [np.asarray(it_.pos) for it_ in gl._items
                 if hasattr(it_, "pos") and it_.pos is not None and len(it_.pos) == RING_SEGMENTS]
        assert len(rings) == 1, (ip, it)
        ctr_fit = (rings[0].min(axis=0) + rings[0].max(axis=0)) / 2  # 采样圆圆心
        assert ctr_fit == pytest.approx(_pair_scene_roles(ip, it).frag_coms[it], abs=0.01)
    gl.set_scene(_pair_scene_roles(-1, -1))  # 角色缺失: b_max 在场也不画环
    rings = [it_ for it_ in gl._items
             if hasattr(it_, "pos") and it_.pos is not None and len(it_.pos) == RING_SEGMENTS]
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


# ---------- 周期晶格网格(2026-10-14: 覆盖表面片段足迹, 取代原点单胞框) ----------

def test_lattice_segments_covers_footprint():
    from gui.view3d.preview_view import lattice_segments
    cell = [[6.0, 0, 0], [0, 6.0, 0], [0, 0, 20.0]]
    pos = np.array([[0.5, 0.5, 0.1], [13.2, 7.4, 2.0]])
    segs, zlo, zhi = lattice_segments(cell, pos)
    assert segs.ndim == 3 and segs.shape[1] == 2 and segs.shape[2] == 3
    assert (zlo, zhi) == (pytest.approx(0.1), pytest.approx(2.0))
    pts = segs.reshape(-1, 3)
    assert pts[:, 0].min() == pytest.approx(0.0)   # floor(0.5/6)=0
    assert pts[:, 0].max() == pytest.approx(18.0)  # ceil(13.2/6)=3 胞
    assert pts[:, 1].max() == pytest.approx(12.0)  # ceil(7.4/6)=2 胞


def test_lattice_segments_cap_and_degenerate():
    from gui.view3d.preview_view import lattice_segments
    cell = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    pos = np.array([[0.5, 0.5, 0.0], [500.5, 500.5, 0.0]])
    segs, _, _ = lattice_segments(cell, pos, cap=4)
    pts = segs.reshape(-1, 3)
    assert pts[:, 0].max() - pts[:, 0].min() == pytest.approx(4.0), "每方向跨度封顶 cap"
    # 奇异横向晶矢: 空段
    segs0, z0, z1 = lattice_segments([[1, 0, 0], [2, 0, 0], [0, 0, 1]], np.zeros((1, 3)))
    assert len(segs0) == 0


def test_lattice_segments_skewed_corners():
    """斜胞: 格点是 a1/a2 平行四边形组合(corner(1,1)=a1+a2)。"""
    from gui.view3d.preview_view import lattice_segments
    a1, a2 = np.array([2.0, 0.0]), np.array([1.0, np.sqrt(3.0)])
    cell = [[a1[0], a1[1], 0], [a2[0], a2[1], 0], [0, 0, 5.0]]
    pos = np.array([[1.0, 0.6, 0.5]])
    segs, _, _ = lattice_segments(cell, pos)
    pts = segs.reshape(-1, 3)[:, :2]
    assert np.any(np.all(np.isclose(pts, a1 + a2), axis=1)), "corner(1,1)=a1+a2 应在场"


# ---------- 入射角辅助线(2026-10-14: θ/χ 弧 + 法线/方位基准, 表面抬层防遮挡) ----------


def _pair_scene(thta="0.35"):
    from conftest import ROOT
    from gui.core import preview, sysdir
    from gui.core.inputwrite import parse_input_file
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")
    vals = dict(parse_input_file(ROOT / "cases" / "case_harmonic_scatter" / "input_qct.txt"))
    vals.update(THTA_MAX=thta, N_THTA="0")
    return preview.assemble_preview(sd, vals, i_proj=0, i_targ=1)


def test_slerp_arc_endpoints_and_midpoint():
    from gui.view3d.preview_view import _slerp_arc
    arc = _slerp_arc((0, 0, 0), (0, 0, 1), (1, 0, 0), 2.0, n=25)
    assert len(arc) == 25
    assert np.allclose(arc[0], [0, 0, 2])
    assert np.allclose(arc[-1], [2, 0, 0])
    mid = 2 * np.array([np.cos(np.pi / 4), 0.0, np.sin(np.pi / 4)])
    assert np.allclose(arc[12], mid, atol=1e-9), "n=25 时索引 12 恰为半程"


def test_guide_elements_pair_anchor_zero():
    from gui.view3d.preview_view import guide_elements
    gd = guide_elements(_pair_scene())
    assert gd is not None
    assert gd["anchor"][2] == pytest.approx(0.0), "气相(pair)锚在 z=0 层"
    assert gd["theta"] == pytest.approx(0.35)
    assert len(gd["theta_arc"]) >= 2 and len(gd["normal_seg"]) == 2


def test_guide_elements_surface_lift():
    """表面范式: 有衬底原子→锚抬至衬底顶面+0.4; 单片段表面(无衬底)→ z=0。"""
    import tempfile
    from pathlib import Path as P
    from gui.core import preview, sysdir
    from gui.view3d.preview_view import guide_elements
    with tempfile.TemporaryDirectory() as td:
        sysd = P(td) / "system"
        sysd.mkdir()
        (sysd / "slab.poscar").write_text('slab\n1.0\n6 0 0\n0 6 0\n0 0 8\nH\n1\nCartesian\n3 3 1.5\n', encoding="utf-8")
        (sysd / "proj.xyz").write_text('1\n\nH 0 0 0\n', encoding="utf-8")
        sd = sysdir.scan_system_dir(sysd)
        vals = {"SCATTERING": "SURFACE", "CELL": "6 0 0, 0 6 0, 0 0 8",
                "E_REL": "10", "R_SEP": "6", "THTA_MAX": "0.5", "N_THTA": "0",
                "DIST_SCHEME": "BOLTZMANN,BOLTZMANN"}
        sc = preview.assemble_preview(sd, vals, i_proj=0)   # 名序: proj.xyz 为片段 0
        gd = guide_elements(sc)
        assert gd is not None
        assert gd["anchor"][2] == pytest.approx(1.5 + 0.4), "衬底顶面 1.5 + 间隙 0.4"
        # 单片段表面(无衬底原子): 锚回 z=0
        solo = P(td) / "solo"
        solo.mkdir()
        xyz = chr(10).join(["1", "", "H 0 0 0", ""])
        (solo / "proj.xyz").write_text(xyz, encoding="utf-8")
        sd2 = sysdir.scan_system_dir(solo)
        sc2 = preview.assemble_preview(sd2, dict(vals, DIST_SCHEME="BOLTZMANN"), i_proj=0)
        gd2 = guide_elements(sc2)
        assert gd2 is not None and gd2["anchor"][2] == pytest.approx(0.0)


def test_scale_label_carries_angles():
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    gl = _Gl3D()
    gl.set_scene(_pair_scene())
    txt = gl._scale_label.text()
    assert "θ=" in txt and "χ=" in txt and "°" in txt, txt
    assert "20.1" in txt, "θ=0.35 rad ≈ 20.1°"


def test_gl_line_items_have_valid_color_tuples():
    """真 GL 看守(2026-10-14 实机报障): GLLinePlotItem 的元组色必须 4 元
    (r,g,b,a)——_rgb01 已含 alpha, 再追加会让 glVertexAttrib4f 抛 TypeError
    (paint 期逐项被吞→元素不可见; offscreen 不渲染 paint 故以此看守)。"""
    from PyQt6.QtWidgets import QApplication
    from pyqtgraph.opengl import GLLinePlotItem
    app = QApplication.instance() or QApplication([])
    gl = _Gl3D()
    gl.set_scene(_pair_scene())
    n = 0
    for it in gl._items:
        if isinstance(it, GLLinePlotItem) and isinstance(it.color, tuple):
            assert len(it.color) == 4, f"线元色组 {len(it.color)} 元: {it.color}"
            n += 1
    assert n >= 3, "场景应含多条元组色线元(键/晶格/辅助线)"
