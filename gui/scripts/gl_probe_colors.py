# -*- coding: utf-8 -*-
"""GL 配色真机探针(2026-10-13 Task 5):真实(非 offscreen)窗口渲染多元素
场景,逐元素像素色读数落盘——供用户在真机目验 3D 原子配色按元素生效。

用法: python gui/scripts/gl_probe_colors.py
读数: gui/logs/gl_probe_colors.log(每方案×每元素: 方案色 / 散点件色 /
      命中像素数 / 判定),截图 gui/logs/gl_probe_colors_shot_<方案>.png
窗口约开 3 秒(每方案 1.2 秒渲染)后自动关闭;仅依赖 PyQt6 + pyqtgraph。
输出路径由脚本位置解析(仓库根/gui/logs),与 cwd 无关。
"""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

LOG = ROOT / "gui" / "logs" / "gl_probe_colors.log"
_log = LOG.open("w", encoding="utf-8")


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    _log.write(s + "\n")
    _log.flush()


def main() -> int:
    from PyQt6.QtCore import QT_VERSION_STR
    from PyQt6.QtWidgets import QApplication

    app = QApplication(sys.argv)
    import pyqtgraph
    p("platform:", app.platformName(), "| Qt", QT_VERSION_STR,
      "| pyqtgraph", pyqtgraph.__version__)
    if app.platformName() == "offscreen":
        p("FATAL: 当前为 offscreen 平台——请在真机(桌面会话)运行本探针")
        return 2

    from gui.core.preview import PreviewScene
    from gui.view3d.preview_view import SCHEMES, _Gl3D

    # 多元素场景:六种符号沿线排开(间距 3 Å,不成键);Zz = 表外符号走回退色
    syms = ["H", "O", "N", "F", "Au", "Zz"]
    scene = PreviewScene(
        symbols=syms, masses=[1.0] * 6,
        positions=[[3.0 * i, 0.0, 0.0] for i in range(6)],
        frag_slices=[(0, 6)], frag_names=["probe"],
        frag_coms=[[7.5, 0.0, 0.0]], frag_vels=[[0.0, 0.0, 0.0]],
    )
    gl = _Gl3D()
    gl.resize(720, 480)
    gl.show()

    def pump(sec: float):
        end = time.time() + sec
        while time.time() < end:
            app.processEvents()
            time.sleep(0.02)

    from pyqtgraph.opengl import GLScatterPlotItem

    def sym_of(it):
        """散点件 → 符号(按其原子位置集合匹配场景分组)。"""
        pts = {tuple(round(v, 6) for v in it.pos[k]) for k in range(len(it.pos))}
        for sym in syms:
            want = {tuple(round(v, 6) for v in scene.positions[k])
                    for k, s in enumerate(syms) if s == sym}
            if pts == want:
                return sym
        return None

    ok_all = True
    try:
        for key, sch in SCHEMES.items():
            gl.set_scheme(key)
            gl.set_scene(scene)
            pump(1.2)  # 真实窗口曝光若干帧(渲染在 paintGL)
            img = gl.view.grab().toImage()
            shot = ROOT / "gui" / "logs" / f"gl_probe_colors_shot_{key}.png"
            img.save(str(shot))
            hist = {sym: 0 for sym in syms}
            for y in range(0, img.height(), 2):
                for x in range(0, img.width(), 2):
                    px = img.pixel(x, y)
                    c = ((px >> 16) & 255, (px >> 8) & 255, px & 255)
                    for sym in syms:
                        w = sch["colors"].get(sym, sch["fallback"])
                        if all(abs(c[k] - w[k]) <= 30 for k in range(3)):
                            hist[sym] += 1
            p(f"--- 方案 {key}({sch['label']}) 背景 {tuple(round(v, 3) for v in gl.view.opts['bgcolor'])} ---")
            for sym in syms:
                w = sch["colors"].get(sym, sch["fallback"])
                it = next((i for i in gl._items
                           if isinstance(i, GLScatterPlotItem) and sym_of(i) == sym), None)
                icol = tuple(round(v, 3) for v in it.color) if it is not None else None
                ok = hist[sym] >= 4 and icol == tuple(
                    round(v / 255.0, 3) for v in w) + (1.0,)
                ok_all = ok_all and ok
                p(f"{sym:>3s} 方案色 {w} 件色 {icol} 命中像素 {hist[sym]:4d} -> {'OK' if ok else 'FAIL'}")
        gl.close()
        p("PROBE", "PASS" if ok_all else "FAIL",
          "| 全部元素球色按方案落屏" if ok_all else "| 存在未按方案落屏的元素(见上行 FAIL)")
        return 0 if ok_all else 1
    except Exception:
        import traceback
        p("FATAL:", traceback.format_exc())
        return 2


if __name__ == "__main__":
    sys.exit(main())
