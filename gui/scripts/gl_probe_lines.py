# -*- coding: utf-8 -*-
"""gl_probe_lines.py - 真机 GL 线元探针(2026-10-14 报障修复后的常备核验):
渲染带晶格/辅助线/轴的场景, 抓像素统计 gl_cell/gl_guide/gl_axis 三令牌色
落屏量(应全部显著>0), 并暴露 paint 期异常。
运行: python gui/scripts/gl_probe_lines.py  (需真窗口, offscreen 自拒)
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

os.environ.setdefault("QT_QPA_PLATFORM", "")  # 不许 offscreen
if os.environ.get("QT_QPA_PLATFORM") == "offscreen":
    print("FATAL: 需要真窗口平台")
    sys.exit(2)

import numpy as np
from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import QApplication

from gui import theme
from gui.core import preview, sysdir
from gui.core.inputwrite import parse_input_file
from gui.view3d.preview_view import _Gl3D, guide_elements, lattice_segments

LOG = ROOT / "gui" / "logs" / "gl_probe_lines.log"
lines = []


def _count(pix, hexcolor, tol=18):
    img = pix.toImage()
    hexc = hexcolor.lstrip("#")
    r0, g0, b0 = (int(hexc[k:k + 2], 16) for k in (0, 2, 4))
    n = 0
    for y in range(0, img.height(), 2):
        for x in range(0, img.width(), 2):
            c = img.pixelColor(x, y)
            if c is None:
                continue
            if (abs(c.red() - r0) < tol and abs(c.green() - g0) < tol
                    and abs(c.blue() - b0) < tol):
                n += 1
    return n


def scene():
    sd = sysdir.scan_system_dir(ROOT / "cases" / "case_spectrum_source"
                                / "h2ag_ebk" / "system")
    vals = dict(parse_input_file(ROOT / "cases" / "case_spectrum_source"
                                 / "h2ag_ebk" / "input_qct.txt"))
    vals.update(SCATTERING="SURFACE", N_THTA="0", THTA_MAX="0.6",
                N_CHI="0", CHI="0.9")
    return preview.assemble_preview(sd, vals, i_proj=0)


def run(gl: _Gl3D, tag: str):
    app.processEvents()
    pix = gl.grab()
    for name, key in (("gl_cell", "gl_cell"), ("gl_guide", "gl_guide"),
                      ("gl_axis", "gl_axis")):
        n = _count(pix, theme.COLORS[key])
        lines.append(f"{tag} {name}: {n}")
        print(lines[-1])


app = QApplication([])
theme.apply_theme(app)
sc = scene()

# A: 现行实现(set_scene 内部: mode='lines' 批段)
gl = _Gl3D()
gl.resize(760, 600)
gl.set_scene(sc)
gl.show()
QTimer.singleShot(1200, lambda: (run(gl, "A_current"), app.quit()))
app.exec()
gl.hide()



LOG.parent.mkdir(exist_ok=True)
with open(LOG, 'w', encoding='utf-8') as f:
    for l in lines:
        f.write(l + chr(10))
print("log ->", LOG)
