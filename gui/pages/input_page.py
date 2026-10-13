# -*- coding: utf-8 -*-
"""input_page.py - 页1:输入文件生成(左参数表单 / 右 3D 起始散射预览)

设计:
  表单键值语义 = "空 = 不写出"(与 Fortran 键表的可省默认对齐;必填键以加粗
  星号标于行标签,不加文字说明)。范式自适应:扫描 SYSTEM_DIR 后按表面散射/
  气相散射启用/禁用范式键(B_MAX/A_LAT/N_AIM/AIM_X/AIM_Y/角色键,条件必填
  的星号随范式动态增减),DIST_SCHEME 逐片段下拉重建,成员参数键按所选词
  动态显隐。校验=core.validate 全部静态规则;
  预览=core.preview 装配(core 公式移植,示意值注记显示在视图下方)。
  左侧分组可折叠(_Section):体系定义/采样方案/输出恒展开,运行控制可折叠
  默认展开,浴与谱/束流参数/成员参数/容器键默认折叠;折叠只藏可见性不禁用
  控件(已填值仍参与校验与写出);选中带成员键的方案词自动展开成员参数组;
  校验出错的键所在折叠组自动展开。
"""

from __future__ import annotations

import shutil
from pathlib import Path

from PyQt6.QtCore import QTimer
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QGridLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QPlainTextEdit, QScrollArea, QSplitter,
    QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout, QWidget,
    QHeaderView, QMessageBox,
)
from PyQt6.QtCore import Qt

from ..core import inputwrite, preview, sysdir, validate
from ..core.keyspec import (
    KEY_TABLE, MEMBER_TABLE, SCHEME_MEMBER_KEYS, ENUMS, REFERENCE_CONTAINER_KEYS,
)
from ..view3d.preview_view import PreviewView


class _Section(QWidget):
    """可折叠参数分组:标题条为 QToolButton(▸/▾ 指示),内容区整体显隐。

    折叠只藏可见性,不禁用控件——已填值仍参与校验与写出(与成员参数显隐
    同语义);不可折叠组标题无箭头且忽略切换。
    """

    def __init__(self, title: str, collapsible: bool = True,
                 expanded: bool = True, parent=None):
        super().__init__(parent)
        self.collapsible = collapsible
        self._btn = QToolButton()
        self._btn.setText(title)
        self._btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._btn.setStyleSheet("QToolButton { border: none; font-weight: 600; padding: 2px 0; }")
        if collapsible:
            self._btn.setCheckable(True)
            self._btn.setCursor(Qt.CursorShape.PointingHandCursor)
            self._btn.clicked.connect(lambda *_: self.toggle())
        self.content = QWidget()
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(2)
        v.addWidget(self._btn)
        v.addWidget(self.content)
        self._expanded = True
        self.set_expanded(expanded)

    @property
    def expanded(self) -> bool:
        return self._expanded

    def set_expanded(self, on: bool):
        self._expanded = on
        if not self.collapsible:
            return
        self._btn.setChecked(on)
        self._btn.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)
        self.content.setVisible(on)

    def toggle(self):
        if self.collapsible:
            self.set_expanded(not self._expanded)


class InputPage(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._w: dict[str, QLineEdit | QComboBox] = {}
        self._checks: dict[str, QCheckBox] = {}
        self._labels: dict[str, QLabel] = {}       # 行标签(星号标注用)
        self._label_base: dict[str, str] = {}      # 标签原文(星号增减的基础)
        self._scheme_combos: list[QComboBox] = []
        self._sysdef: sysdir.SysDef | None = None
        self._container_table: QTableWidget | None = None
        self._info: QLabel = QLabel("未扫描体系文件夹")
        self._notes: QLabel = QLabel("")
        self._notes.setWordWrap(True)
        self._verdict = QPlainTextEdit(readOnly=True)
        self._verdict.setMaximumHeight(120)
        self._view = PreviewView()
        self._sections: dict[str, _Section] = {}
        self._section_of_key: dict[str, _Section] = {}
        self._cur_sec: _Section | None = None
        self._build_ui()
        self._debounce = QTimer(self, singleShot=True, interval=350)
        self._debounce.timeout.connect(self.refresh_preview)

    # ---------- UI 构建 ----------
    def _build_ui(self):
        split = QSplitter(Qt.Orientation.Horizontal, self)
        left = QWidget()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        lay = QVBoxLayout(left)
        legend = QLabel(
            "注:<b>*</b> 加粗标记的字段不可留空;其余输入框留空即不写出该键。"
            "组标题带箭头者可点击折叠/展开(折叠不改变已填值的效力)。"
        )
        legend.setStyleSheet("color: #555;")
        lay.addWidget(legend)
        grid_row = 0

        def group(sid: str, title: str, collapsible: bool = False,
                  expanded: bool = True) -> tuple[_Section, QGridLayout]:
            sec = _Section(title, collapsible=collapsible, expanded=expanded)
            self._sections[sid] = sec
            self._cur_sec = sec
            lay.addWidget(sec)
            gl = QGridLayout(sec.content)
            gl.setVerticalSpacing(4)
            nonlocal grid_row
            grid_row = 0
            return sec, gl

        def add_row(gl: QGridLayout, label: str, w: QWidget, hint: str = "", key: str = ""):
            nonlocal grid_row
            lab = QLabel(label)
            if key:
                self._labels[key] = lab
                self._label_base[key] = label
            gl.addWidget(lab, grid_row, 0)
            gl.addWidget(w, grid_row, 1)
            if hint:
                gl.addWidget(QLabel(hint), grid_row, 2)
            grid_row += 1

        def line(key: str, ph: str = "") -> QLineEdit:
            e = QLineEdit()
            e.setPlaceholderText(ph)
            e.textChanged.connect(self._touch)
            self._w[key] = e
            self._section_of_key[key] = self._cur_sec
            return e

        def combo(key: str, items: list[str], editable: bool = False) -> QComboBox:
            c = QComboBox()
            c.addItems(items)
            c.setEditable(editable)
            if editable:
                c.setCurrentText("")
                c.lineEdit().setPlaceholderText("(默认)")
                c.lineEdit().textChanged.connect(self._touch)
            c.currentIndexChanged.connect(self._touch)
            c.currentTextChanged.connect(self._touch)
            self._w[key] = c
            self._section_of_key[key] = self._cur_sec
            return c

        def check(key: str) -> QCheckBox:
            c = QCheckBox()
            c.stateChanged.connect(self._touch)
            self._checks[key] = c
            self._section_of_key[key] = self._cur_sec
            return c

        # ---- 体系定义(先行:范式决定其余字段形态)----
        g, gl = group("sysdef", "体系定义(SYSTEM_DIR)")
        sd_row = QHBoxLayout()
        self._sd_line = QLineEdit()
        self._sd_line.setPlaceholderText("体系文件夹路径(含 .xyz / POSCAR 反应物文件)")
        b1 = QPushButton("浏览…")
        b1.clicked.connect(self._browse_sysdir)
        b2 = QPushButton("扫描")
        b2.clicked.connect(self._rescan)
        sd_row.addWidget(self._sd_line)
        sd_row.addWidget(b1)
        sd_row.addWidget(b2)
        wrap = QWidget()
        wrap.setLayout(sd_row)
        add_row(gl, "SYSTEM_DIR", wrap, key="SYSTEM_DIR")
        self._info.setStyleSheet("color: #555;")
        gl.addWidget(self._info, grid_row, 0, 1, 3)
        grid_row += 1
        add_row(gl, "E_REL [kcal/mol]", line("E_REL", "碰撞能"), key="E_REL")
        add_row(gl, "R_SEP [Å]", line("R_SEP", "初始间距"), key="R_SEP")
        add_row(gl, "B_MAX [Å]", line("B_MAX", "碰撞盘半径(仅气相散射)"))
        add_row(gl, "PROJECTILE", combo("PROJECTILE", [""]), key="PROJECTILE")
        add_row(gl, "TARGET", combo("TARGET", [""]), key="TARGET")
        add_row(gl, "A_LAT [Å]", line("A_LAT", "晶格常数(表面散射)"), key="A_LAT")
        add_row(gl, "N_AIM", combo("N_AIM", ["", "0", "1", "2"]), "0原点/1胞内/2定点(仅表面)")
        add_row(gl, "AIM_X [Å]", line("AIM_X", "N_AIM=2 定点 x(POSCAR 帧)"))
        add_row(gl, "AIM_Y [Å]", line("AIM_Y", "N_AIM=2 定点 y"))
        add_row(gl, "RANDOM_ORIENT", check("RANDOM_ORIENT"), "随机取向阶段(默认 F)")

        # ---- 采样方案(逐片段;扫描后重建;星号=每片段均须选择)----
        g, gl = group("scheme", "采样方案 DIST_SCHEME *(逐片段)")
        self._scheme_host = QWidget()
        self._scheme_grid = QGridLayout(self._scheme_host)
        gl.addWidget(self._scheme_host)

        # ---- 运行控制(DT 必填在内:可折叠但默认展开)----
        g, gl = group("run", "运行控制", collapsible=True)
        add_row(gl, "TITLE", line("TITLE", "标题(≤80 字符)"))
        add_row(gl, "N_TRAJ", line("N_TRAJ", "1(轨迹数)"))
        add_row(gl, "ISEED", line("ISEED", "1(种子)"))
        add_row(gl, "MAX_STEPS", line("MAX_STEPS", "100000"))
        add_row(gl, "DT [10 fs]", line("DT", "积分步长"), key="DT")
        add_row(gl, "INTEGRATOR", combo("INTEGRATOR", [""] + ENUMS["integ"]), "默认 VERLET")
        add_row(gl, "ORDER", line("ORDER", "1(radau 阶)"))
        add_row(gl, "DT_MIN / DT_MAX", line("DT_MIN", "0.0"), "自适应界(仅 radau;另一格 DT_MAX")
        self._dtmax = QLineEdit()
        self._dtmax.setPlaceholderText("0.0")
        self._dtmax.textChanged.connect(self._touch)
        self._section_of_key["DT_MAX"] = self._cur_sec
        gl.addWidget(self._dtmax, grid_row - 1, 2)
        add_row(gl, "ELEC_METHOD", line("ELEC_METHOD", "ADIABATIC(≤32)"))
        add_row(gl, "N_SURF", line("N_SURF", "0(折叠单绝热面)"))
        add_row(gl, "PES_MAIN", line("PES_MAIN", "1(主 PES 名)"))
        add_row(gl, "GWRITE_LEVEL", combo("GWRITE_LEVEL", ["", "0", "1", "2", "5"]), "默认 2")
        add_row(gl, "OBSERVABLES", line("OBSERVABLES", "记录列名逗号清单(≤256)"))
        eke = QHBoxLayout()
        for k in ("E_KIN", "E_POT", "E_TOT"):
            box = QCheckBox(k)
            box.stateChanged.connect(self._touch)
            self._checks[k] = box
            eke.addWidget(box)
        wrap2 = QWidget()
        wrap2.setLayout(eke)
        add_row(gl, "能量列", wrap2, "默认全 F")

        # ---- 浴与谱 ----
        g, gl = group("bath", "浴与谱", collapsible=True, expanded=False)
        add_row(gl, "BATH", combo("BATH", [""] + ENUMS["bath"]), "默认 NONE")
        add_row(gl, "N_BATH", line("N_BATH", "0(浴步数)"))
        add_row(gl, "DT_BATH [10 fs]", line("DT_BATH", "0.0(N_BATH>0 时须>0)"))
        add_row(gl, "T_BATH [K] / NU_COLL", line("T_BATH", "ANDERSEN 时须>0"))
        self._nucoll = QLineEdit()
        self._nucoll.setPlaceholderText("NU_COLL")
        self._nucoll.textChanged.connect(self._touch)
        self._w["NU_COLL"] = self._nucoll
        self._section_of_key["NU_COLL"] = self._cur_sec
        gl.addWidget(self._nucoll, grid_row - 1, 2)
        add_row(gl, "SPECTRUM_SOURCE", combo("SPECTRUM_SOURCE", [""] + ENUMS["spectrum"]), "默认 COMPUTE")
        add_row(gl, "SPECTRUM_FILE", line("SPECTRUM_FILE", "MANUAL 谱表路径"))

        # ---- 束流参数(入射通道)----
        g, gl = group("beam", "束流参数(入射通道)", collapsible=True, expanded=False)
        add_row(gl, "N_E_REL", combo("N_E_REL", ["", "0", "1", "2"]), "0热/1定值(默认)/2密度加权")
        add_row(gl, "T_TRANS [K]", line("T_TRANS", "N_E_REL=0 时须>0"))
        add_row(gl, "V_WIDTH [Å/(10fs)]", line("V_WIDTH", "N_E_REL=2 时须>0"))
        add_row(gl, "N_THTA", combo("N_THTA", ["", "0", "1"]), "0固定(默认)/1采样")
        add_row(gl, "THTA_MAX [rad]", line("THTA_MAX", "极角;固定与采样共用界"))
        add_row(gl, "N_CHI", combo("N_CHI", ["", "0", "1"]), "0固定(默认)/1均匀")
        add_row(gl, "CHI [rad]", line("CHI", "固定方位角(N_CHI=0)"))
        add_row(gl, "N_B", combo("N_B", ["", "0", "1"]), "0固定 b(默认)/1碰撞盘")

        # ---- 成员参数(动态显隐;束流/落点键已有专属控件,不重复建)----
        g, gl = group("member", "成员参数(随所选方案显示)", collapsible=True, expanded=False)
        self._member_host = QWidget()
        self._member_grid = QGridLayout(self._member_host)
        self._member_keys: set[str] = set()
        for mk in MEMBER_TABLE:
            if mk.key in self._w:  # N_E_REL/N_AIM/AIM_X 等在前述组已有专属控件
                continue
            e = QLineEdit()
            e.setPlaceholderText(mk.hint)
            e.textChanged.connect(self._touch)
            self._w[mk.key] = e
            self._member_keys.add(mk.key)
            self._section_of_key[mk.key] = self._cur_sec
            r = self._member_grid.rowCount()
            self._member_grid.addWidget(QLabel(mk.key), r, 0)
            self._member_grid.addWidget(e, r, 1)
        gl.addWidget(self._member_host)

        # ---- 容器键 ----
        g, gl = group("container", "容器键(活动容器的 *_keys 词表;reference: K_HARM/Q_TERM/E_SPLIT)",
                      collapsible=True, expanded=False)
        self._container_table = QTableWidget(0, 2)
        self._container_table.setHorizontalHeaderLabels(["键", "值(空=不写出)"])
        self._container_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        for k in REFERENCE_CONTAINER_KEYS:
            self._container_table.insertRow(self._container_table.rowCount())
            self._container_table.setItem(self._container_table.rowCount() - 1, 0, QTableWidgetItem(k))
            self._container_table.setItem(self._container_table.rowCount() - 1, 1, QTableWidgetItem(""))
        self._container_table.setMinimumHeight(110)
        gl.addWidget(self._container_table, 0, 0, 1, 2)

        # ---- 动作 ----
        g, gl = group("output", "输出")
        row = QHBoxLayout()
        b_val = QPushButton("校验")
        b_val.clicked.connect(self.run_validation)
        b_out = QPushButton("写出 input_qct.txt…")
        b_out.clicked.connect(self._write_file)
        self._copy_check = QCheckBox("同时复制 system/ 到输出目录")
        row.addWidget(b_val)
        row.addWidget(b_out)
        row.addWidget(self._copy_check)
        wrap3 = QWidget()
        wrap3.setLayout(row)
        gl.addWidget(wrap3, 0, 0, 1, 2)
        lay.addStretch(1)

        # 右侧:预览 + 注记 + 校验结果
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.addWidget(QLabel("起始散射构型(示意——内部几何为平衡坐标,抽样量取代表值)"))
        rl.addWidget(self._view, 1)
        rl.addWidget(self._notes)
        rl.addWidget(QLabel("校验结果"))
        rl.addWidget(self._verdict)

        split.addWidget(scroll)
        split.addWidget(right)
        split.setSizes([560, 640])
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.addWidget(split)
        # 无条件必填键的静态星号(DIST_SCHEME 的星号在组标题上)
        for k in ("SYSTEM_DIR", "E_REL", "R_SEP", "DT"):
            self._mark(k, True)

    def _mark(self, key: str, on: bool):
        """行标签的必填星号增减(加粗键名 + *;标签原文含单位后缀时星号置于键名后)。"""
        lab = self._labels.get(key)
        base = self._label_base.get(key, "")
        if lab is None or not base:
            return
        parts = base.split(" ", 1)
        if on:
            starred = f"<b>{parts[0]} *</b>" + (f" {parts[1]}" if len(parts) > 1 else "")
        else:
            starred = base
        lab.setText(starred)

    # ---------- 交互 ----------
    def _touch(self, *_):
        self._debounce.start()

    def _browse_sysdir(self):
        d = QFileDialog.getExistingDirectory(self, "选择 system/ 体系文件夹")
        if d:
            self._sd_line.setText(d)
            self._rescan()

    def _rescan(self):
        path = self._sd_line.text().strip()
        self._sysdef = sysdir.scan_system_dir(path) if path else None
        self._rebuild_dynamic()
        self._debounce.start()

    def _rebuild_dynamic(self):
        sd = self._sysdef
        # 清逐片段方案下拉
        while self._scheme_grid.count():
            it = self._scheme_grid.takeAt(0)
            w = it.widget()
            if w:
                w.deleteLater()
        self._scheme_combos = []
        if sd is not None and sd.ok:
            words = ENUMS["dist"]
            for i, f in enumerate(sd.frags):
                c = QComboBox()
                c.addItems([""] + words)
                c.currentIndexChanged.connect(self._touch)
                c.currentTextChanged.connect(self._touch)
                # 方案词变更即时刷新成员参数显隐(原先只在重扫描时刷新, 变更不生效)
                c.currentTextChanged.connect(self._refresh_member_visibility)
                self._scheme_combos.append(c)
                r = self._scheme_grid.rowCount()
                self._scheme_grid.addWidget(QLabel(f"片段 {i + 1}:{f.name}({f.nat} 原子)"), r, 0)
                self._scheme_grid.addWidget(c, r, 1)
            surf = sd.surface_model == 1
            n = len(sd.frags)
            # 角色下拉
            for key in ("PROJECTILE", "TARGET"):
                c = self._w[key]
                c.blockSignals(True)
                c.clear()
                c.addItems([""] + sd.frag_names if not (surf and key == "TARGET") else [""])
                c.blockSignals(False)
            # 范式启停
            self._w["B_MAX"].setEnabled(not surf)
            self._w["B_MAX"].setPlaceholderText("碰撞盘半径(仅气相散射)" if not surf else "表面散射不支持")
            for k in ("A_LAT", "N_AIM", "AIM_X", "AIM_Y"):
                self._w[k].setEnabled(surf or k != "N_AIM")
            self._w["N_AIM"].setEnabled(surf)
            self._w["TARGET"].setEnabled(not surf and n > 1)
            self._w["PROJECTILE"].setEnabled(surf or n > 1)
            # 条件必填星号随范式增减(表面:PROJECTILE+A_LAT;气相多片段:双角色)
            self._mark("PROJECTILE", surf or n > 1)
            self._mark("TARGET", (not surf) and n > 1)
            self._mark("A_LAT", surf)
            para = "表面散射(带胞)" if surf else f"气相散射({n} 片段)"
            box = f"; box {sd.box_lx:.4f}×{sd.box_ly:.4f} Å, skew {sd.skew * 57.2958:.1f}°" if surf else ""
            self._info.setText(f"✓ {para}{box} — 片段: {', '.join(sd.frag_names)}")
            # 范式切换清理:本范式禁写的键残留值自动清空(禁写字段的旧值不该变成雷)
            illegal_when_surf = ["B_MAX"]
            illegal_when_mol = ["A_LAT", "AIM_X", "AIM_Y"]
            for k in (illegal_when_surf if surf else illegal_when_mol):
                w = self._w[k]
                if isinstance(w, QLineEdit) and w.text().strip():
                    w.clear()
                elif isinstance(w, QComboBox) and w.currentText().strip():
                    w.setCurrentIndex(0)
            if surf:
                self._w["N_AIM"].setCurrentIndex(0)
            if n == 1:  # 单片段气相:两角色键都不合法
                for k in ("PROJECTILE", "TARGET"):
                    self._w[k].setCurrentIndex(0)
        else:
            errs = sd.errors[0] if sd is not None and sd.errors else "未扫描"
            self._info.setText(f"⚠ {errs}")
        self._refresh_member_visibility()

    def _refresh_member_visibility(self):
        selected: set[str] = set()
        for c in self._scheme_combos:
            t = c.currentText().strip()
            if t:
                selected.add(t)
        visible = set()
        for w in selected:
            visible.update(SCHEME_MEMBER_KEYS.get(w, []))
        for key in self._member_keys:
            self._w[key].setVisible(key in visible)
        for i in range(self._member_grid.count()):
            item = self._member_grid.itemAt(i)
            w = item.widget()
            if isinstance(w, QLabel):
                # 标签文本即键名(可能带单位后缀)
                key = w.text().split()[0]
                w.setVisible(key in visible)
        # 选中带成员键的方案词后自动展开(否则新显行藏在折叠组里不可见)
        if visible and not self._sections["member"].expanded:
            self._sections["member"].set_expanded(True)

    # ---------- 数据汇聚 ----------
    def collect(self) -> tuple[dict[str, str], dict[str, str]]:
        values: dict[str, str] = {}
        for k, w in self._w.items():
            if isinstance(w, QLineEdit):
                t = w.text().strip()
                if t:
                    values[k] = t
            else:  # QComboBox
                t = w.currentText().strip()
                if t and k != "PROJECTILE" and k != "TARGET":
                    values[k] = t
        # 角色键按原文(区分大小写按文件名匹配,大小写不敏感)
        for k in ("PROJECTILE", "TARGET"):
            t = self._w[k].currentText().strip()
            if t:
                values[k] = t
        # 复选框
        for k, c in self._checks.items():
            values[k] = "T" if c.isChecked() else "F"
        # 复合行:DT_MIN 与 DT_MAX、THTA 等
        if self._dtmax.text().strip():
            values["DT_MAX"] = self._dtmax.text().strip()
        # DIST_SCHEME
        words = [c.currentText().strip() for c in self._scheme_combos]
        if all(words) and words:
            values["DIST_SCHEME"] = ",".join(words)
        values.setdefault("SYSTEM_DIR", self._sd_line.text().strip())
        # 容器键表
        cont: dict[str, str] = {}
        t = self._container_table
        for r in range(t.rowCount()):
            k = (t.item(r, 0).text().strip() if t.item(r, 0) else "").upper()
            v = t.item(r, 1).text().strip() if t.item(r, 1) else ""
            if k and v:
                cont[k] = v
        return values, cont

    def _roles(self) -> tuple[int, int]:
        sd = self._sysdef
        i_proj = i_targ = -1
        if sd is not None and sd.ok:
            p = self._w["PROJECTILE"].currentText().strip().lower()
            t = self._w["TARGET"].currentText().strip().lower()
            names = [f.name.lower() for f in sd.frags]
            if p in names:
                i_proj = names.index(p)
            if t in names:
                i_targ = names.index(t)
        return i_proj, i_targ

    # ---------- 动作 ----------
    def refresh_preview(self):
        sd = self._sysdef
        if sd is None or not sd.ok:
            self._view.clear()
            self._notes.setText("")
            return
        values, _ = self.collect()
        i_proj, i_targ = self._roles()
        try:
            scene = preview.assemble_preview(sd, values, i_proj=i_proj, i_targ=i_targ)
        except Exception as exc:  # 预览失败不阻塞表单
            self._notes.setText(f"预览暂不可用:{exc}")
            return
        self._view.set_scene(scene)
        self._notes.setText(" ; ".join(scene.notes))

    def run_validation(self):
        values, cont = self.collect()
        issues = validate.validate(values, self._sysdef, list(cont.keys()))
        self._expand_for_issues(issues, cont)
        if not issues:
            self._verdict.setPlainText("✓ 校验通过(全部静态规则)")
        else:
            self._verdict.setPlainText("\n".join(f"✗ {i}" for i in issues))

    def _expand_for_issues(self, issues, cont_keys: dict[str, str]):
        """校验出错的键所在折叠组自动展开(错误不该藏在折叠里);容器类问题归容器键组。"""
        cont_names = {k.upper() for k in cont_keys}
        for i in issues:
            sec = self._section_of_key.get(i.key)
            if sec is None and (not i.key or i.key in cont_names):
                sec = self._sections.get("container")
            if sec is not None and sec.collapsible and not sec.expanded:
                sec.set_expanded(True)

    def _write_file(self):
        values, cont = self.collect()
        issues = validate.validate(values, self._sysdef, list(cont.keys()))
        if issues:
            self.run_validation()
            QMessageBox.warning(self, "校验未通过", "存在校验错误,先修正再写出(见校验结果)。")
            return
        out_dir = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if not out_dir:
            return
        sysdir_val = "system" if self._copy_check.isChecked() else self._sd_line.text().strip()
        values = dict(values)
        values["SYSTEM_DIR"] = sysdir_val
        header = ["# 由 VENUS2026 GUI 生成(input_qct.txt)"]
        title = values.get("TITLE", "").strip()
        if title:
            header.append(f"# {title}")
        try:
            p = inputwrite.write_input_file(
                out_dir + "/input_qct.txt", values, cont, header_lines=header
            )
            if self._copy_check.isChecked() and self._sysdef is not None:
                src = self._sysdef.folder
                dst = str(Path(out_dir) / "system")
                if not Path(dst).exists():
                    shutil.copytree(src, dst)
            QMessageBox.information(self, "已写出", f"输入文件已写出:\n{p}")
        except Exception as exc:
            QMessageBox.critical(self, "写出失败", str(exc))
