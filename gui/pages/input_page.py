# -*- coding: utf-8 -*-
"""input_page.py - 页1:输入文件生成(左参数表单 / 右 3D 起始散射预览)

设计:
  表单键值语义 = "空 = 不写出"(与 Fortran 键表的可省默认对齐;必填键以加粗
  星号标于行标签,不加文字说明)。行标签显示物理量(PHYSICAL_LABELS),KEY 名
  入 tooltip。条件显隐:被其余设置判"不需要"的行整行隐藏并清空(_apply_cond_
  visibility;隐藏键不进 collect,绝不写出),规则权威 = core.validate 移植的
  Fortran 规则 + propagator.f90 的积分器×阶域;ORDER 随积分器变下拉。成员参
  数按所选词动态显隐(词基显隐不清空、旧值可收集——与条件显隐分开记账)。
  校验=core.validate 全部静态规则;
  预览=core.preview 装配(core 公式移植,示意值注记显示在视图下方)。
  左侧分组可折叠(_Section):体系定义/采样方案/输出恒展开,运行控制可折叠
  默认展开,浴与谱/束流参数/成员参数/容器键默认折叠;折叠只藏可见性不禁用
  控件(已填值仍参与校验与写出);选中带成员键的方案词自动展开成员参数组;
  校验出错的键所在折叠组自动展开。
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

from PyQt6.QtCore import QSize, QTimer
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QFrame, QGridLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QPlainTextEdit, QScrollArea, QSizePolicy,
    QSplitter, QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout,
    QWidget, QHeaderView, QMessageBox,
)
from PyQt6.QtCore import Qt

from .. import theme
from ..core import inputwrite, preview, sysdir, validate
from ..core.keyspec import (
    KEY_TABLE, MEMBER_TABLE, SCHEME_MEMBER_KEYS, ENUMS, REFERENCE_CONTAINER_KEYS,
)
from ..view3d.preview_view import PreviewView

# 行标签显示物理量(2026-10-13 用户裁决:不在输入框左侧显示 KEY 名);KEY 名入 tooltip
PHYSICAL_LABELS = {
    "SYSTEM_DIR": "体系文件夹",
    "E_REL": "碰撞能 [eV]",
    "R_SEP": "初始间距 [Å]",
    "B_MAX": "碰撞参数 b [Å]",
    "PROJECTILE": "弹丸片段",
    "TARGET": "靶片段",
    "A_LAT": "晶格常数 [Å]",
    "N_AIM": "落点模式",
    "AIM_X": "定点落点 x [Å]",
    "AIM_Y": "定点落点 y [Å]",
    "RANDOM_ORIENT": "随机取向",
    "TITLE": "标题",
    "N_TRAJ": "轨迹数",
    "ISEED": "随机种子",
    "MAX_STEPS": "最大步数",
    "DT": "积分步长 [fs]",
    "INTEGRATOR": "积分器",
    "ORDER": "积分阶",
    "DT_MIN": "自适应步长下界 [fs]",
    "DT_MAX": "步长上界 [fs]",
    "ELEC_METHOD": "电子方法名",
    "GWRITE_LEVEL": "输出详略级别",
    "OBSERVABLES": "记录列清单",
    "BATH": "热浴",
    "N_BATH": "浴步数",
    "DT_BATH": "浴步长 [fs]",
    "T_BATH": "浴温度 [K]",
    "NU_COLL": "碰撞频率",
    "SPECTRUM_SOURCE": "谱来源",
    "SPECTRUM_FILE": "谱表文件路径",
    "N_E_REL": "碰撞能模式",
    "T_TRANS": "平动温度 [K]",
    "V_WIDTH": "速度展宽 [Å/fs]",
    "N_THTA": "极角模式",
    "THTA_MAX": "极角 [rad]",
    "N_CHI": "方位角模式",
    "CHI": "固定方位角 [rad]",
    "N_B": "瞄准模式",
}

# 积分器 → 合法阶域(权威 src/dynamics/propagator.f90:126-165/186/223-224)
ORDER_ITEMS = {
    "VERLET": ["1", "2"],
    "SYMPLECTIC": ["4", "6", "8"],
    "RADAU": [str(i) for i in range(1, 16)],
}

# 模式类下拉: (值, 含义) — 选项显示"值 · 含义", 写出取原值(itemData)
MODE_OPTS = {
    "N_E_REL": [("0", "热分布"), ("1", "定值(默认)"), ("2", "密度加权")],
    "N_THTA": [("0", "固定(默认)"), ("1", "sinθ 采样")],
    "N_CHI": [("0", "固定(默认)"), ("1", "均匀采样")],
    "N_B": [("0", "固定 b(默认)"), ("1", "碰撞盘")],
    "N_AIM": [("0", "原点"), ("1", "胞内均匀"), ("2", "定点")],
    "N_OSC": [("0", "热"), ("1", "定态")],
    "N_E_BAR": [("0", "定值"), ("1", "热")],
}


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
        self._btn.setObjectName("SectionHeader")  # 样式由全局 QSS 收口
        self._btn.setText(title)
        self._btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._btn.setIconSize(QSize(14, 14))
        if collapsible:
            self._btn.setCheckable(True)
            self._btn.setCursor(Qt.CursorShape.PointingHandCursor)
            self._btn.clicked.connect(lambda *_: self.toggle())
        self.content = QWidget()
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(theme.SP_S)
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
        self._btn.setIcon(theme.icon("chevron-down" if on else "chevron-right"))
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
        self._info.setWordWrap(True)  # 扫描结果长行不撑最小宽
        self._info.setProperty("class", "caption")
        self._notes = QLabel("")
        self._notes.setWordWrap(True)
        self._notes.setProperty("class", "caption")
        self._verdict = QPlainTextEdit(readOnly=True)
        self._verdict.setMaximumHeight(120)
        self._verdict.setProperty("state", "idle")  # ok/error 时语义色(QSS)
        self._verdict_summary = QLabel("")  # 校验语义摘要行(● 通过 / ● N 处问题)
        self._verdict_summary.setProperty("class", "verdict-summary")
        self._verdict_summary.setProperty("state", "idle")
        self._view = PreviewView()
        self._sections: dict[str, _Section] = {}
        self._section_of_key: dict[str, _Section] = {}
        self._cell_of_key: dict[str, QWidget] = {}        # 键 → 单元对(标签+控件)
        self._rows: list[tuple[QWidget, list[str]]] = []  # 行容器 → 行内键(整行塌缩判定)
        self._pending: list[str] = []                     # row() 取走的键缓冲
        self._cond_hidden: set[str] = set()               # 条件隐藏(collect 跳过,绝不写出)
        self._word_hidden: set[str] = set()               # 成员词未选而隐藏(旧语义,可收集)
        self._cur_sec: _Section | None = None
        self._cur_lay: QVBoxLayout | None = None
        self._build_ui()
        self._debounce = QTimer(self, singleShot=True, interval=350)
        self._debounce.timeout.connect(self.refresh_preview)

    # ---------- UI 构建 ----------
    def _build_ui(self):
        split = QSplitter(Qt.Orientation.Horizontal, self)
        left = QWidget()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        # 横向滚动条按需出现(2026-10-13 窗口适配裁决: 用户拉窄分栏时作后备,
        # 默认分栏份额下内容最小宽 502-568px < 可视宽, 平时不出现)
        scroll.setWidget(left)
        self._scroll = scroll
        lay = QVBoxLayout(left)
        lay.setContentsMargins(theme.SP_M, theme.SP_M, theme.SP_M, theme.SP_M)
        lay.setSpacing(theme.SP_M)
        legend = QLabel(
            "注:<b>*</b> 加粗标记的字段不可留空;其余输入框留空即不写出该键。"
            "组标题带箭头者可点击折叠/展开(折叠不改变已填值的效力);"
            "不需要的输入框按其余设置自动隐藏并清空(各键名见悬停提示)。"
        )
        legend.setWordWrap(True)  # 长图例不撑最小宽
        legend.setProperty("class", "caption")
        lay.addWidget(legend)

        def group(sid: str, title: str, collapsible: bool = False,
                  expanded: bool = True) -> _Section:
            sec = _Section(title, collapsible=collapsible, expanded=expanded)
            self._sections[sid] = sec
            self._cur_sec = sec
            lay.addWidget(sec)
            self._cur_lay = QVBoxLayout(sec.content)
            self._cur_lay.setSpacing(theme.SP_S)
            self._cur_lay.setContentsMargins(
                theme.SP_M, theme.SP_S, theme.SP_M, theme.SP_S)  # 组头下内容缩进
            return sec

        def cell(key: str, w: QWidget, hint: str = "", label: str | None = None) -> QWidget:
            """单元对 = 行标签(物理量)+控件;旁注并入悬停;登记键→单元/段映射。"""
            text = label if label is not None else PHYSICAL_LABELS.get(key, key)
            lab = QLabel(text)
            tip = f"写出键: {key}" + (f"\n{hint}" if hint else "")
            lab.setToolTip(tip)
            w.setToolTip(tip)
            self._labels[key] = lab
            self._label_base[key] = text
            c = QWidget()
            h = QHBoxLayout(c)
            h.setContentsMargins(0, 0, 0, 0)
            h.setSpacing(4)
            h.addWidget(lab)
            h.addWidget(w, 1)
            self._cell_of_key[key] = c
            self._section_of_key[key] = self._cur_sec
            self._pending.append(key)
            return c

        def row(*cells: QWidget):
            """行容器: 1-3 个单元对并排;隐藏某键的单元对时其余左移,整行键全
            隐藏时行容器塌缩(_apply_row_visibility)。"""
            r = QWidget()
            h = QHBoxLayout(r)
            h.setContentsMargins(0, 0, 0, 0)
            h.setSpacing(12)
            for c in cells:
                h.addWidget(c, 1)
            self._cur_lay.addWidget(r)
            self._rows.append((r, list(self._pending)))
            self._pending.clear()

        def line(key: str, ph: str = "") -> QLineEdit:
            e = QLineEdit()
            e.setPlaceholderText(ph)
            e.textChanged.connect(self._touch)
            self._w[key] = e
            return e

        def combo(key: str, items: list[str], editable: bool = False) -> QComboBox:
            c = QComboBox()
            for it in items:               # 普通枚举: 显示 = 值(存 itemData)
                c.addItem(it, it)
            c.setEditable(editable)
            if editable:
                c.setCurrentText("")
                c.lineEdit().setPlaceholderText("(默认)")
                c.lineEdit().textChanged.connect(self._touch)
            c.currentIndexChanged.connect(self._touch)
            c.currentTextChanged.connect(self._touch)
            self._w[key] = c
            return c

        def mcombo(key: str) -> QComboBox:
            """模式下拉: 首项空白=默认不写出;选项显示"值 · 含义", 写出取原值(itemData)。"""
            c = combo(key, [""])
            for v, desc in MODE_OPTS[key]:
                c.addItem(f"{v} · {desc}", v)
            return c

        def check(key: str) -> QCheckBox:
            c = QCheckBox()
            c.stateChanged.connect(self._touch)
            self._checks[key] = c
            return c

        # ---- 体系定义(先行:范式决定其余字段形态)----
        group("sysdef", "体系定义(SYSTEM_DIR)")
        sd_row = QHBoxLayout()
        self._sd_line = QLineEdit()
        self._sd_line.setPlaceholderText("体系文件夹路径(含 .xyz / POSCAR 反应物文件)")
        b1 = QPushButton("浏览…")
        b1.setIcon(theme.icon("folder"))
        b1.clicked.connect(self._browse_sysdir)
        b2 = QPushButton("扫描")
        b2.setIcon(theme.icon("search"))
        b2.clicked.connect(self._rescan)
        sd_row.addWidget(self._sd_line)
        sd_row.addWidget(b1)
        sd_row.addWidget(b2)
        wrap = QWidget()
        wrap.setLayout(sd_row)
        row(cell("SYSTEM_DIR", wrap))
        self._cur_lay.addWidget(self._info)
        row(cell("E_REL", line("E_REL", "如 8.7"), "必填"),
            cell("R_SEP", line("R_SEP", "如 8.0"), "必填"))
        row(cell("B_MAX", line("B_MAX", "如 2.0"), "固定 b 与碰撞盘采样共用, 仅气相散射"),
            cell("RANDOM_ORIENT", check("RANDOM_ORIENT"), "默认关"))
        row(cell("PROJECTILE", combo("PROJECTILE", [""]), "扫描后按片段名选"),
            cell("TARGET", combo("TARGET", [""]), "气相双片段必选"))
        row(cell("A_LAT", line("A_LAT", "须与胞公度"), "表面散射必填"),
            cell("N_AIM", mcombo("N_AIM"), "仅表面散射"))
        row(cell("AIM_X", line("AIM_X", "x(POSCAR 帧)"), "落点模式=2 时"),
            cell("AIM_Y", line("AIM_Y", "y"), "落点模式=2 时"))

        # ---- 采样方案(逐片段;扫描后重建;星号=每片段均须选择)----
        group("scheme", "采样方案 *(逐片段)")
        self._sections["scheme"]._btn.setToolTip("写出键: DIST_SCHEME(每片段一词)")
        self._scheme_host = QWidget()
        self._scheme_grid = QGridLayout(self._scheme_host)
        self._cur_lay.addWidget(self._scheme_host)

        # ---- 运行控制(DT 必填在内:可折叠但默认展开)----
        group("run", "运行控制", collapsible=True)
        row(cell("TITLE", line("TITLE", "≤80 字符")))
        row(cell("N_TRAJ", line("N_TRAJ", "1")),
            cell("ISEED", line("ISEED", "1")),
            cell("MAX_STEPS", line("MAX_STEPS", "100000")))
        row(cell("DT", line("DT", "如 0.5"), "必填"),
            cell("INTEGRATOR", combo("INTEGRATOR", [""] + ENUMS["integ"]), "默认 VERLET"),
            cell("ORDER", combo("ORDER", [""] + ORDER_ITEMS["VERLET"]),
                 "VERLET 1/2·SYMPLECTIC 4/6/8·RADAU 1-15"))
        row(cell("DT_MIN", line("DT_MIN", "下界"), "仅 RADAU"),
            cell("DT_MAX", line("DT_MAX", "上界"), "仅 RADAU"))
        row(cell("ELEC_METHOD", line("ELEC_METHOD", "ADIABATIC(≤32)"), "默认 ADIABATIC"),
            cell("GWRITE_LEVEL", combo("GWRITE_LEVEL", ["", "0", "1", "2", "5"]), "默认 2"))
        row(cell("OBSERVABLES", line("OBSERVABLES", "记录列名逗号清单(≤256)")))

        # ---- 浴与谱 ----
        group("bath", "浴与谱", collapsible=True, expanded=False)
        row(cell("BATH", combo("BATH", [""] + ENUMS["bath"]), "默认 NONE"),
            cell("SPECTRUM_SOURCE", combo("SPECTRUM_SOURCE", [""] + ENUMS["spectrum"]),
                 "默认 COMPUTE"))
        row(cell("N_BATH", line("N_BATH", "0"), "ANDERSEN 时"),
            cell("DT_BATH", line("DT_BATH", "0.0"), "N_BATH>0 时须>0"))
        row(cell("T_BATH", line("T_BATH", "K"), "ANDERSEN 必填"),
            cell("NU_COLL", line("NU_COLL", "每步碰撞概率"), "ANDERSEN 必填"))
        row(cell("SPECTRUM_FILE", line("SPECTRUM_FILE", "MANUAL 谱表路径"), "MANUAL 时必填"))

        # ---- 束流参数(入射通道)----
        group("beam", "束流参数(入射通道)", collapsible=True, expanded=False)
        row(cell("N_E_REL", mcombo("N_E_REL")),
            cell("T_TRANS", line("T_TRANS", "K"), "碰撞能模式=0 时须>0"),
            cell("V_WIDTH", line("V_WIDTH", "Å/fs"), "碰撞能模式=2 时须>0"))
        row(cell("N_THTA", mcombo("N_THTA")),
            cell("THTA_MAX", line("THTA_MAX", "rad"), "固定值即采样界"))
        row(cell("N_CHI", mcombo("N_CHI")),
            cell("CHI", line("CHI", "rad"), "方位角模式=0 时"))
        row(cell("N_B", mcombo("N_B")))

        # ---- 成员参数(按方案词打包成行, 每行≤3 键;束流/落点键已有专属控件)----
        group("member", "成员参数(随所选方案显示)", collapsible=True, expanded=False)
        self._member_keys: set[str] = set()
        member_of = {m.key: m for m in MEMBER_TABLE}
        for word, mkeys in SCHEME_MEMBER_KEYS.items():
            batch: list[QWidget] = []
            for k in mkeys:
                mk = member_of[k]
                if mk.key in self._w:  # N_E_REL/N_AIM/AIM_X 等在前述组已有专属控件
                    continue
                # 标签 = hint 剥去尾部括注的物理量名;约束留在占位符
                m = re.search(r"\(([^)]*)\)\s*$", mk.hint)
                phys = (mk.hint[:m.start()] if m else mk.hint).strip()
                if mk.key in MODE_OPTS:  # 二态模式键: 值·含义下拉
                    w = mcombo(mk.key)
                else:
                    w = QLineEdit()
                    w.setPlaceholderText(m.group(1) if m else mk.hint)
                    w.textChanged.connect(self._touch)
                    self._w[mk.key] = w
                self._member_keys.add(mk.key)
                batch.append(cell(mk.key, w, label=phys))
                if len(batch) == 3:
                    row(*batch)
                    batch = []
            if batch:
                row(*batch)

        # ---- 容器键 ----
        g = group("container", "容器键", collapsible=True, expanded=False)
        g._btn.setToolTip("活动容器的 *_keys 词表;reference 预填 K_HARM/Q_TERM/E_SPLIT")
        self._container_table = QTableWidget(0, 2)
        self._container_table.setHorizontalHeaderLabels(["键", "值(空=不写出)"])
        self._container_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        for k in REFERENCE_CONTAINER_KEYS:
            self._container_table.insertRow(self._container_table.rowCount())
            self._container_table.setItem(self._container_table.rowCount() - 1, 0, QTableWidgetItem(k))
            self._container_table.setItem(self._container_table.rowCount() - 1, 1, QTableWidgetItem(""))
        self._container_table.setMinimumHeight(110)
        self._cur_lay.addWidget(self._container_table)

        # ---- 动作 ----
        group("output", "输出")
        act_row = QHBoxLayout()
        act_row.setSpacing(theme.SP_S)
        b_val = QPushButton("校验")
        b_val.setIcon(theme.icon("check"))
        b_val.clicked.connect(self.run_validation)
        b_out = QPushButton("写出 input_qct.txt…")
        b_out.setIcon(theme.icon("save"))
        b_out.setProperty("class", "primary")  # 实心强调主按钮(校验=次级描边)
        b_out.clicked.connect(self._write_file)
        act_row.addWidget(b_val)
        act_row.addWidget(b_out)
        act_row.addStretch(1)
        wrap3 = QWidget()
        wrap3.setLayout(act_row)
        self._cur_lay.addWidget(wrap3)
        self._copy_check = QCheckBox("同时复制 system/ 到输出目录")
        self._cur_lay.addWidget(self._copy_check)  # 独占一行, 不撑组最小宽
        lay.addStretch(1)

        # 右侧:预览(卡片框) + 注记 + 校验结果(语义摘要行 + 结果框)
        right = QWidget()
        rl = QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.setSpacing(theme.SP_S)
        cap = QLabel("起始散射构型(示意——内部几何为平衡坐标,抽样量取代表值)")
        cap.setProperty("class", "caption")
        cap.setWordWrap(True)
        rl.addWidget(cap)
        view_card = QFrame()
        view_card.setProperty("class", "card")
        vcv = QVBoxLayout(view_card)
        vcv.setContentsMargins(theme.SP_S, theme.SP_S, theme.SP_S, theme.SP_S)
        vcv.addWidget(self._view)
        rl.addWidget(view_card, 1)
        rl.addWidget(self._notes)
        lab_v = QLabel("校验结果")
        lab_v.setProperty("class", "caption")
        rl.addWidget(lab_v)
        rl.addWidget(self._verdict_summary)
        rl.addWidget(self._verdict)

        split.addWidget(scroll)
        split.addWidget(right)
        split.setSizes([560, 640])
        outer = QVBoxLayout(self)
        outer.setContentsMargins(theme.SP_L, theme.SP_L, theme.SP_L, theme.SP_L)
        outer.addWidget(split)
        # 无条件必填键的静态星号(DIST_SCHEME 的星号在组标题上)
        for k in ("SYSTEM_DIR", "E_REL", "R_SEP", "DT"):
            self._mark(k, True)
        # 决定者变更即刷新条件显隐;INTEGRATOR 另触发积分阶选项重建
        for k in ("N_AIM", "INTEGRATOR", "BATH", "SPECTRUM_SOURCE",
                  "N_E_REL", "N_CHI", "N_OSC", "N_E_BAR"):
            self._w[k].currentTextChanged.connect(self._refresh_dynamic_rows)
        self._w["INTEGRATOR"].currentTextChanged.connect(self._rebuild_order_items)
        self._refresh_dynamic_rows()  # 初始态:未扫描范式键全藏, 表单最简

    def _mark(self, key: str, on: bool):
        """行标签的必填星号增减(加粗物理量名 + *;单位后缀 [..] 不进粗体)。"""
        lab = self._labels.get(key)
        base = self._label_base.get(key, "")
        if lab is None or not base:
            return
        parts = base.split(" [", 1)
        if on:
            starred = f"<b>{parts[0]} *</b>" + (f" [{parts[1]}" if len(parts) > 1 else "")
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
                c.currentTextChanged.connect(self._refresh_dynamic_rows)
                self._scheme_combos.append(c)
                r = self._scheme_grid.rowCount()
                self._scheme_grid.addWidget(QLabel(f"片段 {i + 1}:{f.name}({f.nat} 原子)"), r, 0)
                self._scheme_grid.addWidget(c, r, 1)
            surf = sd.surface_model == 1
            n = len(sd.frags)
            # 角色下拉(范式合法性由条件显隐统一处理, 词表都给全)
            for key in ("PROJECTILE", "TARGET"):
                c = self._w[key]
                c.blockSignals(True)
                c.clear()
                c.addItems([""] + sd.frag_names)
                c.blockSignals(False)
            # 条件必填星号随范式增减(表面:PROJECTILE+A_LAT;气相多片段:双角色)
            self._mark("PROJECTILE", surf or n > 1)
            self._mark("TARGET", (not surf) and n > 1)
            self._mark("A_LAT", surf)
            para = "表面散射(带胞)" if surf else f"气相散射({n} 片段)"
            box = f"; box {sd.box_lx:.4f}×{sd.box_ly:.4f} Å, skew {sd.skew * 57.2958:.1f}°" if surf else ""
            self._info.setText(f"✓ {para}{box} — 片段: {', '.join(sd.frag_names)}")
        else:
            errs = sd.errors[0] if sd is not None and sd.errors else "未扫描"
            self._info.setText(f"⚠ {errs}")
        self._refresh_dynamic_rows()

    def _refresh_dynamic_rows(self):
        """动态行总装:成员词基显隐 + 条件显隐(两套记账合成行可见性)。"""
        self._refresh_member_visibility()
        self._apply_cond_visibility()

    def _refresh_member_visibility(self):
        """成员参数词基显隐:未选词的成员行隐藏(不清空、可收集——旧语义)。"""
        selected: set[str] = set()
        for c in self._scheme_combos:
            t = c.currentText().strip()
            if t:
                selected.add(t)
        visible = set()
        for w in selected:
            visible.update(SCHEME_MEMBER_KEYS.get(w, []))
        for key in self._member_keys:
            if key in visible:
                self._word_hidden.discard(key)
            else:
                self._word_hidden.add(key)
            self._apply_row_visibility(key)
        # 选中带成员键的方案词后自动展开(否则新显行藏在折叠组里不可见)
        if visible and not self._sections["member"].expanded:
            self._sections["member"].set_expanded(True)

    def _apply_cond_visibility(self):
        """条件显隐:被其余设置判"不需要"的行整行隐藏并清空(collect 跳过,绝不写出)。
        规则权威 = core.validate 移植的 Fortran 规则;积分器×阶域 = propagator.f90。"""
        sd = self._sysdef
        surf = sd is not None and sd.ok and sd.surface_model == 1
        n = len(sd.frags) if sd is not None and sd.ok else 0

        def cur(k: str) -> str:
            w = self._w.get(k)
            if not isinstance(w, QComboBox):
                return ""
            d = w.currentData()
            return str(d if d is not None else w.currentText()).strip()

        show: dict[str, bool] = {}

        def want(k: str, on: bool):
            show[k] = show.get(k, True) and on

        # 范式窗口(表面:弹丸+晶格+落点;气相对:双角色+碰撞参数;单片段:仅碰撞参数)
        want("PROJECTILE", surf or n > 1)
        want("TARGET", (not surf) and n > 1)
        want("B_MAX", not surf)
        want("A_LAT", surf)
        want("N_AIM", surf)
        # 落点模式 2(定点)才需要坐标
        want("AIM_X", surf and cur("N_AIM") == "2")
        want("AIM_Y", surf and cur("N_AIM") == "2")
        # 自适应步长界仅 RADAU(dt_min 只被 radau_step 消费)
        radau = cur("INTEGRATOR") == "RADAU"
        want("DT_MIN", radau)
        want("DT_MAX", radau)
        # 热浴 ANDERSEN 才有浴参数
        ba = cur("BATH") == "ANDERSEN"
        for k in ("N_BATH", "DT_BATH", "T_BATH", "NU_COLL"):
            want(k, ba)
        # 谱来源 MANUAL 才要谱表文件
        want("SPECTRUM_FILE", cur("SPECTRUM_SOURCE") == "MANUAL")
        # 束流:碰撞能 0热→平动温度/2加权→速度展宽;方位角采样→不需固定值
        want("T_TRANS", cur("N_E_REL") == "0")
        want("V_WIDTH", cur("N_E_REL") == "2")
        want("CHI", cur("N_CHI") != "1")
        # 成员内二态模式(振子 0热→温度/1定态→能级;反应坐标 1热→温度)
        want("T_OSC", cur("N_OSC") != "1")
        want("N_LEVEL", cur("N_OSC") == "1")
        want("T_BAR", cur("N_E_BAR") == "1")

        for k, on in show.items():
            if on:
                self._cond_hidden.discard(k)
            elif k not in self._cond_hidden:
                self._cond_hidden.add(k)
                w = self._w.get(k)  # 隐藏即清空(与范式切换清空同规约)
                if isinstance(w, QLineEdit):
                    w.clear()
                elif isinstance(w, QComboBox):
                    w.setCurrentIndex(0)
            self._apply_row_visibility(k)

    def _apply_row_visibility(self, key: str):
        """单元对可见性 = 条件显隐 ∧ 词基显隐;同行其余键左移补位;
        行容器仅当整行键全隐藏时塌缩。"""
        on = key not in self._cond_hidden and key not in self._word_hidden
        c = self._cell_of_key.get(key)
        if c is not None:
            c.setVisible(on)
        for r, keys in self._rows:
            r.setVisible(any(k not in self._cond_hidden and k not in self._word_hidden
                             for k in keys))

    def row_hidden(self, key: str) -> bool:
        """键的当前可见性(单元对或所在行容器隐藏即隐藏)。"""
        c = self._cell_of_key.get(key)
        return c is None or c.isHidden() or c.parentWidget().isHidden()

    def _same_row(self, a: str, b: str) -> bool:
        """两键是否打包在同一行(行布局契约的测试入口)。"""
        return any(a in keys and b in keys for _, keys in self._rows)

    def set_choice(self, key: str, value: str):
        """按下拉的数据值选择(模式类下拉显示为"值 · 含义", 显示文本非值)。"""
        c = self._w[key]
        for i in range(c.count()):
            if c.itemData(i) == value:
                c.setCurrentIndex(i)
                return

    def _rebuild_order_items(self):
        """积分阶选项随积分器重建(VERLET 1/2·SYMPLECTIC 4/6/8·RADAU 1-15)。"""
        integ = self._w["INTEGRATOR"].currentText().strip() or "VERLET"
        items = [""] + ORDER_ITEMS.get(integ, ORDER_ITEMS["VERLET"])
        c = self._w["ORDER"]
        keep = c.currentData()
        c.blockSignals(True)
        c.clear()
        for it in items:
            c.addItem(it, it)
        c.setCurrentText(keep if keep in items else "")
        c.blockSignals(False)

    # ---------- 数据汇聚 ----------
    def collect(self) -> tuple[dict[str, str], dict[str, str]]:
        values: dict[str, str] = {}
        for k, w in self._w.items():
            if k in self._cond_hidden:  # 条件隐藏的键绝不写出(隐藏即清空, 此为双保险)
                continue
            if isinstance(w, QLineEdit):
                t = w.text().strip()
                if t:
                    values[k] = t
            else:  # QComboBox: 取数据值(模式类下拉显示"值 · 含义"而数据为原值)
                d = w.currentData()
                t = str(d if d is not None else w.currentText()).strip()
                if t and k != "PROJECTILE" and k != "TARGET":
                    values[k] = t
        # 角色键按原文(区分大小写按文件名匹配,大小写不敏感)
        for k in ("PROJECTILE", "TARGET"):
            if k in self._cond_hidden:
                continue
            t = self._w[k].currentText().strip()
            if t:
                values[k] = t
        # 复选框
        for k, c in self._checks.items():
            values[k] = "T" if c.isChecked() else "F"
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
            self._set_verdict_state("ok", 0)
        else:
            self._verdict.setPlainText("\n".join(f"✗ {i}" for i in issues))
            self._set_verdict_state("error", len(issues))

    def _set_verdict_state(self, state: str, n: int):
        """校验语义色: 结果框边框/底色与摘要行按 ok/error 变色(QSS; 文本不变)。"""
        self._verdict.setProperty("state", state)
        self._verdict_summary.setProperty("state", state)
        self._verdict_summary.setText("● 校验通过" if state == "ok" else f"● {n} 处问题")
        for w in (self._verdict, self._verdict_summary):
            w.style().unpolish(w)  # 属性变更后重刷样式
            w.style().polish(w)

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
