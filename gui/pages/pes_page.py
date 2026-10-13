# -*- coding: utf-8 -*-
"""pes_page.py - 页2:自定义 PES 容器检查 → 程序补全 → 隔离树编译

设计:
  三段流:①[扫描] core.pes_scan 十名+八中止,结果表逐名判定,提取 keys
  词表/方法名/模块变量;②[生成] core.pes_gen 按设定补全可选槽位(term 三变体/
  classify 两变体/keys+load 接线),方法组五名不生成(真实物理,界内说明);
  ③[编译验证] core.sandbox 隔离副本编译(QThread 后台流式日志落盘),
  参照容器与被检容器并存对编译合法,绝不触碰工作树。
"""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import QThread, pyqtSignal, Qt
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QGridLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox,
    QPlainTextEdit, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QScrollArea, QVBoxLayout, QWidget,
)

from ..core import pes_gen, pes_scan, sandbox
from ..core.pes_scan import ALL_NAMES
from .. import theme

_WANT = {n: ("subroutine" if n in pes_scan.SUB_NAMES else "function") for n in ALL_NAMES}


class _CompileWorker(QThread):
    line = pyqtSignal(str)
    done = pyqtSignal(object)

    def __init__(self, container: Path, fc: str):
        super().__init__()
        self.container = container
        self.fc = fc
        self.result = None

    def run(self):
        try:
            res = sandbox.run_isolated_compile(
                self.container, fc=self.fc, on_line=self.line.emit
            )
            self.result = res
            self.done.emit(res)
        except Exception as exc:
            self.done.emit(exc)


class PesPage(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: pes_scan.ScanResult | None = None
        self._worker: _CompileWorker | None = None
        self._build_ui()

    # ---------- UI ----------
    def _build_ui(self):
        # 三段内容共处滚动区(最小尺寸塌缩——不把整页内容顶成窗口最小高)
        outer = QVBoxLayout(self)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        host = QWidget()
        outer_v = QVBoxLayout(host)
        outer_v.setContentsMargins(0, 0, 0, 0)
        scroll.setWidget(host)
        outer.addWidget(scroll)

        # ① 扫描(组标题括注细节入组悬停——2026-10-13 文案精简裁决)
        g1 = QGroupBox("① 容器扫描")
        g1.setToolTip("纯势容器契约:十名 + 八中止。\n候选通常在 systems/inactive/ 下,亦可任选路径。")
        l1 = QGridLayout(g1)
        self._dir = QLineEdit()
        self._dir.setPlaceholderText("容器文件夹路径")
        b_br = QPushButton("浏览…")
        b_br.setIcon(theme.icon("folder"))
        b_br.clicked.connect(self._browse)
        b_scan = QPushButton("扫描")
        b_scan.setIcon(theme.icon("search"))
        b_scan.clicked.connect(self.scan)
        l1.addWidget(QLabel("容器文件夹"), 0, 0)
        l1.addWidget(self._dir, 0, 1)
        l1.addWidget(b_br, 0, 2)
        l1.addWidget(b_scan, 0, 3)
        self._quick = QComboBox()
        self._quick.addItem("(请选择)")
        self._quick.setToolTip("列出 systems/inactive/ 下的候选容器,选中即扫描")
        self._fill_quick()
        self._quick.currentIndexChanged.connect(self._quick_pick)
        l1.addWidget(QLabel("快选"), 1, 0)
        l1.addWidget(self._quick, 1, 1, 1, 3)
        self._table = QTableWidget(0, 5)
        self._table.setHorizontalHeaderLabels(["预留名", "期望形", "在场", "实际", "归属模块/文件"])
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._table.setMinimumHeight(240)
        l1.addWidget(self._table, 2, 0, 1, 4)
        self._scan_info = QLabel("")
        self._scan_info.setWordWrap(True)
        l1.addWidget(self._scan_info, 3, 0, 1, 4)
        outer_v.addWidget(g1)

        # ② 生成器
        g2 = QGroupBox("② 程序补全")
        g2.setToolTip("可行模板;整文件覆盖式 auto_slots.f90,删除即撤销。")
        l2 = QGridLayout(g2)
        self._gen_term = QCheckBox("term 终止判定")
        self._term_var = QComboBox()
        self._term_var.addItems(list(pes_gen.TERM_VARIANTS))
        self._term_r = QDoubleSpinBox()
        self._term_r.setRange(0.1, 1e5)
        self._term_r.setDecimals(3)
        self._term_r.setValue(5.0)
        self._term_r.setSuffix(" Å")
        l2.addWidget(self._gen_term, 0, 0)
        l2.addWidget(self._term_var, 0, 1)
        l2.addWidget(QLabel("R"), 0, 2)
        l2.addWidget(self._term_r, 0, 3)
        l2.addWidget(QLabel("变体说明见下拉悬停"), 0, 4)

        self._gen_cls = QCheckBox("classify 结果分类")
        self._cls_var = QComboBox()
        self._cls_var.addItems(list(pes_gen.CLASSIFY_VARIANTS))
        self._cls_r = QDoubleSpinBox()
        self._cls_r.setRange(0.1, 1e5)
        self._cls_r.setDecimals(3)
        self._cls_r.setValue(3.0)
        self._cls_r.setSuffix(" Å")
        self._cls_label = QLineEdit("bound")
        l2.addWidget(self._gen_cls, 1, 0)
        l2.addWidget(self._cls_var, 1, 1)
        l2.addWidget(QLabel("R"), 1, 2)
        l2.addWidget(self._cls_r, 1, 3)
        l2.addWidget(self._cls_label, 1, 4)

        self._gen_kl = QCheckBox("keys+load 参数接线")
        self._kl_mod = QComboBox()
        self._kl_mod.setEditable(True)
        self._kl_table = QTableWidget(0, 2)
        self._kl_table.setHorizontalHeaderLabels(["容器键名", "目标 real(8) 变量"])
        self._kl_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._kl_table.setMinimumHeight(100)
        b_add = QPushButton("+ 绑定")
        b_add.clicked.connect(lambda: self._kl_table.insertRow(self._kl_table.rowCount()))
        b_del = QPushButton("− 删行")
        b_del.clicked.connect(lambda: self._kl_table.removeRow(self._kl_table.currentRow()))
        l2.addWidget(self._gen_kl, 2, 0)
        l2.addWidget(QLabel("目标模块"), 2, 1)
        l2.addWidget(self._kl_mod, 2, 2)
        l2.addWidget(b_add, 2, 3)
        l2.addWidget(b_del, 2, 4)
        l2.addWidget(self._kl_table, 3, 0, 1, 5)
        note = QLabel("方法组五名不生成(真实物理,作者提供);keys+load 需作者模块的公共 real(8) 变量,扫描后预填。")
        note.setToolTip(
            "未提供时绝热面缺省即可跑。keys+load 需作者模块的公共 real(8) 变量(扫描后自动预填)。"
        )
        note.setWordWrap(True)
        l2.addWidget(note, 4, 0, 1, 5)
        b_gen = QPushButton("生成/覆盖 auto_slots.f90")
        b_gen.setIcon(theme.icon("wand"))
        b_gen.clicked.connect(self.generate)
        l2.addWidget(b_gen, 5, 0, 1, 2)
        outer_v.addWidget(g2)

        # ③ 编译
        g3 = QGroupBox("③ 编译验证")
        g3.setToolTip("隔离副本编译,不动工作树;参照容器与被检容器并存对编译合法。")
        l3 = QGridLayout(g3)
        self._fc = QComboBox()
        self._fc.addItems(["gfortran", "ifx"])
        b_mk = QPushButton("开始编译")
        b_mk.setIcon(theme.icon("hammer"))
        b_mk.setProperty("class", "primary")  # 本页主动作(实心强调)
        b_mk.clicked.connect(self.compile)
        self._b_clean = QPushButton("清理沙盒")
        # 判词走 Task 3 简洁文案; 日志随沙盒删除的细节入悬停(build.log 写在沙盒内)
        self._b_clean.setToolTip("删除隔离沙盒副本;编译日志文件一并删除")
        self._b_clean.clicked.connect(self._cleanup)
        self._b_clean.setEnabled(False)
        l3.addWidget(QLabel("编译器"), 0, 0)
        l3.addWidget(self._fc, 0, 1)
        l3.addWidget(b_mk, 0, 2)
        l3.addWidget(self._b_clean, 0, 3)
        self._log = QPlainTextEdit(readOnly=True)
        self._log.setMinimumHeight(220)
        self._log.setProperty("class", "log")  # 等宽字体由全局 QSS 收口
        l3.addWidget(self._log, 1, 0, 1, 4)
        self._mk_verdict = QLabel("")
        self._mk_verdict.setWordWrap(True)
        l3.addWidget(self._mk_verdict, 2, 0, 1, 4)
        outer_v.addWidget(g3)
        outer_v.addStretch(1)

    def _fill_quick(self):
        try:
            root = sandbox.repo_root() / "systems" / "inactive"
            for d in sorted(root.iterdir()):
                if d.is_dir():
                    self._quick.addItem(d.name)
        except OSError:
            pass

    def _quick_pick(self, idx: int):
        if idx <= 0:
            return
        name = self._quick.itemText(idx)
        self._dir.setText(str(sandbox.repo_root() / "systems" / "inactive" / name))
        self.scan()

    # ---------- ① 扫描 ----------
    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "选择容器文件夹")
        if d:
            self._dir.setText(d)
            self.scan()

    def scan(self):
        path = self._dir.text().strip()
        if not path:
            return
        res = pes_scan.scan_container(path)
        self._scan = res
        self._table.setRowCount(len(ALL_NAMES))
        for i, name in enumerate(ALL_NAMES):
            d = res.found(name)
            self._table.setItem(i, 0, QTableWidgetItem(name))
            self._table.setItem(i, 1, QTableWidgetItem(_WANT[name]))
            if d is None:
                it = QTableWidgetItem("—")
                it.setForeground(Qt.GlobalColor.gray)
                self._table.setItem(i, 2, it)
                self._table.setItem(i, 3, QTableWidgetItem("(缺)"))
                self._table.setItem(i, 4, QTableWidgetItem(""))
            else:
                it = QTableWidgetItem("✓")
                it.setForeground(Qt.GlobalColor.darkGreen)
                self._table.setItem(i, 2, it)
                self._table.setItem(i, 3, QTableWidgetItem(
                    "subroutine" if d.kind == "s" else "function"))
                self._table.setItem(i, 4, QTableWidgetItem(f"{d.module} ({d.file})"))
        if not res.is_tier1:
            self._scan_info.setText("ℹ 显式注册容器:reg 权威,十名忽略。")
        elif res.ok:
            extras = []
            if res.keys_words:
                extras.append(f"keys 词表: {', '.join(res.keys_words)}")
            if res.method_name_value:
                extras.append(f"方法名: {res.method_name_value}")
            self._scan_info.setText("✓ 契约满足。" + (";".join(extras)))
        else:
            self._scan_info.setText("✗ " + " / ".join(res.errors))
        # keys+load 表单预填:目标模块 = pes 归属模块,变量 = 模块 real(8) 清单
        self._kl_mod.clear()
        if res.is_tier1:
            for mod in sorted(res.module_vars):
                self._kl_mod.addItem(mod)
            pes_mod = res.found("pes").module if res.found("pes") else ""
            if pes_mod:
                self._kl_mod.setCurrentText(pes_mod)
        # 变体 tooltip
        for i, (k, v) in enumerate(pes_gen.TERM_VARIANTS.items()):
            self._term_var.setItemData(i, v, Qt.ItemDataRole.ToolTipRole)
        for i, (k, v) in enumerate(pes_gen.CLASSIFY_VARIANTS.items()):
            self._cls_var.setItemData(i, v, Qt.ItemDataRole.ToolTipRole)

    # ---------- ② 生成 ----------
    def generate(self):
        if self._scan is None or not self._scan.is_tier1:
            QMessageBox.warning(self, "先扫描", "请先扫描一个纯势容器候选文件夹。")
            return
        path = Path(self._dir.text().strip())
        term = cls = kl = None
        if self._gen_term.isChecked():
            term = pes_gen.TermSpec(self._term_var.currentText(), self._term_r.value())
        if self._gen_cls.isChecked():
            var = self._cls_var.currentText()
            cls = pes_gen.ClassifySpec(var, r=self._cls_r.value(), label=self._cls_label.text())
        if self._gen_kl.isChecked():
            mod = self._kl_mod.currentText().strip()
            if not mod:
                QMessageBox.warning(self, "缺目标模块", "keys+load 需要目标模块名。")
                return
            binds = []
            for r in range(self._kl_table.rowCount()):
                k = self._kl_table.item(r, 0)
                v = self._kl_table.item(r, 1)
                k = k.text().strip().upper() if k else ""
                v = v.text().strip() if v else ""
                if k and v:
                    binds.append((k, v))
            if not binds:
                QMessageBox.warning(self, "缺绑定", "keys+load 至少需要一条键-变量绑定。")
                return
            kl = pes_gen.KeysLoadSpec(module=mod, bindings=binds)
        if term is None and cls is None and kl is None:
            QMessageBox.information(self, "未选择", "请勾选至少一个要生成的槽位。")
            return
        try:
            p = pes_gen.generate_auto_slots(path, term=term, classify=cls, keys_load=kl)
        except ValueError as exc:
            QMessageBox.critical(self, "生成被拒", str(exc))
            return
        QMessageBox.information(self, "已生成", f"已写出 {p}\n自动重扫中…")
        self.scan()

    # ---------- ③ 编译 ----------
    def compile(self):
        if self._worker is not None and self._worker.isRunning():
            return
        path = Path(self._dir.text().strip())
        if not path.is_dir():
            QMessageBox.warning(self, "路径无效", "容器文件夹不存在。")
            return
        self._log.clear()
        self._mk_verdict.setText("编译中…(隔离副本)")
        self._worker = _CompileWorker(path, self._fc.currentText())
        self._worker.line.connect(self._log.appendPlainText)
        self._worker.done.connect(self._compile_done)
        self._worker.start()

    def _compile_done(self, res):
        if isinstance(res, Exception):
            self._mk_verdict.setText(f"✗ 启动失败:{res}")
            return
        if res.ok:
            self._mk_verdict.setText(
                f"✓ 编译通过。二进制:{res.sandbox}/src/build/{self._fc.currentText()}/venus.e;"
                f"日志:{res.log}"
            )
        else:
            self._mk_verdict.setText(f"✗ 编译失败(exit {res.returncode})。日志:{res.log};沙盒保留供诊断")
        self._b_clean.setEnabled(True)
        self._last = res

    def _cleanup(self):
        res = getattr(self, "_last", None)
        if res is not None:
            sandbox.cleanup(res)
            self._b_clean.setEnabled(False)
            self._mk_verdict.setText("沙盒已清理。")
