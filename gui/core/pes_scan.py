# -*- coding: utf-8 -*-
"""pes_scan.py - Tier 1 纯势容器十名扫描(scripts/prebuild.sh:85-168 的 Python 移植)

设计:
  auto_scan_file 的 awk 语义逐条移植:剥 '!' 注释尾→去行首空白→小写副本
  匹配→(subroutine|function) 关键字可在前缀后("logical function term(sta)"
  惯用形)→'end' 前缀排除("end subroutine pes" 不算)→实参括号剥除不强求
  (无括号定义也计数——宁可响亮报错不静默错织)→模块归属 = 定义前最近的
  module 语句(排除 "module procedure")。八条有名中止的判定顺序与文案
  基调与 prebuild 一致(文案换成中文陈述,判据键一致)。
  附加提取(非契约,供 GUI 展示):keys() 词表字面量、method_name() 返回
  字面量、模块级 real(8) 公共变量清单(生成器表单预填用)。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# 十名与期望形(prebuild.sh:118-119, 142-143)
SUB_NAMES = ("pes", "load", "prop", "match", "sync", "init")         # 期望 subroutine
FUNC_NAMES = ("keys", "term", "classify", "method_name")             # 期望 function
ALL_NAMES = ("pes", "keys", "load", "term", "classify",
             "prop", "match", "sync", "init", "method_name")         # 核对顺序

_MODULE_RE = re.compile(r"^module\s+[a-z_]")
_DEF_RE = re.compile(r"(subroutine|function)\s+[a-z_]")
_END_TAIL_RE = re.compile(r"end\s*$")


@dataclass
class NameDef:
    """一个预留名的定义(一次出现)。"""
    kind: str      # 's'(subroutine) / 'f'(function)
    name: str
    module: str    # 归属模块名(原始大小写;空 = 模块外)
    file: str


@dataclass
class ScanResult:
    """容器文件夹扫描结论。"""
    folder: str
    is_tier1: bool = True          # False = 含 reg_*.f90(Tier 2,非本检查器对象)
    has_f90: bool = False
    defs: dict[str, list[NameDef]] = field(default_factory=dict)   # 名 → 出现列表
    errors: list[str] = field(default_factory=list)                # 中止级问题(中文)
    notes: list[str] = field(default_factory=list)
    keys_words: list[str] = field(default_factory=list)            # 提取:keys() 词表
    method_name_value: str = ""                                    # 提取:注册名
    module_vars: dict[str, list[str]] = field(default_factory=list)  # 提取:模块 real(8) 公共变量名

    def found(self, name: str) -> NameDef | None:
        lst = self.defs.get(name, [])
        return lst[0] if lst else None

    @property
    def ok(self) -> bool:
        return self.is_tier1 and self.has_f90 and not self.errors and self.found("pes") is not None

    @property
    def prop_match_pair(self) -> bool:
        return self.found("prop") is not None and self.found("match") is not None


def scan_source(text: str, fname: str = "<src>") -> list[tuple[str, str, str]]:
    """一个源文本 → [(kind, name, module)](prebuild auto_scan_file 的 awk 移植)。"""
    out: list[tuple[str, str, str]] = []
    module = ""
    for raw in text.splitlines():
        line = re.sub(r"!.*", "", raw)          # 剥注释尾(含字符串里 '!' 的 awk 同款误剥)
        line = line.lstrip()
        lline = line.lower()
        if _MODULE_RE.match(lline):
            words = line.split()
            if len(words) >= 2 and words[1].lower() != "procedure":
                module = words[1]
        m = _DEF_RE.search(lline)
        if not m:
            continue
        if m.start() > 0 and _END_TAIL_RE.search(lline[: m.start()]):
            continue                             # "end subroutine pes" 排除
        words = lline[m.start():].split()
        kw, nm = words[0], words[1]
        nm = re.sub(r"\(.*", "", nm)             # 实参括号剥除不强求
        if nm in SUB_NAMES or nm in FUNC_NAMES:
            out.append(("s" if kw == "subroutine" else "f", nm, module))
    return out


def _extract_keys_words(folder: Path, defs: dict[str, list[NameDef]]) -> list[str]:
    """keys() 体里的字符串字面量(展示用;非契约检查)。"""
    d = defs.get("keys", [None])[0] if defs.get("keys") else None
    if d is None:
        return []
    text = (folder / d.file).read_text(encoding="utf-8", errors="replace")
    m = re.search(r"function\s+keys\s*\(\s*\)(.*?)(?:end function|$)", text, re.S | re.I)
    if not m:
        return []
    return [w for w in re.findall(r"'([A-Za-z0-9_]+)'", m.group(1))]


def _extract_method_name(folder: Path, defs: dict[str, list[NameDef]]) -> str:
    d = defs.get("method_name", [None])[0] if defs.get("method_name") else None
    if d is None:
        return ""
    text = (folder / d.file).read_text(encoding="utf-8", errors="replace")
    m = re.search(r"function\s+method_name\s*\(\s*\)(.*?)(?:end function|$)", text, re.S | re.I)
    if not m:
        return ""
    lits = re.findall(r"'([A-Za-z0-9_.\-]+)'", m.group(1))
    return lits[0] if lits else ""


_VAR_RE = re.compile(
    r"^\s*real\s*\(\s*8\s*\)\s*::\s*([A-Za-z]\w*(?:\s*,\s*[A-Za-z]\w*)*)"
    r"(\s*=\s*[-+0-9.eEdD]+)?\s*$"
)


def _extract_module_vars(folder: Path) -> dict[str, list[str]]:
    """各模块 real(8) 标量变量名(带初始化者优先)——生成器表单预填用。"""
    out: dict[str, list[str]] = {}
    for f in sorted(folder.glob("*.f90")):
        module = ""
        for raw in f.read_text(encoding="utf-8", errors="replace").splitlines():
            line = re.sub(r"!.*", "", raw)
            low = line.lstrip().lower()
            if _MODULE_RE.match(low):
                words = line.split()
                if len(words) >= 2 and words[1].lower() != "procedure":
                    module = words[1]
                continue
            m = _VAR_RE.match(line)
            if m and module:
                names = [t.strip() for t in m.group(1).split(",")]
                out.setdefault(module, []).extend(names)
    return out


def scan_container(folder: str | Path) -> ScanResult:
    """扫描容器文件夹:候选判定 → 十名扫描 → 八中止(镜像 prebuild 顺序)。"""
    p = Path(folder)
    res = ScanResult(folder=str(p))
    if not p.is_dir():
        res.errors.append(f"文件夹不存在: {p}")
        return res
    f90s = sorted(p.glob("*.f90"))
    res.has_f90 = bool(f90s)
    if not f90s:
        res.errors.append("文件夹无 .f90 文件(纯势容器须暴露预留名,势能面组至少 pes)")
        return res
    if any(f.name.startswith("reg_") for f in f90s):
        res.is_tier1 = False
        res.notes.append("含 reg_*.f90:显式 reg 容器(Tier 2),reg 是权威,十名全忽略")
        return res

    defs: dict[str, list[NameDef]] = {n: [] for n in ALL_NAMES}
    for f in f90s:
        for kind, nm, module in scan_source(f.read_text(encoding="utf-8", errors="replace"), f.name):
            defs[nm].append(NameDef(kind, nm, module, f.name))
    res.defs = defs

    # 逐名核对:计数/形/归属(prebuild.sh:142-152)
    for nm in ALL_NAMES:
        want = "s" if nm in SUB_NAMES else "f"
        got = defs[nm]
        if not got:
            continue
        if len(got) >= 2:
            res.errors.append(f"预留名 '{nm}' 定义了 {len(got)} 次(只允许一次)")
            continue
        d = got[0]
        if d.kind != want:
            res.errors.append(
                f"预留名 '{nm}' 须为 {'subroutine' if want == 's' else 'function'},"
                f"实际是 {'subroutine' if d.kind == 's' else 'function'}"
            )
            continue
        if not d.module:
            res.errors.append(f"预留名 '{nm}' 定义在模块之外(织入需要模块归属)")
            continue

    # 配对规则(prebuild.sh:153-166)
    n_pm = sum(1 for n in ("prop", "match") if defs[n])
    if n_pm == 1:
        res.errors.append("prop 与 match 须成对定义(半个包不构成方法)")
    if n_pm == 2:
        if not defs["method_name"]:
            res.errors.append("prop+match 已定义但缺 method_name(注册名是作者数据)")
    else:
        if defs["sync"]:
            res.errors.append("sync 已定义但 prop+match 未定义(孤立的 sync 不构成包)")
        if defs["init"]:
            res.errors.append("init 已定义但 prop+match 未定义(孤立的 init 不构成包)")
        if defs["method_name"]:
            res.errors.append("method_name 已定义但 prop+match 未定义")
    if not defs["pes"]:
        res.errors.append(
            "容器已安装(文件夹在场且无 reg_*.f90)却什么都没接:请暴露恰好一个 'pes' 入口,"
            "或写 reg_*.f90 完整容器"
        )

    # 展示性提取
    if not res.errors:
        res.keys_words = _extract_keys_words(p, defs)
        res.method_name_value = _extract_method_name(p, defs)
    res.module_vars = _extract_module_vars(p)
    return res
