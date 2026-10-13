# -*- coding: utf-8 -*-
"""sysdir.py - system/ 体系文件夹解析(src/utils/specio.f90 + interface/sysdef_interface.f90 的移植)

设计:
  .xyz / POSCAR 解析、元素质量表、名序折叠排序、单胞规则(散射类自 SCATTERING 键取,2026-10-13 裁决)
  逐条移植自 Fortran 侧;坐标全程笛卡尔 Å
  (POSCAR Direct 分数坐标解析时一次转笛卡尔,速度块读跳)。
  质量表数值逐字取自 specio.f90 的 DATA 表(程序物理常数,不得重推)。
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

# 元素质量表 [amu](specio.f90 tbl_symb/tbl_mass 逐字转录;D 同位素 2.014 是别名)
ELEM_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
ELEM_MASSES = [
    1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011, 14.007, 15.999,
    18.998403163, 20.1797, 22.98976928, 24.305, 26.9815385, 28.085,
    30.973761998, 32.06, 35.45, 39.948, 39.0983, 40.078, 44.955908,
    47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934,
    63.546, 65.38, 69.723, 72.630, 74.921595, 78.971, 79.904, 83.798,
    85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, 98.0, 101.07,
    102.90550, 106.42, 107.8682, 112.414, 114.818, 118.710, 121.760,
    127.60, 126.90447, 131.293, 132.90545196, 137.327, 138.90547,
    140.116, 140.90766, 144.242, 145.0, 150.36, 151.964, 157.25,
    158.92535, 162.500, 164.93033, 167.259, 168.93422, 173.045,
    174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.966569, 200.592, 204.38, 207.2, 208.98040, 209.0,
    210.0, 222.0, 223.0, 226.0, 227.0, 232.0377, 231.03588, 238.02891,
    237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0,
    259.0, 266.0, 267.0, 268.0, 269.0, 270.0, 269.0, 278.0, 281.0,
    282.0, 285.0, 286.0, 289.0, 290.0, 293.0, 294.0, 294.0,
]
ELEM_MASS = dict(zip(ELEM_SYMBOLS, ELEM_MASSES))
ELEM_MASS["D"] = 2.014  # 氘别名(specio.f90 elem_mass)


def elem_mass(symbol: str) -> float:
    """元素符号 → 质量 [amu](大小写规范化;未知名抛 ValueError)。"""
    s = symbol.strip()
    if not s:
        raise ValueError("空元素符号")
    norm = s[0].upper() + s[1:].lower()
    if norm not in ELEM_MASS:
        raise ValueError(f"未知元素符号: {symbol}")
    return ELEM_MASS[norm]


@dataclass
class Frag:
    """一个反应物(=一个体系文件)。"""
    file: str                 # 文件名
    name: str                 # 反应物名(扩展名剥离)
    symbols: list[str]        # 每原子元素符号
    coords: list[list[float]]  # 每原子笛卡尔坐标 [Å] (n,3)
    cell: list[list[float]] | None = None  # 晶胞 3x3(POSCAR 才有)

    @property
    def nat(self) -> int:
        return len(self.symbols)

    @property
    def masses(self) -> list[float]:
        return [elem_mass(s) for s in self.symbols]


@dataclass
class SysDef:
    """system/ 文件夹解析产物(sysdef buffer + list_atoms 派生量的合并视图)。"""
    folder: str
    frags: list[Frag] = field(default_factory=list)
    has_cell: bool = False                   # 扫描侧事实:任一 POSCAR 带胞(信息性;
                                              # 散射类自 SCATTERING 键取,2026-10-13 裁决)
    cell: list[list[float]] | None = None   # 3x3 晶矢 [Å](行存;扫描信息,落点框以 CELL 键为准)
    box_lx: float = 0.0                     # 胞 x 边长 [Å](config_atoms 余弦定理;扫描派生)
    box_ly: float = 0.0
    skew: float = 0.0                       # 胞斜角 [rad]
    errors: list[str] = field(default_factory=list)

    @property
    def frag_names(self) -> list[str]:
        return [f.name for f in self.frags]

    def frag_by_name(self, name: str) -> Frag | None:
        low = name.strip().lower()
        for f in self.frags:
            if f.name.lower() == low:
                return f
        return None

    @property
    def ok(self) -> bool:
        return not self.errors and bool(self.frags)


def _tokens(line: str) -> list[str]:
    return line.replace("\t", " ").split()


def parse_xyz(text: str, fname: str = "<xyz>") -> Frag:
    """解析 .xyz:原子数行 + 注释行 + "符号 x y z" 行(specio read_xyz 移植)。"""
    lines = text.splitlines()
    if not lines:
        raise ValueError(f"{fname}: 空文件(缺原子数行)")
    try:
        n = int(_tokens(lines[0])[0])
    except (ValueError, IndexError):
        raise ValueError(f"{fname}: 原子数行不是整数: {lines[0]}") from None
    if n <= 0:
        raise ValueError(f"{fname}: 原子数非正")
    if len(lines) < 2:
        raise ValueError(f"{fname}: 缺注释行")
    symbols: list[str] = []
    coords: list[list[float]] = []
    for i in range(n):
        row = 2 + i
        if row >= len(lines):
            raise ValueError(f"{fname}: 第 {row + 1} 行缺失(坐标块提前结束)")
        toks = _tokens(lines[row])
        if len(toks) < 4:
            raise ValueError(f"{fname}: 第 {row + 1} 行过短: {lines[row]}")
        sym = toks[0]
        if len(sym) > 4:
            raise ValueError(f"{fname}: 第 {row + 1} 行符号超 4 字符: {sym}")
        try:
            xyz = [float(toks[1]), float(toks[2]), float(toks[3])]
        except ValueError:
            raise ValueError(f"{fname}: 第 {row + 1} 行坐标非数值: {lines[row]}") from None
        symbols.append(sym)
        coords.append(xyz)
    name = Path(fname).stem
    return Frag(file=fname, name=name, symbols=symbols, coords=coords)


def parse_poscar(text: str, fname: str = "<poscar>") -> tuple[Frag, list[list[float]]]:
    """解析 POSCAR(不带胞体系文件)。返回 (Frag, cell 3x3)。

    specio read_poscar 移植:注释行/普适缩放因子(负值=目标体积规则)/三条晶矢
    /元素行+计数行/[Selective dynamics]/Direct|Cartesian/坐标行;尾部速度块
    读跳不解释;Direct 分数坐标 = fr · cell(matmul(fr, cell),行存晶矢线性组合)。
    """
    lines = text.splitlines()

    def need(i: int, what: str) -> str:
        if i >= len(lines):
            raise ValueError(f"{fname}: 缺{what}(第 {i + 1} 行之后文件结束)")
        return lines[i]

    # 1. 缩放因子(负值 = 目标体积规则)
    try:
        s = float(_tokens(need(1, "缩放因子行"))[0])
    except (ValueError, IndexError):
        raise ValueError(f"{fname}: 缩放因子行非数值: {lines[1] if len(lines) > 1 else ''}") from None
    if s == 0.0:
        raise ValueError(f"{fname}: 缩放因子为零")

    # 2. 三条晶矢(0 基行 2-4:注释 0/缩放 1/晶矢 2,3,4)
    vec = []
    for i in range(3):
        toks = _tokens(need(2 + i, f"晶矢行 {i + 1}"))
        if len(toks) < 3:
            raise ValueError(f"{fname}: 晶矢行 {i + 1} 过短: {lines[2 + i]}")
        try:
            vec.append([float(t) for t in toks[:3]])
        except ValueError:
            raise ValueError(f"{fname}: 晶矢行 {i + 1} 非数值: {lines[2 + i]}") from None
    if s < 0.0:
        det = (vec[0][0] * (vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1])
               - vec[0][1] * (vec[1][0] * vec[2][2] - vec[1][2] * vec[2][0])
               + vec[0][2] * (vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0]))
        if det == 0.0:
            raise ValueError(f"{fname}: 晶矢奇异(负缩放因子=目标体积)")
        s_eff = (abs(s) / abs(det)) ** (1.0 / 3.0)
    else:
        s_eff = s
    cell = [[s_eff * vec[r][c] for c in range(3)] for r in range(3)]

    # 3. 元素行 + 计数行(按块展开)
    sym_toks = _tokens(need(5, "元素符号行"))
    if not sym_toks:
        raise ValueError(f"{fname}: 元素符号行不可读")
    for t in sym_toks:
        if len(t) > 4:
            raise ValueError(f"{fname}: 元素符号超 4 字符: {t}")
    cnt_toks = _tokens(need(6, "计数行"))
    if len(cnt_toks) != len(sym_toks):
        raise ValueError(f"{fname}: 计数行字段数与元素行不符: {lines[6]}")
    try:
        cnts = [int(t) for t in cnt_toks]
    except ValueError:
        raise ValueError(f"{fname}: 计数行非整数: {lines[6]}") from None
    if any(c < 0 for c in cnts):
        raise ValueError(f"{fname}: 计数行出现负数")
    nat = sum(cnts)
    if nat <= 0:
        raise ValueError(f"{fname}: 原子数为零")
    symbols: list[str] = []
    for sym, c in zip(sym_toks, cnts):
        symbols.extend([sym] * c)

    # 4. 模式标签(可先有 Selective dynamics 行)
    idx = 7
    tag = _tokens(need(idx, "模式标签行"))
    if tag and tag[0][0] in "Ss":
        idx += 1
        tag = _tokens(need(idx, "模式标签行"))
    if not tag:
        raise ValueError(f"{fname}: 模式标签行不可读")
    head = tag[0][0]
    if head in "Dd":
        direct = True
    elif head in "Cc":
        direct = False
    else:
        raise ValueError(f"{fname}: 未知模式标签(需 Direct 或 Cartesian): {lines[idx]}")
    idx += 1

    # 5. 坐标行(Direct → 行存晶矢线性组合;Cartesian → 乘普适因子)
    coords: list[list[float]] = []
    for i in range(nat):
        toks = _tokens(need(idx + i, f"坐标行 {i + 1}"))
        if len(toks) < 3:
            raise ValueError(f"{fname}: 坐标行 {i + 1} 过短: {lines[idx + i]}")
        try:
            fr = [float(t) for t in toks[:3]]
        except ValueError:
            raise ValueError(f"{fname}: 坐标行 {i + 1} 非数值: {lines[idx + i]}") from None
        if direct:
            coords.append([
                fr[0] * cell[0][c] + fr[1] * cell[1][c] + fr[2] * cell[2][c]
                for c in range(3)
            ])
        else:
            coords.append([s_eff * v for v in fr])

    # 6. 速度块:存在即读跳(行数不足 nat 时 Fortran 侧报软错;GUI 宽容——只警告用)
    name = Path(fname).stem
    frag = Frag(file=fname, name=name, symbols=symbols, coords=coords, cell=cell)
    return frag, cell


def classify_entry(entry: str) -> tuple[str, str]:
    """文件夹条目 → (格式 'xyz'|'poscar', 反应物名);不可分类抛 ValueError。

    sysdef_interface classify_entry 移植:扩展名折叠(.xyz/.poscar)优先,
    裸名 POSCAR* 前缀次之;基名 = 反应物名。
    """
    fname = entry.strip()
    if "." in fname:
        base, ext = fname[: fname.rfind(".")], fname[fname.rfind(".") + 1:].lower()
    else:
        base, ext = fname, ""
    if ext == "xyz":
        fmt = "xyz"
    elif ext == "poscar":
        fmt = "poscar"
    else:
        fmt = "poscar" if len(base) >= 6 and base[:6].lower() == "poscar" else ""
    if not fmt or not base:
        raise ValueError(
            f"未知体系文件类型(需 .xyz / .poscar / POSCAR*): {entry}"
            if base else f"空反应物名(纯扩展名文件名): {entry}"
        )
    return fmt, base


def scan_system_dir(folder: str | Path) -> SysDef:
    """扫描 system/ 文件夹 → SysDef(名序折叠排序/单胞规则/范式/box 派生)。"""
    p = Path(folder)
    out = SysDef(folder=str(p))
    if not p.is_dir():
        out.errors.append(f"体系文件夹不存在: {p}")
        return out
    entries = sorted(q.name for q in p.iterdir() if q.is_file())
    mats: list[tuple[str, str, str]] = []  # (entry, fmt, name)
    for e in entries:
        try:
            fmt, name = classify_entry(e)
        except ValueError as exc:
            out.errors.append(str(exc))
            continue
        mats.append((e, fmt, name))
    if not mats and not out.errors:
        out.errors.append(f"体系文件夹无体系文件: {p}")
        if out.errors and not mats:
            return out
    if not mats:
        return out
    # 名序(大小写折叠)排序 = 片段顺序(sysdef 插入排序的稳定等价)
    mats.sort(key=lambda t: t[2].lower())
    for i in range(1, len(mats)):
        if mats[i][2].lower() == mats[i - 1][2].lower():
            out.errors.append(
                f"反应物名重复(大小写折叠): {mats[i - 1][0]} 与 {mats[i][0]}"
            )
    any_cell = False
    cell_file = ""
    for entry, fmt, name in mats:
        try:
            text = (p / entry).read_text(encoding="utf-8", errors="replace")
            if fmt == "xyz":
                frag = parse_xyz(text, entry)
            else:
                frag, cell = parse_poscar(text, entry)
                if any_cell:
                    out.errors.append(
                        f"多个带胞文件({cell_file} 与 {entry}):单胞是系统级唯一事实"
                    )
                else:
                    any_cell = True
                    cell_file = entry
                    out.cell = cell
            # 质量查表(未知符号即错,镜像 elem_mass 的硬停语义)
            try:
                _ = frag.masses
            except ValueError as exc:
                out.errors.append(f"{entry}: {exc}")
            out.frags.append(frag)
        except ValueError as exc:
            out.errors.append(str(exc))
    out.has_cell = any_cell
    if any_cell and out.cell is not None:
        c = out.cell
        box_lx = math.sqrt(sum(v * v for v in c[0]))
        box_ly = math.sqrt(sum(v * v for v in c[1]))
        if box_lx > 0.0 and box_ly > 0.0:
            out.box_lx = box_lx
            out.box_ly = box_ly
            cos_skew = (c[0][0] * c[1][0] + c[0][1] * c[1][1] + c[0][2] * c[1][2])
            cos_skew = max(-1.0, min(1.0, cos_skew / (box_lx * box_ly)))
            out.skew = math.acos(cos_skew)
    return out
