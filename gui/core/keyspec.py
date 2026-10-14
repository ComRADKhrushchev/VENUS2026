# -*- coding: utf-8 -*-
"""keyspec.py - 输入键表镜像(src/input.f90 key_tab/dispatch_tbl 的 Python 侧单一副本)

设计:
  与 Fortran 键表逐行对应(名称/块/类型/枚举族/必填/默认/槽宽/条件键),由
  tests/test_keyspec.py 的镜像同步测试对 src/input.f90 源文本比对防漂移。
  通道级成员键(dispatch_tbl 13 行: 入射束流定律与浴物种)附类型与约束摘要;
  方案词自身的参数只经 DIST_SCHEME 括号(2026-10-13 裁定: 词[名=值, ...],
  无顶层成员键、无回退)—— SCHEME_PARAMS 为逐词括号词表,约束语义移植自
  各 samp_*.f90 的 init 校验(写出前静态可判部分)。
  本模块只有数据与单值解析,跨键约束在 validate.py。
"""

from dataclasses import dataclass, field

# 值类型(与 input.f90 kind_* 对应)
KIND_INT = "int"          # 整数标量(禁逗号列表)
KIND_REAL = "real"        # 实数标量
KIND_CHAR = "char"        # 字符串(原样保留)
KIND_LOG = "log"          # 逻辑 T/F/TRUE/FALSE/.TRUE./.FALSE./1/0
KIND_LEVEL = "level"      # 整数限于 {0,1,2,5}
KIND_ENUM = "enum"        # 枚举词(标量)
KIND_ELIST = "elist"      # 逗号分隔枚举词数组
KIND_CELL = "cell"        # 3x3 晶胞:三逗号组各三实数(CELL 专用)

BLK_FLOW = "flow"         # 块一:运行控制
BLK_SYSTEM = "system"     # 块二:体系定义

LOGIC_WORDS = {"T", "F", "TRUE", "FALSE", ".TRUE.", ".FALSE.", "1", "0"}
LEVEL_SET = {0, 1, 2, 5}

# 枚举注册表(input.f90 enum_*_tbl)
ENUM_TASK = ["TRAJECTORY"]
ENUM_INTEG = ["VERLET", "SYMPLECTIC", "RADAU"]
ENUM_DIST = [
    "BOLTZMANN", "EBK", "ROTATION", "J", "SURFACE_OSCILLATOR",
    "NORMALMODE", "BARRIER_EXCITATION", "GLO_TARGET", "THERMALIZE", "WIGNER",
]
ENUM_BATH = ["NONE", "ANDERSEN"]
ENUM_SPECTRUM = ["COMPUTE", "MANUAL"]
ENUM_SCATTERING = ["GAS", "SURFACE"]
ENUMS = {
    "task": ENUM_TASK,
    "integ": ENUM_INTEG,
    "dist": ENUM_DIST,
    "bath": ENUM_BATH,
    "spectrum": ENUM_SPECTRUM,
    "scattering": ENUM_SCATTERING,
}


@dataclass(frozen=True)
class KeySpec:
    """一行程序级键(镜像 key_spec_t)。"""
    name: str
    blk: str
    kind: str
    enum: str = ""            # ENUMS 的键名(enum_*_tbl)
    required: bool = False
    default: str = ""
    slotw: int = 0            # 槽宽(0 = 无固定槽)
    cond: bool = False        # 范式条件键(永不注入默认、永不无条件必填)


# 程序级键表(顺序 = key_tab 顺序 = 写出顺序;默认值逐字取自 input.f90:106-138)
KEY_TABLE: list[KeySpec] = [
    KeySpec("TASK", BLK_FLOW, KIND_ENUM, "task", False, "TRAJECTORY"),
    KeySpec("N_TRAJ", BLK_FLOW, KIND_INT, "", False, "1"),
    KeySpec("TITLE", BLK_FLOW, KIND_CHAR, "", False, "", 80),
    KeySpec("ISEED", BLK_FLOW, KIND_INT, "", False, "1"),
    KeySpec("MAX_STEPS", BLK_FLOW, KIND_INT, "", False, "100000"),
    KeySpec("INTEGRATOR", BLK_FLOW, KIND_ENUM, "integ", False, "VERLET"),
    KeySpec("DT", BLK_FLOW, KIND_REAL, "", True, ""),
    KeySpec("ORDER", BLK_FLOW, KIND_INT, "", False, "1"),
    KeySpec("DT_MIN", BLK_FLOW, KIND_REAL, "", False, "0.0"),
    KeySpec("DT_MAX", BLK_FLOW, KIND_REAL, "", False, "0.0"),
    KeySpec("ELEC_METHOD", BLK_FLOW, KIND_CHAR, "", False, "ADIABATIC", 32),
    KeySpec("GWRITE_LEVEL", BLK_FLOW, KIND_LEVEL, "", False, "2"),
    KeySpec("OBSERVABLES", BLK_FLOW, KIND_CHAR, "", False, "", 256),
    KeySpec("SYSTEM_DIR", BLK_SYSTEM, KIND_CHAR, "", True, "", 128),
    KeySpec("E_REL", BLK_SYSTEM, KIND_REAL, "", True, ""),
    KeySpec("B_MAX", BLK_SYSTEM, KIND_REAL, "", False, "0.0"),
    KeySpec("R_SEP", BLK_SYSTEM, KIND_REAL, "", True, ""),
    KeySpec("DIST_SCHEME", BLK_SYSTEM, KIND_ELIST, "dist", True, ""),
    KeySpec("RANDOM_ORIENT", BLK_SYSTEM, KIND_LOG, "", False, "F"),
    KeySpec("BATH", BLK_SYSTEM, KIND_ENUM, "bath", False, "NONE"),
    KeySpec("N_BATH", BLK_SYSTEM, KIND_INT, "", False, "0"),
    KeySpec("DT_BATH", BLK_SYSTEM, KIND_REAL, "", False, "0.0"),
    KeySpec("SCATTERING", BLK_SYSTEM, KIND_ENUM, "scattering", False, "GAS"),
    KeySpec("CELL", BLK_SYSTEM, KIND_CELL, "", False, "", cond=True),
    KeySpec("SPECTRUM_SOURCE", BLK_SYSTEM, KIND_ENUM, "spectrum", False, "COMPUTE"),
    KeySpec("SPECTRUM_FILE", BLK_SYSTEM, KIND_CHAR, "", False, "", 128),
    KeySpec("PROJECTILE", BLK_SYSTEM, KIND_CHAR, "", False, "", 64, cond=True),
    KeySpec("TARGET", BLK_SYSTEM, KIND_CHAR, "", False, "", 64, cond=True),
]

KEY_NAMES = {k.name for k in KEY_TABLE}


@dataclass(frozen=True)
class MemberKey:
    """一行成员参数键(镜像 dispatch_tbl,附 GUI 侧类型与约束提示)。"""
    key: str
    member: str          # 归属成员(dispatch 表原文,如 'incident_surface/pair')
    kind: str            # GUI 侧解析类型(int/real)
    hint: str = ""       # 中文说明


# 通道级成员参数键(input.f90 dispatch_tbl 13 行顺序;类型与约束摘自各成员 init)
MEMBER_TABLE: list[MemberKey] = [
    MemberKey("N_E_REL", "incident_surface/pair", KIND_INT, "碰撞能模式 0=热分布/1=定值/2=密度加权"),
    MemberKey("T_TRANS", "incident_surface/pair", KIND_REAL, "平动温度 [K](N_E_REL=0 时须 >0)"),
    MemberKey("V_WIDTH", "incident_surface/pair", KIND_REAL, "速度展宽 [Å/fs](N_E_REL=2 时须 >0)"),
    MemberKey("N_THTA", "incident_surface/pair", KIND_INT, "极角模式 0=固定/1=按 sinθ 采样"),
    MemberKey("THTA_MAX", "incident_surface/pair", KIND_REAL, "极角界/固定值 [rad](0..π/2)"),
    MemberKey("N_CHI", "incident_surface/pair", KIND_INT, "方位角模式 0=固定/1=均匀采样"),
    MemberKey("CHI", "incident_surface/pair", KIND_REAL, "固定方位角 [rad](N_CHI=0)"),
    MemberKey("N_B", "incident_surface/pair", KIND_INT, "瞄准模式 0=固定 b=B_MAX/1=碰撞盘采样"),
    MemberKey("N_AIM", "incident_surface/pair", KIND_INT, "落点模式 0=原点/1=胞内均匀/2=定点(仅表面)"),
    MemberKey("AIM_X", "incident_surface", KIND_REAL, "定点落点 x [Å, CELL 帧](N_AIM=2,仅表面)"),
    MemberKey("AIM_Y", "incident_surface", KIND_REAL, "定点落点 y [Å, CELL 帧](N_AIM=2,仅表面)"),
    MemberKey("T_BATH", "bath_andersen", KIND_REAL, "Andersen 浴温度 [K](BATH=ANDERSEN 时须 >0)"),
    MemberKey("NU_COLL", "bath_andersen", KIND_REAL, "Andersen 碰撞频率(同上,须 >0)"),
]

MEMBER_KEY_NAMES = {m.key for m in MEMBER_TABLE}

@dataclass(frozen=True)
class SchemeParam:
    """一个方案词的括号参数(2026-10-13 裁定: 只经 DIST_SCHEME 括号取参)。"""
    name: str
    kind: str            # int/real
    hint: str = ""


# 方案词 → 括号参数词表(空表 = 该词无参数,写成裸词);required=False 的参数
# 缺省时走成员字段默认(与 Fortran 成员默认一致)
SCHEME_PARAMS: dict[str, list[SchemeParam]] = {
    "BOLTZMANN": [SchemeParam("T_VIB", KIND_REAL, "振动温度 [K](须 >0)")],
    "EBK": [SchemeParam("N_VIB", KIND_INT, "振动量子数(≥0,缺省 0)"),
            SchemeParam("N_ROT", KIND_INT, "转动量子数(≥0,缺省 0)")],
    "ROTATION": [SchemeParam("J_ROT", KIND_INT, "固定转动角动量量子数 J(≥0,缺省 0)")],
    "J": [SchemeParam("T_ROT", KIND_REAL, "转动温度 [K](须 >0)")],
    "SURFACE_OSCILLATOR": [SchemeParam("N_OSC", KIND_INT, "能量律 0=热/1=定态(缺省 0)"),
                           SchemeParam("T_OSC", KIND_REAL, "振子温度 [K](缺省 300)"),
                           SchemeParam("N_LEVEL", KIND_INT, "定态能级(≥0,缺省 0)"),
                           SchemeParam("K_OSC", KIND_REAL, "力常数 [eV/Å²](缺省 1.30092582)")],
    "NORMALMODE": [SchemeParam("E_VIB", KIND_REAL, "振动能 [eV](须 >0,缺省 0.43364194)")],
    "BARRIER_EXCITATION": [SchemeParam("E_STAB", KIND_REAL, "鞍点稳定模分能 [eV](≥0)"),
                           SchemeParam("N_E_BAR", KIND_INT, "反应坐标能量律 0=定值/1=热"),
                           SchemeParam("E_BAR", KIND_REAL, "反应坐标能量 [eV](≥0)"),
                           SchemeParam("T_BAR", KIND_REAL, "反应坐标温度 [K](N_E_BAR=1 须 >0)")],
    "GLO_TARGET": [SchemeParam("T_GLO", KIND_REAL, "GLO 靶温 [K](须 >0)"),
                   SchemeParam("GAMMA_GLO", KIND_REAL, "GLO 摩擦系数(须 >0)"),
                   SchemeParam("W_GHOST", KIND_REAL, "ghost 频率 [1/fs](须 >0)")],
    "THERMALIZE": [SchemeParam("T_EQ", KIND_REAL, "平衡温度 [K](须 >0)"),
                   SchemeParam("N_EQ", KIND_INT, "平衡步数(≥1;各片段须一致)"),
                   SchemeParam("DT_EQ", KIND_REAL, "平衡步长 [fs](须 >0;各片段须一致)")],
    "WIGNER": [SchemeParam("E_SCALE", KIND_REAL, "能量标度(零点能份额,须 >0,缺省 1)")],
}


def _top_split(raw: str) -> list[str]:
    """按括号深度 0 的逗号切分(括号内逗号是参数分隔符)。"""
    items, depth, cur = [], 0, []
    for ch in raw:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == "," and depth == 0:
            items.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    items.append("".join(cur).strip())
    return [it for it in items if it != ""]


def parse_dist_scheme(raw: str) -> list[tuple[str, dict[str, str]]]:
    """解析 DIST_SCHEME 值: 词或 词[名=值, ...](逐片段)。

    返回 [(词, {参数名: 原文值}), ...];形式错误抛 ValueError(中文消息,
    与 src/input.f90 parse_elist 的命名错误一一对应)。大小写不敏感,词与
    参数名折为大写。
    """
    if raw.count("[") != raw.count("]"):
        raise ValueError("DIST_SCHEME: 括号不配对(每个 [ 须在片段项内闭合)")
    out: list[tuple[str, dict[str, str]]] = []
    for item in _top_split(raw):
        if item == "":
            raise ValueError("DIST_SCHEME: 空片段项")
        br = item.find("[")
        params: dict[str, str] = {}
        if br == 0:
            raise ValueError(f'DIST_SCHEME: 片段项 "{item}" 括号前无方案词')
        if br > 0:
            word = item[:br].strip().upper()
            tail = item[br + 1:].strip()
            if not tail.endswith("]"):
                raise ValueError(f'DIST_SCHEME: 片段项 "{item}" 闭合括号后有残留文本')
            body = tail[:-1].strip()
            if body == "":
                raise ValueError(f'DIST_SCHEME: 空括号 "{item}"(无参数的词写裸词)')
            if "[" in body or "]" in body:
                raise ValueError(f'DIST_SCHEME: 片段项 "{item}" 括号嵌套(参数为平铺 名=值)')
            seen: set[str] = set()
            for pair in body.split(","):
                pair = pair.strip()
                if pair == "":
                    raise ValueError(f'DIST_SCHEME: 片段项 "{item}" 有空参数项')
                if "=" not in pair:
                    raise ValueError(f'DIST_SCHEME: 参数 "{pair}" 不是 名=值 形式')
                name, val = pair.split("=", 1)
                name, val = name.strip().upper(), val.strip()
                if name == "":
                    raise ValueError(f'DIST_SCHEME: 参数 "{pair}" 无参数名')
                if val == "":
                    raise ValueError(f'DIST_SCHEME: 参数 {name} 取值缺失(该键无系统无关默认)')
                if not all(c.isalnum() or c == "_" for c in name):
                    raise ValueError(f'DIST_SCHEME: 参数名 "{name}" 含非法字符(字母/数字/下划线)')
                if name in seen:
                    raise ValueError(f'DIST_SCHEME: 参数 {name} 在同一片段内重复(一载体一真相源)')
                seen.add(name)
                params[name] = val
        else:
            word = item.upper()
        if word not in ENUM_DIST:
            raise ValueError(f'DIST_SCHEME: 非法方案词 "{word}"(合法: {", ".join(ENUM_DIST)})')
        out.append((word, params))
    if not out:
        raise ValueError("DIST_SCHEME: 空值列表")
    return out


def format_dist_scheme(items: list[tuple[str, dict[str, str]]]) -> str:
    """构造 DIST_SCHEME 值: 逐片段 词[名=值, ...](无参数片段写裸词)。"""
    parts = []
    for word, params in items:
        if params:
            inner = ", ".join(f"{k}={v}" for k, v in params.items())
            parts.append(f"{word}[{inner}]")
        else:
            parts.append(word)
    return ", ".join(parts)


# 束流键(入射通道参数;成员名 dispatch 表记 'incident_surface/pair')
BEAM_KEYS = [
    "N_E_REL", "T_TRANS", "V_WIDTH", "N_THTA", "THTA_MAX", "N_CHI", "CHI",
    "N_B", "N_AIM", "AIM_X", "AIM_Y",
]

# 容器键通道上限(input.f90:69-71, 278-306)
CONTAINER_KEY_MAX = 128
CONTAINER_KEY_LEN = 32

# reference 容器的键词表(reference_box.f90 keys();GUI 预填表用,非契约)
REFERENCE_CONTAINER_KEYS = ["K_HARM", "Q_TERM", "E_SPLIT"]


def parse_scalar(kind: str, enum: str, raw: str):
    """按 kind 解析单值;失败抛 ValueError(消息中文,面向表单)。

    与 input.f90 的类型域检查对应:整数禁逗号列表;level 限 {0,1,2,5};
    逻辑词表;枚举词大小写不敏感;posreal 解析为实数且 >0。
    """
    s = raw.strip()
    if kind == KIND_INT:
        if "," in s:
            raise ValueError(f"应为单个整数,收到列表: {raw}")
        try:
            return int(s)
        except ValueError:
            raise ValueError(f"应为整数: {raw}") from None
    if kind == KIND_REAL:
        try:
            return float(s.replace("D", "E").replace("d", "e"))
        except ValueError:
            raise ValueError(f"应为实数: {raw}") from None
    if kind == KIND_CELL:
        groups = [g.strip() for g in s.split(",")]
        if len(groups) != 3 or not all(groups):
            raise ValueError(
                f"应为三行晶胞矢量(三逗号组各三实数, 如 '6 0 0, 0 6 0, 0 0 20'): {raw}")
        cell = []
        for i, g in enumerate(groups, 1):
            toks = g.split()
            if len(toks) != 3:
                raise ValueError(f"晶胞第 {i} 行应为三个数, 收到 {len(toks)} 个: '{g}'")
            try:
                cell.append([float(t.replace("D", "E").replace("d", "e")) for t in toks])
            except ValueError:
                raise ValueError(f"晶胞第 {i} 行含非数值: '{g}'") from None
        return cell
    if kind == KIND_LEVEL:
        v = parse_scalar(KIND_INT, "", raw)
        if v not in LEVEL_SET:
            raise ValueError(f"GWRITE_LEVEL 限 {{0,1,2,5}}: {raw}")
        return v
    if kind == KIND_LOG:
        u = s.upper()
        if u not in LOGIC_WORDS:
            raise ValueError(f"应为逻辑词 T/F/TRUE/FALSE/.TRUE./.FALSE./1/0: {raw}")
        return u in ("T", "TRUE", ".TRUE.", "1")
    if kind == KIND_ENUM:
        u = s.upper()
        words = ENUMS.get(enum, [])
        if u not in words:
            raise ValueError(f"应为注册词 {'/'.join(words)}: {raw}")
        return u
    if kind == KIND_ELIST:
        # 括号形式(2026-10-13 裁定): 逐片段 词[名=值, ...],委托专用解析
        items = parse_dist_scheme(s)
        return [w for w, _p in items]
    return s  # char 原样
