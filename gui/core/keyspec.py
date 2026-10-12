# -*- coding: utf-8 -*-
"""keyspec.py - 输入键表镜像(src/input.f90 key_tab/dispatch_tbl 的 Python 侧单一副本)

设计:
  与 Fortran 键表逐行对应(名称/块/类型/枚举族/必填/默认/槽宽/条件键),由
  tests/test_keyspec.py 的镜像同步测试对 src/input.f90 源文本比对防漂移。
  成员参数键(dispatch_tbl 35 行)附类型与"选中才生效"的约束摘要,约束
  语义移植自各 samp_*.f90 / beam_laws.f90 的 init 校验(写出前静态可判部分)。
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
KIND_POSREAL = "posreal"  # 正实数(A_LAT 专用)

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
ENUMS = {
    "task": ENUM_TASK,
    "integ": ENUM_INTEG,
    "dist": ENUM_DIST,
    "bath": ENUM_BATH,
    "spectrum": ENUM_SPECTRUM,
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
    KeySpec("N_SURF", BLK_FLOW, KIND_INT, "", False, "0"),
    KeySpec("PES_MAIN", BLK_FLOW, KIND_CHAR, "", False, "1", 32),
    KeySpec("GWRITE_LEVEL", BLK_FLOW, KIND_LEVEL, "", False, "2"),
    KeySpec("OBSERVABLES", BLK_FLOW, KIND_CHAR, "", False, "", 256),
    KeySpec("E_KIN", BLK_FLOW, KIND_LOG, "", False, "F"),
    KeySpec("E_POT", BLK_FLOW, KIND_LOG, "", False, "F"),
    KeySpec("E_TOT", BLK_FLOW, KIND_LOG, "", False, "F"),
    KeySpec("SYSTEM_DIR", BLK_SYSTEM, KIND_CHAR, "", True, "", 128),
    KeySpec("E_REL", BLK_SYSTEM, KIND_REAL, "", True, ""),
    KeySpec("B_MAX", BLK_SYSTEM, KIND_REAL, "", False, "0.0"),
    KeySpec("R_SEP", BLK_SYSTEM, KIND_REAL, "", True, ""),
    KeySpec("DIST_SCHEME", BLK_SYSTEM, KIND_ELIST, "dist", True, ""),
    KeySpec("RANDOM_ORIENT", BLK_SYSTEM, KIND_LOG, "", False, "F"),
    KeySpec("BATH", BLK_SYSTEM, KIND_ENUM, "bath", False, "NONE"),
    KeySpec("N_BATH", BLK_SYSTEM, KIND_INT, "", False, "0"),
    KeySpec("DT_BATH", BLK_SYSTEM, KIND_REAL, "", False, "0.0"),
    KeySpec("SPECTRUM_SOURCE", BLK_SYSTEM, KIND_ENUM, "spectrum", False, "COMPUTE"),
    KeySpec("SPECTRUM_FILE", BLK_SYSTEM, KIND_CHAR, "", False, "", 128),
    KeySpec("PROJECTILE", BLK_SYSTEM, KIND_CHAR, "", False, "", 64, cond=True),
    KeySpec("TARGET", BLK_SYSTEM, KIND_CHAR, "", False, "", 64, cond=True),
    KeySpec("A_LAT", BLK_SYSTEM, KIND_POSREAL, "", False, "", cond=True),
]

KEY_NAMES = {k.name for k in KEY_TABLE}


@dataclass(frozen=True)
class MemberKey:
    """一行成员参数键(镜像 dispatch_tbl,附 GUI 侧类型与约束提示)。"""
    key: str
    member: str          # 归属成员(dispatch 表原文,如 'incident_surface/pair')
    kind: str            # GUI 侧解析类型(int/real)
    hint: str = ""       # 中文说明


# 成员参数键(input.f90:183-218 dispatch_tbl 顺序;类型与约束摘自各成员 init)
MEMBER_TABLE: list[MemberKey] = [
    MemberKey("T_VIB_A", "boltzmann", KIND_REAL, "Boltzmann 振动温度 [K](选中 boltzmann 时须 >0)"),
    MemberKey("T_VIB_B", "boltzmann", KIND_REAL, "Boltzmann 第二通道温度 [K](保留未用)"),
    MemberKey("N_E_REL", "incident_surface/pair", KIND_INT, "碰撞能模式 0=热分布/1=定值/2=密度加权"),
    MemberKey("T_TRANS", "incident_surface/pair", KIND_REAL, "平动温度 [K](N_E_REL=0 时须 >0)"),
    MemberKey("V_WIDTH", "incident_surface/pair", KIND_REAL, "速度展宽 [Å/(10 fs)](N_E_REL=2 时须 >0)"),
    MemberKey("N_THTA", "incident_surface/pair", KIND_INT, "极角模式 0=固定/1=按 sinθ 采样"),
    MemberKey("THTA_MAX", "incident_surface/pair", KIND_REAL, "极角界/固定值 [rad](0..π/2)"),
    MemberKey("N_CHI", "incident_surface/pair", KIND_INT, "方位角模式 0=固定/1=均匀采样"),
    MemberKey("CHI", "incident_surface/pair", KIND_REAL, "固定方位角 [rad](N_CHI=0)"),
    MemberKey("N_B", "incident_surface/pair", KIND_INT, "瞄准模式 0=固定 b=B_MAX/1=碰撞盘采样"),
    MemberKey("N_AIM", "incident_surface/pair", KIND_INT, "落点模式 0=原点/1=胞内均匀/2=定点(仅表面)"),
    MemberKey("AIM_X", "incident_surface", KIND_REAL, "定点落点 x [Å, POSCAR 帧](N_AIM=2,仅表面)"),
    MemberKey("AIM_Y", "incident_surface", KIND_REAL, "定点落点 y [Å, POSCAR 帧](N_AIM=2,仅表面)"),
    MemberKey("T_ROT", "j", KIND_REAL, "转动温度 [K](选中 j 时须 >0)"),
    MemberKey("J_ROT", "rotation", KIND_INT, "固定转动角动量量子数 J(≥0)"),
    MemberKey("N_OSC", "surface_oscillator", KIND_INT, "表面振子模式 0=热/1=定态"),
    MemberKey("T_OSC", "surface_oscillator", KIND_REAL, "振子温度 [K](N_OSC=0 时须 >0)"),
    MemberKey("N_LEVEL", "surface_oscillator", KIND_INT, "定态能级(≥0,N_OSC=1)"),
    MemberKey("K_OSC", "surface_oscillator", KIND_REAL, "振子力常数 [eV/Å²](须 >0)"),
    MemberKey("E_VIB", "normalmode", KIND_REAL, "振动能 [kcal/mol](须 >0)"),
    MemberKey("E_STAB", "barrier_excitation", KIND_REAL, "鞍点稳定模分能 [kcal/mol](≥0)"),
    MemberKey("N_E_BAR", "barrier_excitation", KIND_INT, "反应坐标能量模式 0=定值/1=热"),
    MemberKey("E_BAR", "barrier_excitation", KIND_REAL, "反应坐标能量 [kcal/mol](≥0)"),
    MemberKey("T_BAR", "barrier_excitation", KIND_REAL, "反应坐标温度 [K](N_E_BAR=1 时须 >0)"),
    MemberKey("T_BATH", "bath_andersen", KIND_REAL, "Andersen 浴温度 [K](BATH=ANDERSEN 时须 >0)"),
    MemberKey("NU_COLL", "bath_andersen", KIND_REAL, "Andersen 碰撞频率(同上,须 >0)"),
    MemberKey("T_GLO", "glo_target", KIND_REAL, "GLO 靶温 [K](须 >0)"),
    MemberKey("GAMMA_GLO", "glo_target", KIND_REAL, "GLO 摩擦系数(须 >0)"),
    MemberKey("W_GHOST", "glo_target", KIND_REAL, "ghost 质量 [amu](须 >0)"),
    MemberKey("T_EQ", "thermalize", KIND_REAL, "平衡温度 [K](须 >0)"),
    MemberKey("N_EQ", "thermalize", KIND_INT, "平衡步数(≥1)"),
    MemberKey("DT_EQ", "thermalize", KIND_REAL, "平衡步长 [10 fs](须 >0)"),
    MemberKey("E_SCALE", "wigner", KIND_REAL, "Wigner 能量标度 [kcal/mol](须 >0)"),
    MemberKey("N_VIB", "ebk", KIND_INT, "EBK 振动量子数(≥0,载体须双原子)"),
    MemberKey("N_ROT", "ebk", KIND_INT, "EBK 转动量子数(≥0,载体须双原子)"),
]

MEMBER_KEY_NAMES = {m.key for m in MEMBER_TABLE}

# DIST_SCHEME 词 → 该词激活的成员参数键(动态显隐用)
SCHEME_MEMBER_KEYS: dict[str, list[str]] = {
    "BOLTZMANN": ["T_VIB_A", "T_VIB_B"],
    "EBK": ["N_VIB", "N_ROT"],
    "ROTATION": ["J_ROT"],
    "J": ["T_ROT"],
    "SURFACE_OSCILLATOR": ["N_OSC", "T_OSC", "N_LEVEL", "K_OSC"],
    "NORMALMODE": ["E_VIB"],
    "BARRIER_EXCITATION": ["E_STAB", "N_E_BAR", "E_BAR", "T_BAR"],
    "GLO_TARGET": ["T_GLO", "GAMMA_GLO", "W_GHOST"],
    "THERMALIZE": ["T_EQ", "N_EQ", "DT_EQ"],
    "WIGNER": ["E_SCALE"],
}

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
    if kind in (KIND_REAL, KIND_POSREAL):
        try:
            v = float(s.replace("D", "E").replace("d", "e"))
        except ValueError:
            raise ValueError(f"应为实数: {raw}") from None
        if kind == KIND_POSREAL and v <= 0.0:
            raise ValueError(f"应为正实数: {raw}")
        return v
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
        out = []
        for tok in s.split(","):
            tok = tok.strip()
            if not tok:
                continue
            u = tok.upper()
            if u not in ENUMS.get(enum, []):
                raise ValueError(f"DIST_SCHEME 含未注册词: {tok}")
            out.append(u)
        return out
    return s  # char 原样
