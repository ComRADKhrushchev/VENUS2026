# -*- coding: utf-8 -*-
"""validate.py - 写盘前静态预检(全部写输入文件时可判的 fatal 规则移植)

设计:
  输入 = 键值字典 + SysDef(system/ 扫描产物)+ 容器键声明表;
  输出 = 错误列表(每条含键名,可定位表单字段)。规则来源:
  input.f90(必填/类型/槽宽/谱规则)、config_atoms.f90(范式窗口/角色)、
  beam_laws.f90 beam_init(模式域/正性/搭配)、inc_surface.f90(B_MAX 拦截/
  A_LAT 公度)、inc_pair.f90(双片段/定点禁)、sampler.f90(方案长度)、
  bath.f90、各 samp_*.f90 成员 init(选中才生效的正性约束)。
"""

from __future__ import annotations

import math

from .keyspec import (
    KEY_TABLE, KEY_NAMES, MEMBER_TABLE, MEMBER_KEY_NAMES, SCHEME_PARAMS,
    parse_scalar, parse_dist_scheme, CONTAINER_KEY_MAX, CONTAINER_KEY_LEN,
)
from .sysdir import SysDef

REQUIRED_KEYS = [s.name for s in KEY_TABLE if s.required]
SLOT_WIDTHS = {s.name: s.slotw for s in KEY_TABLE if s.slotw}


class Issue:
    """一条校验问题(定位到键,便于表单高亮)。"""

    def __init__(self, key: str, msg: str):
        self.key = key
        self.msg = msg

    def __repr__(self):
        return f"[{self.key}] {self.msg}"

    def __str__(self):
        return f"[{self.key}] {self.msg}"


def _num(values: dict, key: str):
    v = values.get(key, "").strip()
    if not v:
        return None
    try:
        return float(v.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def _int(values: dict, key: str):
    v = values.get(key, "").strip()
    if not v:
        return None
    try:
        return int(v)
    except ValueError:
        return None


def _written(values: dict, key: str) -> bool:
    """条件/可选键"是否写出了非空值"。"""
    return bool(values.get(key, "").strip())


def validate(
    values: dict[str, str],
    sysdef: SysDef | None,
    container_declared: list[str] | None = None,
) -> list[Issue]:
    """跑全部静态规则;返回 Issue 列表(空 = 可安全写出)。"""
    issues: list[Issue] = []
    values = {k.strip().upper(): v for k, v in values.items()}

    def err(key: str, msg: str):
        issues.append(Issue(key, msg))

    # ---- 1. 未知键(input.f90 parse:键表 ∪ dispatch ∪ 容器声明)----
    declared = [k.strip().upper() for k in (container_declared or []) if k.strip()]
    known = KEY_NAMES | MEMBER_KEY_NAMES | set(declared)
    for k in values:
        if k not in known:
            err(k, "未知键(不在键表/成员分派表/容器声明内)")

    # ---- 2. 容器键声明域(input.f90 input_declare_keys)----
    if len(declared) > CONTAINER_KEY_MAX:
        err("", f"容器键声明超上限 {CONTAINER_KEY_MAX}")
    if len(set(declared)) != len(declared):
        err("", "容器键声明重复")
    for k in declared:
        if not k:
            continue
        if len(k) > CONTAINER_KEY_LEN:
            err(k, f"容器键超 {CONTAINER_KEY_LEN} 字符")
        if k in KEY_NAMES or k in MEMBER_KEY_NAMES:
            err(k, "容器键与程序键/成员键撞名")

    # ---- 3. 必填五键 + 单值类型/槽宽(input.f90)----
    for r in REQUIRED_KEYS:
        if not values.get(r, "").strip():
            err(r, "取值缺失(该键无系统无关默认)")
    for spec in KEY_TABLE:
        v = values.get(spec.name, "").strip()
        if not v:
            continue
        try:
            parse_scalar(spec.kind, spec.enum, v)
        except ValueError as exc:
            err(spec.name, str(exc))
        if spec.slotw and len(v) > spec.slotw:
            err(spec.name, f"值超槽宽 {spec.slotw} 字符")
    for mk in MEMBER_TABLE:
        v = values.get(mk.key, "").strip()
        if not v:
            continue
        try:
            parse_scalar(mk.kind, "", v)
        except ValueError as exc:
            err(mk.key, str(exc))

    # ---- 4. 谱规则(input.f90:793-798;spectrum_interface)----
    spec_src = (values.get("SPECTRUM_SOURCE", "") or "COMPUTE").strip().upper() or "COMPUTE"
    if _written(values, "SPECTRUM_FILE") and spec_src != "MANUAL":
        err("SPECTRUM_FILE", "SPECTRUM_SOURCE≠MANUAL 时不得写谱表文件")

    # ---- 5. system/ 与散射类相关(config_atoms / sysdef;范式来自 SCATTERING 键,
    #      2026-10-13 裁决:主控键决定,不再由扫描推导)----
    if sysdef is not None and sysdef.ok:
        scat = (values.get("SCATTERING", "") or "GAS").strip().upper() or "GAS"
        surf = scat == "SURFACE"
        names = [f.name.lower() for f in sysdef.frags]
        n_frag = len(sysdef.frags)
        proj = values.get("PROJECTILE", "").strip()
        targ = values.get("TARGET", "").strip()
        if surf:
            if _written(values, "B_MAX"):
                err("B_MAX", "表面散射不支持 B_MAX(横向落点由 N_AIM 落点模式专有,2026-10-01 裁决)")
            if targ:
                err("TARGET", "表面散射不支持 TARGET(仅 PROJECTILE)")
            if not proj:
                err("PROJECTILE", "表面散射须提供 PROJECTILE")
            elif proj.lower() not in names:
                err("PROJECTILE", f"PROJECTILE 未命中片段(可用: {', '.join(sysdef.frag_names)})")
            if not _written(values, "CELL"):
                err("CELL", "表面散射须提供晶胞 CELL(三行晶矢, 如 '6 0 0, 0 6 0, 0 0 20')")
        else:
            if _written(values, "CELL"):
                err("CELL", "气相散射不支持 CELL(表面散射专用键)")
            if n_frag == 1:
                if proj:
                    err("PROJECTILE", "单片段体系不支持 PROJECTILE")
                if targ:
                    err("TARGET", "单片段体系不支持 TARGET")
            elif n_frag > 1:
                if not proj:
                    err("PROJECTILE", "多片段气相散射须提供 PROJECTILE")
                if not targ:
                    err("TARGET", "多片段气相散射须提供 TARGET")
                if proj and proj.lower() not in names:
                    err("PROJECTILE", f"PROJECTILE 未命中片段(可用: {', '.join(sysdef.frag_names)})")
                if targ and targ.lower() not in names:
                    err("TARGET", f"TARGET 未命中片段(可用: {', '.join(sysdef.frag_names)})")
                if proj and targ and proj.lower() == targ.lower():
                    err("TARGET", "PROJECTILE 与 TARGET 须指向不同片段")
        # DIST_SCHEME 长度 = 片段数(sampler.f90:193-198);括号形式
        ds = values.get("DIST_SCHEME", "")
        if ds.strip():
            try:
                items = parse_dist_scheme(ds)
            except ValueError as exc:
                err("DIST_SCHEME", str(exc))
                items = []
            if items and len(items) != n_frag:
                err("DIST_SCHEME", f"DIST_SCHEME 片段数({len(items)})须等于片段数({n_frag})")
            # ebk 载体须双原子(samp_ebk.f90:135)
            for i, (w, _pars) in enumerate(items):
                if w == "EBK" and i < n_frag and sysdef.frags[i].nat != 2:
                    err("N_VIB", f"EBK 载体 {sysdef.frags[i].name} 非双原子(片段 {i + 1})")
    else:
        if sysdef is not None and sysdef.errors:
            err("SYSTEM_DIR", f"system/ 扫描有错: {sysdef.errors[0]}")

    # ---- 5b. 积分器×阶域(dynamics/propagator.f90:126-165/186/223-224:
    #      VERLET 1(Verlet)/2(Beeman)·SYMPLECTIC 4/6/8·RADAU 1-15)----
    integ = (values.get("INTEGRATOR", "") or "VERLET").strip().upper() or "VERLET"
    order_domains = {
        "VERLET": {1, 2}, "SYMPLECTIC": {4, 6, 8}, "RADAU": set(range(1, 16)),
    }
    if integ in order_domains:
        domain = order_domains[integ]
        order = _int(values, "ORDER")
        if order is None:  # 未写出 → 程序默认 1(对 SYMPLECTIC 非法)
            if 1 not in domain:
                err("ORDER", f"{integ} 须写 ORDER ∈ {sorted(domain)}(缺省默认 1 对其非法)")
        elif order not in domain:
            err("ORDER", f"{integ} 的 ORDER 须 ∈ {sorted(domain)}")

    # ---- 6. 束流键(beam_laws beam_init + inc_pair/inc_surface 差异)----
    n_e_rel = _int(values, "N_E_REL")
    if n_e_rel is not None and n_e_rel not in (0, 1, 2):
        err("N_E_REL", "须为 0(热)/1(定值)/2(密度加权)")
    if n_e_rel == 0 and (_num(values, "T_TRANS") or 0.0) <= 0.0:
        err("T_TRANS", "N_E_REL=0 须 T_TRANS > 0")
    if n_e_rel == 2 and (_num(values, "V_WIDTH") or 0.0) <= 0.0:
        err("V_WIDTH", "N_E_REL=2 须 V_WIDTH > 0")
    n_thta = _int(values, "N_THTA")
    thta = _num(values, "THTA_MAX")
    if n_thta is not None and n_thta not in (0, 1):
        err("N_THTA", "须为 0(固定)/1(采样)")
    if thta is not None:
        if n_thta == 1:
            if not (0.0 < thta <= 0.5 * math.pi):
                err("THTA_MAX", "采样极角界须在 (0, π/2]")
        else:
            if not (0.0 <= thta <= 0.5 * math.pi):
                err("THTA_MAX", "固定极角须在 [0, π/2]")
    n_chi = _int(values, "N_CHI")
    if n_chi is not None and n_chi not in (0, 1):
        err("N_CHI", "须为 0(固定)/1(采样)")
    n_b = _int(values, "N_B")
    if n_b is not None and n_b not in (0, 1):
        err("N_B", "须为 0(固定 b)/1(碰撞盘)")
    n_aim = _int(values, "N_AIM")
    if n_aim is not None and n_aim not in (0, 1, 2):
        err("N_AIM", "须为 0(原点)/1(胞内)/2(定点)")
    scat = (values.get("SCATTERING", "") or "GAS").strip().upper() or "GAS"
    if n_aim == 1 and scat != "SURFACE":
        err("N_AIM", "落点模式 1(胞内)须 SCATTERING = SURFACE")
    if n_aim == 2:
        if not (_written(values, "AIM_X") and _written(values, "AIM_Y")):
            err("AIM_X", "N_AIM=2 须同时给 AIM_X 与 AIM_Y")
        if scat != "SURFACE":
            err("N_AIM", "落点模式 2(定点)仅表面散射可用")
    e_rel = _num(values, "E_REL")
    if n_e_rel in (1, 2) and e_rel is not None and e_rel <= 0.0:
        err("E_REL", "碰撞能模式 1/2 须 E_REL > 0")
    b_max = _num(values, "B_MAX")
    r_sep = _num(values, "R_SEP")
    pairish = sysdef is not None and sysdef.ok and scat != "SURFACE" and len(sysdef.frags) > 1
    if pairish:
        if b_max is not None and b_max <= 0.0:
            err("B_MAX", "气相散射碰撞盘半径须 > 0")
        if _written(values, "AIM_X") or _written(values, "AIM_Y"):
            err("AIM_X", "AIM_X/AIM_Y 为表面散射专用键(气相散射不支持)")
        if n_aim == 2:
            err("N_AIM", "pair 通道禁用定点落点(N_AIM=2)")
    if b_max is not None and b_max == 0.0 and n_b == 1:
        err("N_B", "B_MAX=0(定点束)须 N_B=0")
    if r_sep is not None and r_sep <= 0.0:
        err("R_SEP", "初始间距须 > 0")

    # ---- 7. 浴(bath.f90:104-115)----
    bath = (values.get("BATH", "") or "NONE").strip().upper() or "NONE"
    n_bath = _int(values, "N_BATH")
    if n_bath is not None and n_bath < 0:
        err("N_BATH", "须 ≥ 0")
    if n_bath is not None and n_bath > 0:
        if (_num(values, "DT_BATH") or 0.0) <= 0.0:
            err("DT_BATH", "N_BATH>0 须 DT_BATH > 0")
    if bath == "ANDERSEN":
        if (_num(values, "T_BATH") or 0.0) <= 0.0:
            err("T_BATH", "BATH=ANDERSEN 须 T_BATH > 0")
        if (_num(values, "NU_COLL") or 0.0) <= 0.0:
            err("NU_COLL", "BATH=ANDERSEN 须 NU_COLL > 0")

    # ---- 8. 方案词括号参数(2026-10-13 裁定: 逐片段取参;samp_*.f90 init 约束)----
    ds = values.get("DIST_SCHEME", "").strip()
    items: list[tuple[str, dict[str, str]]] = []
    if ds:
        try:
            items = parse_dist_scheme(ds)
        except ValueError as exc:
            err("DIST_SCHEME", str(exc))

    # 各词的必给参数(缺省走成员字段默认;标注"须 >0"默认 0 的参数必须显式给)
    word_required = {
        "BOLTZMANN": {"T_VIB"}, "J": {"T_ROT"},
        "GLO_TARGET": {"T_GLO", "GAMMA_GLO", "W_GHOST"},
        "THERMALIZE": {"T_EQ", "N_EQ", "DT_EQ"},
    }

    def _pf(params: dict, name: str, i: int):
        """片段参数的数值解释(非数值报错,缺失/非数值返回 None)。"""
        raw = params.get(name)
        if raw is None:
            return None
        try:
            return float(raw.replace("D", "E").replace("d", "e"))
        except ValueError:
            err(name, f"片段 {i + 1} 的 {name} 不是数值: '{raw}'")
            return None

    def _pi(params: dict, name: str, i: int):
        raw = params.get(name)
        if raw is None:
            return None
        try:
            return int(raw)
        except ValueError:
            err(name, f"片段 {i + 1} 的 {name} 不是整数: '{raw}'")
            return None

    for i, (word, params) in enumerate(items):
        tag = f"片段 {i + 1}({word})"
        legal = {sp.name for sp in SCHEME_PARAMS.get(word, [])}
        for name in params:
            if name not in legal:
                err(name, f"{tag}: {name} 不是 {word} 的参数(合法: {', '.join(sorted(legal)) or '无'})")
        for name in word_required.get(word, set()):
            if name not in params:
                err(name, f"{tag}: 须给 {name}(其成员默认不合法)")
        if word == "BOLTZMANN":
            v = _pf(params, "T_VIB", i)
            if v is not None and v <= 0.0:
                err("T_VIB", f"{tag}: 须 T_VIB > 0")
        elif word == "J":
            v = _pf(params, "T_ROT", i)
            if v is not None and v <= 0.0:
                err("T_ROT", f"{tag}: 须 T_ROT > 0")
        elif word == "ROTATION":
            v = _pi(params, "J_ROT", i)
            if v is not None and v < 0:
                err("J_ROT", f"{tag}: 须 ≥ 0")
        elif word == "EBK":
            for k in ("N_VIB", "N_ROT"):
                v = _pi(params, k, i)
                if v is not None and v < 0:
                    err(k, f"{tag}: 须 ≥ 0")
        elif word == "SURFACE_OSCILLATOR":
            n_osc = _pi(params, "N_OSC", i)
            if n_osc is not None and n_osc not in (0, 1):
                err("N_OSC", f"{tag}: 须为 0(热)/1(定态)")
            if n_osc == 0:
                v = _pf(params, "T_OSC", i)
                if v is not None and v <= 0.0:
                    err("T_OSC", f"{tag}: N_OSC=0 须 T_OSC > 0(缺省 300 合法)")
            if n_osc == 1:
                v = _pi(params, "N_LEVEL", i)
                if v is not None and v < 0:
                    err("N_LEVEL", f"{tag}: 须 ≥ 0")
            v = _pf(params, "K_OSC", i)
            if v is not None and v <= 0.0:
                err("K_OSC", f"{tag}: 须 > 0")
        elif word == "NORMALMODE":
            v = _pf(params, "E_VIB", i)
            if v is not None and v <= 0.0:
                err("E_VIB", f"{tag}: 须 > 0")
        elif word == "BARRIER_EXCITATION":
            v = _pf(params, "E_STAB", i)
            if v is not None and v < 0.0:
                err("E_STAB", f"{tag}: 须 ≥ 0")
            n_e_bar = _pi(params, "N_E_BAR", i)
            if n_e_bar is not None and n_e_bar not in (0, 1):
                err("N_E_BAR", f"{tag}: 须为 0/1")
            v = _pf(params, "E_BAR", i)
            if v is not None and v < 0.0:
                err("E_BAR", f"{tag}: 须 ≥ 0")
            if n_e_bar == 1:
                v = _pf(params, "T_BAR", i)
                if v is not None and v <= 0.0:
                    err("T_BAR", f"{tag}: N_E_BAR=1 须 T_BAR > 0")
        elif word == "GLO_TARGET":
            for k in ("T_GLO", "GAMMA_GLO", "W_GHOST"):
                v = _pf(params, k, i)
                if v is not None and v <= 0.0:
                    err(k, f"{tag}: 须 > 0")
        elif word == "THERMALIZE":
            v = _pf(params, "T_EQ", i)
            if v is not None and v <= 0.0:
                err("T_EQ", f"{tag}: 须 T_EQ > 0")
            n_eq = _pi(params, "N_EQ", i)
            if n_eq is not None and n_eq < 1:
                err("N_EQ", f"{tag}: 须 ≥ 1")
            v = _pf(params, "DT_EQ", i)
            if v is not None and v <= 0.0:
                err("DT_EQ", f"{tag}: 须 DT_EQ > 0")
        elif word == "WIGNER":
            v = _pf(params, "E_SCALE", i)
            if v is not None and v <= 0.0:
                err("E_SCALE", f"{tag}: 须 > 0")

    # THERMALIZE 共享一次平衡跑: N_EQ/DT_EQ 各片段须一致(samp_thermalize init)
    th = [(i, p) for i, (w, p) in enumerate(items) if w == "THERMALIZE"]
    for k in ("N_EQ", "DT_EQ"):
        vals = {p.get(k) for _, p in th if k in p}
        if len(vals) > 1:
            err(k, f"THERMALIZE 各片段的 {k} 须一致(共享一次平衡跑): {sorted(vals)}")

    # ---- 9. MANUAL 谱表(bath 同级的显式检查)----
    if spec_src == "MANUAL" and not _written(values, "SPECTRUM_FILE"):
        err("SPECTRUM_FILE", "SPECTRUM_SOURCE=MANUAL 须给谱表文件")

    return issues
