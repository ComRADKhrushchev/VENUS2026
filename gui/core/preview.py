# -*- coding: utf-8 -*-
"""preview.py - 起始散射构型装配(beam_laws.f90 / inc_surface.f90 / inc_pair.f90 公式移植)

设计:
  静态 3D 预览的计算核:片段内部几何取体系文件平衡坐标(不抽样振动/
  转动),COM 放置/速度赋值逐式移植 Fortran 侧——
    u_inc  = (−sinθ·cosχ, −sinθ·sinχ, −cosθ)                 (inc_surface.f90:126)
    p_site = (x0 + b·cosφ + R_SEP·sinθ·cosχ,
              y0 + b·sinφ + R_SEP·sinθ·sinχ, R_SEP)           (:127-129)
    pair 速度劈分 w1 = m2/(m1+m2)、w2 = m1/(m1+m2)            (inc_pair.f90:115-118)
    v = √(2·E_REL/(μ·dyn_ev)), dyn_ev = 103.6426965           (beam_laws.f90, eV/fs 制)
  抽样量(N_B/N_AIM/N_THTA/N_CHI 采样时)以**代表值**代入并在 notes 里
  标注,面板须显示"示意值"字样;入射通道选择 = samp_incident_word 规则
  (表面→surface;分子双片段→pair;单片段→single 不放置)。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .keyspec import KIND_CELL, parse_scalar
from .sysdir import SysDef, Frag

DYN_EV = 103.6426965  # 1 amu·Å²/fs² 折 eV(consts.f90 dyn_ev,镜像不得重推)


def incident_word(surface_model: int, n_frag: int) -> str:
    """入射通道选择(sampler.f90:132-146 规则;surface_model 来自 SCATTERING 键)。"""
    if surface_model != 0:
        return "incident_surface"
    return "incident_pair" if n_frag > 1 else "incident_single"


def cell_from_values(values: dict[str, str]):
    """解析 CELL 键(三行晶矢);缺席或坏值返回 None(2026-10-13 裁决:胞来自输入)。"""
    raw = values.get("CELL", "").strip()
    if not raw:
        return None
    try:
        cell = parse_scalar(KIND_CELL, "", raw)
    except ValueError:
        return None
    return cell


def box_from_cell(cell) -> tuple[float, float, float]:
    """面内盒标量(box_lx/box_ly/skew)——余弦定理,镜像 config_atoms 装配步 4。"""
    a1, a2 = cell[0], cell[1]
    lx = math.sqrt(sum(c * c for c in a1))
    ly = math.sqrt(sum(c * c for c in a2))
    cos_skew = max(-1.0, min(1.0, sum(x * y for x, y in zip(a1, a2)) / (lx * ly)))
    return lx, ly, math.acos(cos_skew)


def _com(frag: Frag) -> list[float]:
    m = frag.masses
    tot = sum(m)
    return [sum(m[i] * frag.coords[i][c] for i in range(frag.nat)) / tot for c in range(3)]


def _shift_com(frag: Frag, com_new: list[float]) -> list[list[float]]:
    """把片段原子平移到指定 COM(beam_com_shift 等价;内部几何不变)。"""
    com = _com(frag)
    return [
        [frag.coords[i][c] + com_new[c] - com[c] for c in range(3)]
        for i in range(frag.nat)
    ]


@dataclass
class BeamView:
    """装配所需的束流代表值(全部显式,便于界面复述)。"""
    theta: float          # 极角 [rad]
    chi: float            # 方位角 [rad]
    b: float              # 瞄准半径 [Å]
    phi: float            # 盘方位角 [rad]
    drop: list[float]     # 落点 (x0, y0) [Å]
    speed: float          # 碰撞速度 [Å/fs]
    u_inc: list[float]    # 入射单位矢量
    p_site: list[float]   # 弹丸 COM 落位 [Å]
    mu: float             # 有效质量(表面=弹丸质量;pair=折合质量) [amu]
    v_proj: list[float]   # 弹丸 COM 速度 [Å/fs]
    v_targ: list[float] | None = None  # 靶 COM 速度(pair 非零)


@dataclass
class PreviewScene:
    """3D 预览场景:全原子坐标 + 注记元素。"""
    symbols: list[str] = field(default_factory=list)
    masses: list[float] = field(default_factory=list)
    positions: list[list[float]] = field(default_factory=list)  # (n,3) [Å]
    frag_slices: list[tuple[int, int]] = field(default_factory=list)  # 每片段 (起,止)
    frag_names: list[str] = field(default_factory=list)
    frag_coms: list[list[float]] = field(default_factory=list)
    frag_vels: list[list[float]] = field(default_factory=list)   # 每片段 COM 速度
    cell: list[list[float]] | None = None     # 3x3 晶矢(表面)
    beam: BeamView | None = None
    b_max: float | None = None                # 碰撞盘半径(pair)
    aim_marker: list[float] | None = None     # 瞄准点 (x0+b·cosφ, y0+b·sinφ, 0)
    notes: list[str] = field(default_factory=list)  # "示意值"注记


def assemble_preview(
    sysdef: SysDef,
    values: dict[str, str],
    i_proj: int = -1,
    i_targ: int = -1,
) -> PreviewScene:
    """从体系定义与输入键值装配初始构型场景。

    values: 输入键 KEY→值字符串(与 inputwrite 同一套);
    i_proj/i_targ: PROJECTILE/TARGET 命中的片段行(0 基;-1 = 未选)。
    """
    def num(key: str, default: float = 0.0) -> float:
        v = values.get(key, "").strip()
        if not v:
            return default
        try:
            return float(v.replace("D", "E").replace("d", "e"))
        except ValueError:
            return default

    def integer(key: str, default: int = 0) -> int:
        v = values.get(key, "").strip()
        if not v:
            return default
        try:
            return int(v)
        except ValueError:
            return default

    scat = (values.get("SCATTERING", "") or "GAS").strip().upper() or "GAS"
    surf = scat == "SURFACE"
    cell = cell_from_values(values) if surf else None
    scene = PreviewScene(cell=cell)
    scene.notes.append("预览为示意:片段内部几何取文件平衡坐标(未抽样振动/转动)")

    n_frag = len(sysdef.frags)
    channel = incident_word(1 if surf else 0, n_frag)
    e_rel = num("E_REL")
    r_sep = num("R_SEP")

    # 各片段先落到自身文件坐标(表面衬底/单片段即最终位置)
    placed: list[list[list[float]]] = [f.coords for f in sysdef.frags]
    scene.frag_names = [f.name for f in sysdef.frags]

    if channel == "incident_single":
        scene.notes.append("单片段通道:入射成员不放置(体系文件坐标原样)")
    elif i_proj < 0 or (channel == "incident_pair" and i_targ < 0):
        scene.notes.append("PROJECTILE/TARGET 角色未解析,预览以文件坐标原样显示")
    else:
        # ---- 束流代表值(抽样量取代表;固定量取输入值)----
        n_thta = integer("N_THTA", 0)
        thta_max = num("THTA_MAX")
        theta = thta_max  # 固定模式即输入值;采样模式以界值为代表值
        if n_thta == 1:
            scene.notes.append("N_THTA=1 采样极角,预览取 θ=THTA_MAX 为代表值")
        chi = num("CHI")
        if integer("N_CHI", 0) == 1:
            chi = 0.0
            scene.notes.append("N_CHI=1 采样方位角,预览取 χ=0 为代表值")
        b_max = num("B_MAX")
        n_b = integer("N_B", 0)
        if channel == "incident_surface":
            b_max = b_max if b_max > 0 else 0.0  # 表面不写 B_MAX(法定 0 = 定点束)
        b = b_max if n_b == 0 else 0.7 * b_max
        phi = chi if n_b == 0 else math.pi / 4
        if n_b == 1:
            scene.notes.append("N_B=1 碰撞盘采样,预览取 b=0.7·B_MAX、φ=45° 为代表值")

        n_aim = integer("N_AIM", 0)
        if n_aim == 1:
            if cell is not None:
                lx, ly, skew = box_from_cell(cell)
                a1 = (lx, 0.0)
                a2 = (ly * math.cos(skew), ly * math.sin(skew))
                drop = [0.5 * a1[0] + 0.5 * a2[0], 0.5 * a1[1] + 0.5 * a2[1]]
            else:
                drop = [0.0, 0.0]
            scene.notes.append("N_AIM=1 胞内均匀采样,预览取胞心为代表值")
        elif n_aim == 2:
            drop = [num("AIM_X"), num("AIM_Y")]
        else:
            drop = [0.0, 0.0]
        if channel == "incident_pair":
            drop = [0.0, 0.0]  # pair 落点恒原点(N_AIM∈{0}语义;2 在 pair 非法)

        # 速度(mode 1 定值;0/2 采样模式预览取 E_REL 对应速度为代表)
        m1 = sum(sysdef.frags[i_proj].masses)
        if channel == "incident_surface":
            mu = m1
        else:
            m2 = sum(sysdef.frags[i_targ].masses)
            mu = m1 * m2 / (m1 + m2)
        speed = math.sqrt(2.0 * e_rel / (mu * DYN_EV)) if (mu > 0 and e_rel > 0) else 0.0
        n_e_rel = integer("N_E_REL", 1)
        if n_e_rel != 1:
            scene.notes.append(f"N_E_REL={n_e_rel} 采样能速,预览取 E_REL 对应速度为代表值")

        sin_t = math.sin(theta)
        u_inc = [-sin_t * math.cos(chi), -sin_t * math.sin(chi), -math.cos(theta)]
        p_site = [
            drop[0] + b * math.cos(phi) + r_sep * sin_t * math.cos(chi),
            drop[1] + b * math.sin(phi) + r_sep * sin_t * math.sin(chi),
            r_sep,
        ]
        v_proj = [speed * c for c in u_inc]
        v_targ = None
        if channel == "incident_pair":
            mt = m1 + m2
            v_proj = [speed * (m2 / mt) * c for c in u_inc]
            v_targ = [-speed * (m1 / mt) * c for c in u_inc]

        placed[i_proj] = _shift_com(sysdef.frags[i_proj], p_site)
        if channel == "incident_pair":
            placed[i_targ] = _shift_com(sysdef.frags[i_targ], [drop[0], drop[1], 0.0])

        scene.beam = BeamView(
            theta=theta, chi=chi, b=b, phi=phi, drop=drop, speed=speed,
            u_inc=u_inc, p_site=p_site, mu=mu, v_proj=v_proj, v_targ=v_targ,
        )
        scene.b_max = b_max if channel == "incident_pair" else None
        aim = [drop[0] + b * math.cos(phi), drop[1] + b * math.sin(phi), 0.0]
        scene.aim_marker = aim

    # 汇全原子行
    off = 0
    for f, pos in zip(sysdef.frags, placed):
        scene.symbols.extend(f.symbols)
        scene.masses.extend(f.masses)
        scene.positions.extend(pos)
        scene.frag_slices.append((off, off + f.nat))
        scene.frag_coms.append(_com(Frag(file=f.file, name=f.name, symbols=f.symbols, coords=pos)))
        off += f.nat
    # 片段 COM 速度(默认零;入射通道赋值)
    scene.frag_vels = [[0.0, 0.0, 0.0] for _ in sysdef.frags]
    if scene.beam is not None:
        scene.frag_vels[i_proj] = scene.beam.v_proj
        if scene.beam.v_targ is not None:
            scene.frag_vels[i_targ] = scene.beam.v_targ
    return scene
