# -*- coding: utf-8 -*-
"""pes_gen.py - 缺失槽位代码生成器(Tier 1 容器可选槽的可行模板)

设计:
  用户澄清(2026-10-12):PES 没写的可选槽(如终止判定 term)由 GUI 按设定
  生成对应 Fortran 代码,仅覆盖可行模板——
    term     三变体: 任意原子逃逸 |q|>R / 前两片段 COM 间距>R / 弹丸 z>R
    classify 两变体: COM 距阈值 → 'bound'/'scattered' / 常量标签
    keys+load        键词表声明 + pull 接线(目标 = 作者模块的 real(8) 公共变量)
  方法组五名(prop/match/sync/init/method_name)不生成——真实物理;绝热面
  缺省即可跑。产物 = 容器文件夹内单文件 gui_slots.f90(单模块 gui_slots,
  删除该文件即撤销);生成前扫描防撞名,生成后应重扫确认。
  参照: auto_sample 范本(systems/inactive/auto_sample)与 reg_reference.f90
  的 pull 接线形。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from .pes_scan import scan_container, ScanResult

GEN_FILE = "gui_slots.f90"
GEN_MODULE = "gui_slots"

TERM_VARIANTS = {
    "any_atom": "任意原子坐标分量越出 |q|>R 球(原点逃逸,auto_sample 同款)",
    "frag_com": "前两片段 COM 间距 > R(气相对分离判据)",
    "proj_z": "弹丸片段 COM 高度 z > R(表面反射/脱附判据)",
}
CLASSIFY_VARIANTS = {
    "com_dist": "前两片段 COM 距 < R → 'bound',否则 'scattered'",
    "const": "常量标签(接线练习桩)",
}


@dataclass
class TermSpec:
    variant: str        # TERM_VARIANTS 键
    r: float            # 判据半径 [Å]


@dataclass
class ClassifySpec:
    variant: str        # CLASSIFY_VARIANTS 键
    r: float = 0.0      # com_dist 阈值 [Å]
    label: str = "bound"


@dataclass
class KeysLoadSpec:
    """keys()+load() 生成参数:键名 → 作者模块的 real(8) 公共变量。"""
    module: str                     # 目标模块名(如 pes 所在模块)
    bindings: list[tuple[str, str]] = field(default_factory=list)  # (键名, 变量名)


def _fmt_r(r: float) -> str:
    return f"{r:.16g}d0"


def _term_body(spec: TermSpec) -> str:
    r = _fmt_r(spec.r)
    if spec.variant == "any_atom":
        return f"""   logical function term(sta)
      use state, only: state_t
      implicit none
      type(state_t), intent(in) :: sta
      term = any(abs(sta%q) > {r})
   end function term"""
    if spec.variant == "frag_com":
        return f"""   logical function term(sta)
      use state, only: state_t
      use config_atoms, only: list_atoms
      implicit none
      type(state_t), intent(in) :: sta
      real(8) :: com_a(3), com_b(3), dv(3)
      integer :: k, a
      term = .false.
      if (size(list_atoms%frag) < 2) return
      com_a = 0.0d0
      do k = 1, list_atoms%frag(1)%nat
         a = list_atoms%frag(1)%list(k)
         com_a = com_a + list_atoms%mass(a)*sta%q(3*a - 2:3*a)
      end do
      com_a = com_a/sum(list_atoms%mass(list_atoms%frag(1)%list))
      com_b = 0.0d0
      do k = 1, list_atoms%frag(2)%nat
         a = list_atoms%frag(2)%list(k)
         com_b = com_b + list_atoms%mass(a)*sta%q(3*a - 2:3*a)
      end do
      com_b = com_b/sum(list_atoms%mass(list_atoms%frag(2)%list))
      dv = com_a - com_b
      term = sqrt(dot_product(dv, dv)) > {r}
   end function term"""
    if spec.variant == "proj_z":
        return f"""   logical function term(sta)
      use state, only: state_t
      use config_atoms, only: list_atoms
      implicit none
      type(state_t), intent(in) :: sta
      real(8) :: com(3)
      integer :: k, a
      term = .false.
      if (list_atoms%i_proj < 1 .or. list_atoms%i_proj > size(list_atoms%frag)) return
      com = 0.0d0
      do k = 1, list_atoms%frag(list_atoms%i_proj)%nat
         a = list_atoms%frag(list_atoms%i_proj)%list(k)
         com = com + list_atoms%mass(a)*sta%q(3*a - 2:3*a)
      end do
      com = com/sum(list_atoms%mass(list_atoms%frag(list_atoms%i_proj)%list))
      term = com(3) > {r}
   end function term"""
    raise ValueError(f"未知 term 变体: {spec.variant}")


def _classify_body(spec: ClassifySpec) -> str:
    if spec.variant == "com_dist":
        r2 = _fmt_r(spec.r * spec.r)
        return f"""   function classify(fin) result(name)
      use final_state, only: fin_t
      use config_atoms, only: list_atoms
      implicit none
      type(fin_t), intent(in) :: fin
      character(len=32) :: name
      real(8) :: com_a(3), com_b(3), dv(3)
      integer :: k, a
      name = 'scattered'
      if (size(list_atoms%frag) < 2) then
         name = 'bound'
         return
      end if
      com_a = 0.0d0
      do k = 1, list_atoms%frag(1)%nat
         a = list_atoms%frag(1)%list(k)
         com_a = com_a + list_atoms%mass(a)*fin%q_raw(3*a - 2:3*a)
      end do
      com_a = com_a/sum(list_atoms%mass(list_atoms%frag(1)%list))
      com_b = 0.0d0
      do k = 1, list_atoms%frag(2)%nat
         a = list_atoms%frag(2)%list(k)
         com_b = com_b + list_atoms%mass(a)*fin%q_raw(3*a - 2:3*a)
      end do
      com_b = com_b/sum(list_atoms%mass(list_atoms%frag(2)%list))
      dv = com_a - com_b
      if (dot_product(dv, dv) < {r2}) name = 'bound'
   end function classify"""
    if spec.variant == "const":
        label = spec.label[:32]
        return f"""   function classify(fin) result(name)
      use final_state, only: fin_t
      implicit none
      type(fin_t), intent(in) :: fin
      character(len=32) :: name
      name = '{label}'
   end function classify"""
    raise ValueError(f"未知 classify 变体: {spec.variant}")


def _keys_load_body(spec: KeysLoadSpec) -> tuple[list[str], str]:
    """返回 (use 行列表, contains 体)。"""
    keys = [k for k, _ in spec.bindings]
    vars_ = [v for _, v in spec.bindings]
    use_lines = [
        "   use input, only: pull",
        f"   use {spec.module}, only: " + ", ".join(vars_),
    ]
    arr = ", ".join(f"'{k}'" for k in keys)
    body = f"""   function keys() result(kk)
      character(len=16) :: kk({len(keys)})
      kk = [character(len=16) :: {arr}]
   end function keys
   subroutine load()
""" + "\n".join(
        f"      call pull('{k}', {v})" for k, v in spec.bindings
    ) + """
   end subroutine load"""
    return use_lines, body


def generate_gui_slots(
    folder: str | Path,
    term: TermSpec | None = None,
    classify: ClassifySpec | None = None,
    keys_load: KeysLoadSpec | None = None,
) -> Path:
    """生成/覆盖 <folder>/gui_slots.f90;被请求槽位若已存在(非本文件)则拒绝。

    返回写入路径;ValueError 说明冲突(调用方重扫后可定位)。
    """
    if term is None and classify is None and keys_load is None:
        raise ValueError("未请求任何槽位")
    scan: ScanResult = scan_container(folder)
    for slot in ("term", "classify", "keys", "load"):
        requested = (slot == "term" and term is not None) or \
                    (slot == "classify" and classify is not None) or \
                    (slot in ("keys", "load") and keys_load is not None)
        if not requested:
            continue
        existing = scan.found(slot)
        if existing is not None and existing.file != GEN_FILE:
            raise ValueError(
                f"槽位 '{slot}' 已由 {existing.file} 提供(生成会造成重定义——"
                f"如需 GUI 版本请先移除原定义)"
            )
    # 既有 gui_slots.f90 的旧定义视同不存在(整文件覆盖)
    scan.defs = {
        n: [d for d in lst if d.file != GEN_FILE] for n, lst in scan.defs.items()
    }

    use_lines: list[str] = ["   implicit none", "   private"]
    publics: list[str] = []
    bodies: list[str] = []
    if term is not None:
        publics.append("term")
        bodies.append(_term_body(term))
    if classify is not None:
        publics.append("classify")
        bodies.append(_classify_body(classify))
    if keys_load is not None:
        if not keys_load.bindings:
            raise ValueError("keys+load 需要至少一条键-变量绑定")
        if len({k for k, _ in keys_load.bindings}) != len(keys_load.bindings):
            raise ValueError("键名重复")
        for k, _ in keys_load.bindings:
            if not k or len(k) > 32:
                raise ValueError(f"键名非法(空或超 32 字符): {k}")
        kl_use, kl_body = _keys_load_body(keys_load)
        use_lines = kl_use + use_lines
        publics.extend(["keys", "load"])
        bodies.append(kl_body)

    head = [
        "!=====================================================================",
        "! gui_slots.f90 - VENUS2026 GUI 生成的容器槽位补件(整文件覆盖式)",
        "!   删除本文件即撤销全部生成内容;手工修改会在下次生成时丢失。",
        "!=====================================================================",
        f"module {GEN_MODULE}",
    ] + use_lines
    head.append("   public :: " + ", ".join(publics))
    text = "\n".join(head) + "\ncontains\n" + "\n\n".join(bodies) + f"\nend module {GEN_MODULE}\n"
    p = Path(folder) / GEN_FILE
    p.write_text(text, encoding="utf-8", newline="\n")
    return p
