# -*- coding: utf-8 -*-
"""pes_scan 测试:auto_sample 全绿 + drill 负面变异 + 对真实 prebuild.sh 的差分核对。"""
import shutil
import subprocess
from pathlib import Path

import pytest
from conftest import ROOT

from gui.core import pes_scan

AUTO_SAMPLE = ROOT / "systems" / "inactive" / "auto_sample"


def test_auto_sample_full_green():
    res = pes_scan.scan_container(AUTO_SAMPLE)
    assert res.is_tier1 and res.has_f90
    assert not res.errors, res.errors
    for name in pes_scan.ALL_NAMES:
        assert res.found(name) is not None, name
    assert res.found("pes").kind == "s"
    assert res.found("term").kind == "f"
    assert res.found("pes").module == "auto_sample_pot"
    assert res.found("load").module == "auto_sample_io"
    assert res.found("method_name").module == "auto_sample_met"
    assert res.prop_match_pair
    assert res.keys_words == ["K_SPRING", "R_TERM"]
    assert res.method_name_value == "auto_tdhf"
    assert "k_spring" in res.module_vars.get("auto_sample_pot", [])


POT_ONLY = """\
module my_pot
   implicit none
   private
   public :: pes
contains
   subroutine pes(q, mass, g, v)
      real(8), intent(in)  :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
      v = 0.0d0
      g = 0.0d0
   end subroutine pes
end module my_pot
"""


def _pot_folder(tmp_path: Path, files: dict[str, str]) -> Path:
    d = tmp_path / "mycon"
    d.mkdir()
    for name, text in files.items():
        (d / name).write_text(text, encoding="utf-8")
    return d


def test_pot_only_green(tmp_path):
    d = _pot_folder(tmp_path, {"pot.f90": POT_ONLY})
    res = pes_scan.scan_container(d)
    assert res.ok, res.errors
    assert res.found("pes").module == "my_pot"
    assert res.found("term") is None
    assert not res.prop_match_pair


def test_reg_file_is_ordinary_source(tmp_path):
    """reg 形态 2026-10-14 退役:reg_*.f90 只是普通源文件,不改变容器判读。"""
    d = _pot_folder(tmp_path, {"pot.f90": POT_ONLY, "reg_x.f90": "module reg_x\nend module\n"})
    res = pes_scan.scan_container(d)
    assert res.is_tier1
    assert res.ok, res.errors
    assert res.found("pes").module == "my_pot"


def test_duplicate_pes(tmp_path):
    extra = (
        "module extra\ncontains\n"
        "   subroutine pes(q, mass, g, v)\n"
        "      real(8), intent(in) :: q(:), mass(:)\n"
        "      real(8), intent(out) :: g(:), v\n"
        "   end subroutine pes\n"
        "end module\n"
    )
    d = _pot_folder(tmp_path, {"a.f90": POT_ONLY, "b.f90": extra})
    res = pes_scan.scan_container(d)
    assert any("'pes'" in e and "定义了 2 次" in e for e in res.errors)


def test_wrong_form_keys(tmp_path):
    text = POT_ONLY.replace("public :: pes", "public :: pes, keys").replace(
        "contains",
        "contains\n   subroutine keys()\n   end subroutine keys",
    )
    d = _pot_folder(tmp_path, {"pot.f90": text})
    res = pes_scan.scan_container(d)
    assert any("'keys'" in e and "须为 function" in e for e in res.errors)


def test_prop_without_match(tmp_path):
    text = POT_ONLY.replace(
        "end module my_pot",
        """   subroutine prop(sta, dt)
      use state, only: state_t
      type(state_t), intent(inout) :: sta
      real(8), intent(in) :: dt
   end subroutine prop
end module my_pot""",
    ).replace("public :: pes", "public :: pes, prop")
    d = _pot_folder(tmp_path, {"pot.f90": text})
    res = pes_scan.scan_container(d)
    assert any("prop 与 match" in e for e in res.errors)


def test_lone_method_name(tmp_path):
    text = POT_ONLY.replace(
        "end module my_pot",
        """   function method_name() result(s)
      character(len=32) :: s
      s = 'x_method'
   end function method_name
end module my_pot""",
    ).replace("public :: pes", "public :: pes, method_name")
    d = _pot_folder(tmp_path, {"pot.f90": text})
    res = pes_scan.scan_container(d)
    assert any("method_name" in e for e in res.errors)


def test_no_pes_wires_nothing(tmp_path):
    d = _pot_folder(tmp_path, {"pot.f90": "module x\n   implicit none\nend module x\n"})
    res = pes_scan.scan_container(d)
    assert any("pes" in e for e in res.errors)
    assert not res.ok


def test_scan_edge_forms(tmp_path):
    """前缀形计入、end 排除、注释尾剥除、模块外定义报错。"""
    text = """\
module m
   implicit none
   private
   public :: pes, term
contains
   subroutine pes(q, mass, g, v)
      real(8), intent(in) :: q(:), mass(:)
      real(8), intent(out) :: g(:), v
   end subroutine pes
   logical function term(sta)
      use state, only: state_t
      type(state_t), intent(in) :: sta
      term = .false.
   end function term
end module m
! subroutine ghost(q)  -- 注释尾,不计
subroutine outside(q)
   real(8), intent(in) :: q(:)
end subroutine outside
"""
    d = _pot_folder(tmp_path, {"pot.f90": text})
    res = pes_scan.scan_container(d)
    # outside 不在预留名内 → 不计;term 前缀形已计入
    assert res.found("term") is not None
    assert res.found("term").kind == "f"
    assert res.ok


# ---- 差分核对:与真实 prebuild.sh 对同一组文件夹的判定一致(bash 可用时)----

_VARIANTS = {
    "valid": {"pot.f90": POT_ONLY},
    "pes_dup": {"a.f90": POT_ONLY, "b.f90": POT_ONLY},
    "bad_keys": {"pot.f90": POT_ONLY.replace("public :: pes", "public :: pes, keys").replace(
        "contains", "contains\n   subroutine keys()\n   end subroutine keys")},
    "no_pes": {"pot.f90": "module x\nend module x\n"},
    "lone_sync": {"pot.f90": POT_ONLY.replace(
        "end module my_pot",
        "   subroutine sync(s)\n      use state, only: elec_t\n      type(elec_t), intent(inout) :: s\n   end subroutine sync\nend module my_pot",
    ).replace("public :: pes", "public :: pes, sync")},
}


@pytest.mark.parametrize("name", sorted(_VARIANTS))
def test_prebuild_differential(name, tmp_path):
    from gui.core.shell import find_bash, posix_path

    bash = find_bash()
    if bash is None:
        pytest.skip("无 Git Bash,跳过与 prebuild.sh 的差分核对")
    tree = tmp_path / f"t_{name}"
    src = tree / "src"
    src.mkdir(parents=True)
    (tree / "bld").mkdir()
    con = tree / "systems" / "mycon"
    con.mkdir(parents=True)
    for fn, text in _VARIANTS[name].items():
        (con / fn).write_text(text, encoding="utf-8")
    prebuild = ROOT / "scripts" / "prebuild.sh"
    proc = subprocess.run(
        ["bash", posix_path(prebuild), posix_path(src), posix_path(tree / "bld")],
        capture_output=True, text=True, timeout=60, executable=bash,
    )
    mine = pes_scan.scan_container(con)
    if name == "valid":
        assert proc.returncode == 0, proc.stderr
        assert mine.ok, mine.errors
    else:
        assert proc.returncode == 1, (name, proc.stdout)
        assert not mine.ok, (name, mine.errors)
        # 同名提及:prebuild 的英文中止名 vs 我的中文错误里同名出现
        msg = (proc.stderr or "") + (proc.stdout or "")
        for word in ("pes", "keys", "sync", "prop", "method_name", "wires nothing"):
            if word in msg:
                assert any(word in e for e in mine.errors), (name, word, mine.errors)
                break
