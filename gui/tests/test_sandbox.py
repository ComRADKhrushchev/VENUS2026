# -*- coding: utf-8 -*-
"""sandbox 慢测(显式跑):auto_sample 隔离树 make 绿 + 生成槽位端到端编译。

默认不进常规批次:pytest gui/tests -m slow 才执行(整树编译 ~分钟级)。
"""
import shutil
import textwrap

import pytest
from conftest import ROOT

from gui.core import pes_gen, sandbox

pytestmark = pytest.mark.slow


def test_auto_sample_isolated_make():
    bash = sandbox.find_bash()
    if bash is None:
        pytest.skip("无 Git Bash")
    res = sandbox.run_isolated_compile(ROOT / "systems" / "inactive" / "auto_sample")
    try:
        assert res.ok, f"隔离编译失败,日志: {res.log}"
        assert (res.sandbox / "src" / "build" / "gfortran" / "venus.e").is_file()
    finally:
        # 日志留档到测试临时区之外?pytest tmp 由 pytest 管;此处直接保留在
        # sandbox 临时目录由 cleanup 删除——断言消息已给出日志路径供诊断
        sandbox.cleanup(res)


POT_WITH_VAR = textwrap.dedent("""\
    module my_pot
       implicit none
       private
       public :: pes, k_a, r_b
       real(8) :: k_a = 1.0d0
       real(8) :: r_b = 0.0d0
       logical :: r_b_found = .false.
    contains
       subroutine pes(q, mass, g, v)
          real(8), intent(in)  :: q(:), mass(:)
          real(8), intent(out) :: g(:), v
          v = 0.5d0*k_a*dot_product(q, q)
          g = k_a*q
       end subroutine pes
    end module my_pot
    """)


def test_generated_slots_compile_end_to_end(tmp_path):
    """pot-only 容器 + GUI 生成的 term/classify/keys+load → 隔离 make 绿。"""
    bash = sandbox.find_bash()
    if bash is None:
        pytest.skip("无 Git Bash")
    con = tmp_path / "mycon"
    con.mkdir()
    (con / "pot.f90").write_text(POT_WITH_VAR, encoding="utf-8")
    pes_gen.generate_gui_slots(
        con,
        term=pes_gen.TermSpec(variant="any_atom", r=4.5),
        classify=pes_gen.ClassifySpec(variant="const", label="bound"),
        keys_load=pes_gen.KeysLoadSpec(
            module="my_pot", bindings=[("K_A", "k_a"), ("R_B", "r_b")]
        ),
    )
    res = sandbox.run_isolated_compile(con)
    try:
        assert res.ok, f"生成槽位容器隔离编译失败,日志: {res.log}"
    finally:
        sandbox.cleanup(res)
