# -*- coding: utf-8 -*-
"""pes_gen 测试:pot-only 容器生成槽位后重扫绿、撞名拒绝、覆盖语义。"""
import pytest

from gui.core import pes_gen, pes_scan
from test_pes_scan import POT_ONLY, _pot_folder


def test_generate_term_and_rescan(tmp_path):
    d = _pot_folder(tmp_path, {"pot.f90": POT_ONLY})
    p = pes_gen.generate_gui_slots(d, term=pes_gen.TermSpec(variant="any_atom", r=2.5))
    assert p.name == "gui_slots.f90"
    res = pes_scan.scan_container(d)
    assert res.ok, res.errors
    assert res.found("term") is not None
    assert res.found("term").file == "gui_slots.f90"
    text = p.read_text(encoding="utf-8")
    assert "any(abs(sta%q) > 2.5d0)" in text


def test_generate_all_three(tmp_path):
    d = _pot_folder(tmp_path, {"pot.f90": POT_ONLY})
    pes_gen.generate_gui_slots(
        d,
        term=pes_gen.TermSpec(variant="proj_z", r=6.0),
        classify=pes_gen.ClassifySpec(variant="com_dist", r=1.2),
        keys_load=pes_gen.KeysLoadSpec(
            module="my_pot", bindings=[("K_A", "k_a")]
        ),
    )
    # keys+load 目标变量须在作者模块在场——补一个变量再生成合法文件
    text = (d / "pot.f90").read_text(encoding="utf-8")
    (d / "pot.f90").write_text(
        text.replace("public :: pes", "public :: pes, k_a\n   real(8) :: k_a = 1.0d0"),
        encoding="utf-8",
    )
    pes_gen.generate_gui_slots(
        d,
        term=pes_gen.TermSpec(variant="proj_z", r=6.0),
        classify=pes_gen.ClassifySpec(variant="com_dist", r=1.2),
        keys_load=pes_gen.KeysLoadSpec(module="my_pot", bindings=[("K_A", "k_a")]),
    )
    res = pes_scan.scan_container(d)
    assert res.ok, res.errors
    for slot in ("term", "classify", "keys", "load"):
        assert res.found(slot) is not None, slot
        assert res.found(slot).file == "gui_slots.f90"
    # keys 词表提取到生成键
    assert "K_A" in res.keys_words
    body = (d / "gui_slots.f90").read_text(encoding="utf-8")
    assert "call pull('K_A', k_a)" in body


def test_conflict_existing_term(tmp_path):
    text = POT_ONLY.replace("public :: pes", "public :: pes, term").replace(
        "end module my_pot",
        "   logical function term(sta)\n      use state, only: state_t\n"
        "      type(state_t), intent(in) :: sta\n      term = .false.\n"
        "   end function term\nend module my_pot",
    )
    d = _pot_folder(tmp_path, {"pot.f90": text})
    with pytest.raises(ValueError, match="term"):
        pes_gen.generate_gui_slots(d, term=pes_gen.TermSpec("any_atom", 1.0))


def test_overwrite_semantics(tmp_path):
    """整文件覆盖式:只请求 term 的再生成会移除先前的 classify(文档化语义)。"""
    d = _pot_folder(tmp_path, {"pot.f90": POT_ONLY})
    pes_gen.generate_gui_slots(d, classify=pes_gen.ClassifySpec("const", label="ok"))
    assert pes_scan.scan_container(d).found("classify") is not None
    pes_gen.generate_gui_slots(d, term=pes_gen.TermSpec("any_atom", 3.0))
    res = pes_scan.scan_container(d)
    assert res.ok
    assert res.found("term") is not None
    assert res.found("classify") is None


def test_empty_request_rejected(tmp_path):
    d = _pot_folder(tmp_path, {"pot.f90": POT_ONLY})
    with pytest.raises(ValueError):
        pes_gen.generate_gui_slots(d)
