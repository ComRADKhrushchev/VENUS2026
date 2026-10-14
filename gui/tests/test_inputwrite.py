# -*- coding: utf-8 -*-
"""inputwrite 测试:在档案例金标复现(键值语义一致)+ 往返 + 未知键拒绝。"""
from pathlib import Path

import pytest
from conftest import ROOT

from gui.core import inputwrite
from gui.core.keyspec import REFERENCE_CONTAINER_KEYS


def _golden(case: str, container_keys: set[str]):
    """解析在档 input_qct.txt → 重写 → 双方注入默认后键值全等。"""
    path = ROOT / "cases" / case / "input_qct.txt"
    parsed = inputwrite.parse_input_file(path)
    prog = {k: v for k, v in parsed.items() if k not in container_keys}
    cont = {k: v for k, v in parsed.items() if k in container_keys}
    text = inputwrite.write_input(prog, cont)
    reparsed = inputwrite.parse_input_text(text)
    a = inputwrite.apply_defaults(parsed)
    b = inputwrite.apply_defaults(reparsed)
    # 容器键不注入默认,直接并入比较集
    for k, v in cont.items():
        a[k] = v
        b[k] = v
    assert a == b, (
        f"{case}: 差异 "
        f"{ {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)} }"
    )
    return text


def test_golden_harmonic_scatter():
    text = _golden("case_harmonic_scatter", set(REFERENCE_CONTAINER_KEYS))
    assert "K_HARM = 0.04" in text
    assert "Q_TERM = 9.5" in text


def test_golden_h3_ebk():
    _golden("case_spectrum_source/h3_ebk", set())


def test_golden_b_layer_andersen():
    _golden("case_b_layer/andersen", set(REFERENCE_CONTAINER_KEYS))


def test_golden_surface_window():
    _golden("case_surface_window", set(REFERENCE_CONTAINER_KEYS))


def test_roundtrip_and_continuation():
    text = "A = 1\nB = x,\n  y\n"
    parsed = inputwrite.parse_input_text(text)
    assert parsed["B"] == "x, y"


def test_unknown_key_rejected():
    with pytest.raises(ValueError):
        inputwrite.write_input({"MODEL": "X"})


def test_include_defaults():
    text = inputwrite.write_input({"DT": "0.05"}, include_defaults=True)
    assert "GWRITE_LEVEL = 2" in text
    assert "RANDOM_ORIENT = F" in text
    assert "PROJECTILE" not in text  # 条件键永不注入默认
