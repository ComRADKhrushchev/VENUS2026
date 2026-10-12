# -*- coding: utf-8 -*-
"""keyspec ↔ src/input.f90 镜像同步测试(防漂移:键表/枚举/分派表逐一比对源文本)。"""
import re
from pathlib import Path

from conftest import ROOT

from gui.core import keyspec


def _src() -> str:
    return (ROOT / "src" / "input.f90").read_text(encoding="utf-8")


def test_key_table_mirror():
    """keyspec.KEY_TABLE 与 input.f90 key_tab(32 行)逐字段一致。"""
    text = _src()
    row = re.compile(
        r"key_spec_t\('(\w+)',\s*blk_(\w+),\s*kind_(\w+),\s*enum_(\w+),\s*"
        r"\.(true|false)\.,\s*'([^']*)'"
        r"(?:,\s*slotw=(\d+))?(?:,\s*cond=\.(true|false)\.)?"
    )
    rows = row.findall(text)
    assert len(rows) == 32, f"key_tab 应 32 行,解析得 {len(rows)}"
    assert len(keyspec.KEY_TABLE) == 32
    for fort, mine in zip(rows, keyspec.KEY_TABLE):
        name, blk, kind, enum, req, default, slotw, cond = fort
        assert mine.name == name
        assert mine.blk == ("flow" if blk == "flow" else "system"), name
        kind_map = {
            "int": "int", "real": "real", "char": "char", "log": "log",
            "level": "level", "enum": "enum", "elist": "elist", "posreal": "posreal",
        }
        assert mine.kind == kind_map[kind], name
        assert (mine.enum or "none") == enum, name
        assert mine.required == (req == "true"), name
        assert mine.default == default, name
        assert mine.slotw == int(slotw or 0), name
        assert mine.cond == (cond == "true"), name


def test_required_five():
    """必填键恰为 DT/SYSTEM_DIR/E_REL/R_SEP/DIST_SCHEME(input.f90:15)。"""
    req = {s.name for s in keyspec.KEY_TABLE if s.required}
    assert req == {"DT", "SYSTEM_DIR", "E_REL", "R_SEP", "DIST_SCHEME"}


def test_enums_mirror():
    """五个枚举注册表的词与序与 input.f90 一致。"""
    text = _src()
    for tbl, want in [
        ("task", keyspec.ENUM_TASK), ("integ", keyspec.ENUM_INTEG),
        ("dist", keyspec.ENUM_DIST), ("bath", keyspec.ENUM_BATH),
        ("spectrum", keyspec.ENUM_SPECTRUM),
    ]:
        m = re.search(
            rf"enum_row_t\), parameter :: {tbl}_tbl\(\d+\) = \[ &(.*?)\]",
            text, re.S,
        )
        assert m, tbl
        words = re.findall(r"enum_row_t\('(\w+)',\s*\d+\)", m.group(1))
        assert words == want, (tbl, words)


def test_dispatch_mirror():
    """MEMBER_TABLE 键与归属成员与 dispatch_tbl(35 行)一致(类型为 GUI 侧补注)。"""
    text = _src()
    rows = re.findall(r"dispatch_row_t\('(\w+)',\s*'([^']+)'", text)
    assert len(rows) == 35
    assert len(keyspec.MEMBER_TABLE) == 35
    for (key, member), mine in zip(rows, keyspec.MEMBER_TABLE):
        assert mine.key == key
        assert mine.member == member


def test_bare_defaults():
    """GWRITE_LEVEL 裸默认 2 / RANDOM_ORIENT 裸默认 F(E26/归档流法规)。"""
    d = {s.name: s.default for s in keyspec.KEY_TABLE}
    assert d["GWRITE_LEVEL"] == "2"
    assert d["RANDOM_ORIENT"] == "F"
    assert d["SPECTRUM_SOURCE"] == "COMPUTE"


def test_parse_scalar_domains():
    ps = keyspec.parse_scalar
    assert ps(keyspec.KIND_LEVEL, "", "2") == 2
    for bad in ("3", "7", "-1"):
        try:
            ps(keyspec.KIND_LEVEL, "", bad)
            raise AssertionError(bad)
        except ValueError:
            pass
    assert ps(keyspec.KIND_LOG, "", ".TRUE.") is True
    assert ps(keyspec.KIND_LOG, "", "f") is False
    assert ps(keyspec.KIND_ENUM, "integ", "radau") == "RADAU"
    assert ps(keyspec.KIND_ELIST, "dist", "boltzmann, ebk") == ["BOLTZMANN", "EBK"]
    try:
        ps(keyspec.KIND_ELIST, "dist", "INCIDENT")
        raise AssertionError
    except ValueError:
        pass
    assert ps(keyspec.KIND_POSREAL, "", "3.0") == 3.0
    try:
        ps(keyspec.KIND_POSREAL, "", "0.0")
        raise AssertionError
    except ValueError:
        pass
