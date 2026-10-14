# -*- coding: utf-8 -*-
"""validate 测试:三范式的正例 + 各族 fatal 规则的负例(写盘前静态预检)。"""
from conftest import ROOT

from gui.core import sysdir, validate


def _pair_sd():
    return sysdir.scan_system_dir(ROOT / "cases" / "case_harmonic_scatter" / "system")


def _surf_sd():
    return sysdir.scan_system_dir(ROOT / "cases" / "case_spectrum_source" / "h2ag_ebk" / "system")


def _pair_values():
    return {
        "SYSTEM_DIR": "system", "E_REL": "200.0", "R_SEP": "8.0", "DT": "0.05",
        "DIST_SCHEME": "THERMALIZE[T_EQ=300.0, N_EQ=25, DT_EQ=0.05], ROTATION[J_ROT=3]",
        "PROJECTILE": "frag_a", "TARGET": "frag_b",
        "B_MAX": "2.0", "N_B": "1",
    }


def test_pair_positive():
    issues = validate.validate(_pair_values(), _pair_sd(), ["K_HARM"])
    assert issues == [], [str(i) for i in issues]


def test_missing_required():
    v = _pair_values()
    del v["DT"]
    issues = validate.validate(v, _pair_sd())
    assert any(i.key == "DT" for i in issues)


def test_surface_rules():
    sd = _surf_sd()
    base = {
        "SYSTEM_DIR": "system", "E_REL": "200.0", "R_SEP": "8.0", "DT": "0.05",
        "DIST_SCHEME": "EBK[N_VIB=1, N_ROT=2]", "PROJECTILE": "H2",
        "SCATTERING": "SURFACE", "CELL": "2.9437 0 0, 1.47185 2.5493 0, 0 0 20",
        "N_B": "0", "N_AIM": "1",
    }
    assert validate.validate(base, sd) == []
    # 写 B_MAX → 拦截(2026-10-01 裁决)
    bad = dict(base, B_MAX="2.0")
    assert any(i.key == "B_MAX" for i in validate.validate(bad, sd))
    # 写 TARGET → 拦截
    bad = dict(base, TARGET="H2")
    assert any(i.key == "TARGET" for i in validate.validate(bad, sd))
    # 缺 CELL
    bad = dict(base)
    del bad["CELL"]
    issues = validate.validate(bad, sd)
    assert any(i.key == "CELL" for i in issues)
    # 气相写 CELL → 拦截
    bad = dict(_pair_values(), CELL="3 0 0, 0 3 0, 0 0 9")
    assert any(i.key == "CELL" for i in validate.validate(bad, _pair_sd()))


def test_role_rules():
    sd = _pair_sd()
    v = _pair_values()
    del v["TARGET"]
    assert any(i.key == "TARGET" for i in validate.validate(v, sd))
    v = dict(_pair_values(), PROJECTILE="frag_b", TARGET="frag_b")
    assert any("不同片段" in i.msg for i in validate.validate(v, sd))
    v = dict(_pair_values(), PROJECTILE="nope", TARGET="frag_b")
    assert any(i.key == "PROJECTILE" for i in validate.validate(v, sd))


def test_dist_scheme_length():
    v = dict(_pair_values(), DIST_SCHEME="THERMALIZE")
    assert any(i.key == "DIST_SCHEME" and "片段数" in i.msg
               for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), DIST_SCHEME="THERMALIZE,INCIDENT")
    assert any(i.key == "DIST_SCHEME" for i in validate.validate(v, _pair_sd()))


def test_beam_rules():
    v = dict(_pair_values(), B_MAX="0.0", N_B="1")
    assert any(i.key == "N_B" for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), B_MAX="0.0", N_B="0")
    assert any(i.key == "B_MAX" for i in validate.validate(v, _pair_sd()))  # pair 须 >0
    v = dict(_pair_values(), N_AIM="2", AIM_X="1.0", AIM_Y="1.0")
    assert any("表面" in i.msg for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), AIM_X="1.0")
    assert any(i.key == "AIM_X" for i in validate.validate(v, _pair_sd()))


def test_bath_and_spectrum_rules():
    v = dict(_pair_values(), BATH="ANDERSEN")
    assert any(i.key in ("T_BATH", "NU_COLL") for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), N_BATH="10", DT_BATH="0.0")
    assert any(i.key == "DT_BATH" for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), SPECTRUM_FILE="spec.tbl")
    assert any(i.key == "SPECTRUM_FILE" for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), SPECTRUM_SOURCE="MANUAL")
    assert any(i.key == "SPECTRUM_FILE" for i in validate.validate(v, _pair_sd()))


def test_member_and_unknown_and_container_keys():
    v = dict(_pair_values(), DIST_SCHEME="BOLTZMANN,ROTATION", T_VIB="0.0")
    assert any(i.key == "T_VIB" for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values(), MODEL="X")
    assert any(i.key == "MODEL" for i in validate.validate(v, _pair_sd()))
    v = dict(_pair_values())
    out = validate.validate(v, _pair_sd(), ["DT"])
    assert any(i.key == "DT" and "撞名" in i.msg for i in out)
    out = validate.validate(v, _pair_sd(), ["X" * 33])
    assert any("32 字符" in i.msg for i in out)


def test_integrator_order_domain():
    """积分器×阶合法域(权威 src/dynamics/propagator.f90:126-165/186/223-224):
    VERLET 1/2·SYMPLECTIC 4/6/8·RADAU 1-15; SYMPLECTIC 缺省(默认1)会 fatal, 须拦。"""
    v = dict(_pair_values(), INTEGRATOR="SYMPLECTIC")
    assert any(i.key == "ORDER" for i in validate.validate(v, _pair_sd()))
    for ok in ("4", "6", "8"):
        assert not any(i.key == "ORDER" for i in validate.validate(dict(v, ORDER=ok), _pair_sd()))
    assert any(i.key == "ORDER" for i in validate.validate(dict(v, ORDER="5"), _pair_sd()))
    v2 = dict(_pair_values(), INTEGRATOR="VERLET", ORDER="3")
    assert any(i.key == "ORDER" for i in validate.validate(v2, _pair_sd()))
    v3 = dict(_pair_values(), INTEGRATOR="VERLET")  # 缺省 → 默认 1 合法
    assert not any(i.key == "ORDER" for i in validate.validate(v3, _pair_sd()))
    v4 = dict(_pair_values(), INTEGRATOR="RADAU", ORDER="15")
    assert not any(i.key == "ORDER" for i in validate.validate(v4, _pair_sd()))
    v5 = dict(_pair_values(), INTEGRATOR="RADAU", ORDER="16")
    assert any(i.key == "ORDER" for i in validate.validate(v5, _pair_sd()))
