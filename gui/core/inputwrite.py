# -*- coding: utf-8 -*-
"""inputwrite.py - input_qct.txt 生成器 + 语法镜像解析器

设计:
  写出键序 = key_tab 序(块一→块二)→ 成员参数键(dispatch 序)→ 容器键;
  行格式 `KEY = VALUE`(与在档案例同款;等号两侧各一空格)。默认值省略
  (可切换全写出);条件键仅写出用户实际填写的行。
  parse_input_text 逐字镜像 input.f90 parse_file(注释/空行跳过、'=' 切分、
  键大写、值 adjustl、行尾 ',' 续行吸收且逗号保留为分隔符),供往返与
  金标对比测试——它不是运行时权威,权威永远在 Fortran 解析器。
"""

from __future__ import annotations

from pathlib import Path

from .keyspec import KEY_TABLE, MEMBER_TABLE, KEY_NAMES, MEMBER_KEY_NAMES


def write_input(
    values: dict[str, str],
    container_keys: dict[str, str] | None = None,
    header_lines: list[str] | None = None,
    include_defaults: bool = False,
) -> str:
    """把键值字典写成 input_qct.txt 文本。

    values: 程序键+成员键的 KEY→值字符串(值原样,已由 validate 校验);
    container_keys: 容器键 KEY→值(追加在成员键之后);
    header_lines: 头部 '#' 注释行(默认一行生成器标记);
    include_defaults: True 时非必填非条件键的缺省以默认值写出。
    """
    values = {k.strip().upper(): (v.strip() if isinstance(v, str) else v) for k, v in values.items()}
    container_keys = {k.strip().upper(): v.strip() for k, v in (container_keys or {}).items()}
    known = KEY_NAMES | MEMBER_KEY_NAMES | set(container_keys)
    unknown = set(values) - known
    if unknown:
        raise ValueError(f"values 含未知键: {sorted(unknown)}")

    lines: list[str] = []
    if header_lines is None:
        lines.append("# 由 VENUS2026 GUI 生成(input_qct.txt)")
    else:
        lines.extend(header_lines)

    for spec in KEY_TABLE:
        v = values.get(spec.name, "")
        if not v and include_defaults and not spec.required and not spec.cond:
            v = spec.default
        if v:
            lines.append(f"{spec.name} = {v}")
    for mk in MEMBER_TABLE:
        v = values.get(mk.key, "")
        if v:
            lines.append(f"{mk.key} = {v}")
    for k in sorted(container_keys):
        v = container_keys[k]
        if v:
            lines.append(f"{k} = {v}")
    return "\n".join(lines) + "\n"


def write_input_file(
    path: str | Path,
    values: dict[str, str],
    container_keys: dict[str, str] | None = None,
    header_lines: list[str] | None = None,
    include_defaults: bool = False,
) -> Path:
    p = Path(path)
    p.write_text(
        write_input(values, container_keys, header_lines, include_defaults),
        encoding="utf-8",
        newline="\n",
    )
    return p


def parse_input_text(text: str) -> dict[str, str]:
    """镜像解析(测试/对比用):KEY→值(键大写、值 strip);错误抛 ValueError。"""
    out: dict[str, str] = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].lstrip()
        kln = i + 1
        i += 1
        if not line.strip() or line[0] in "#!":
            continue
        eq = line.find("=")
        if eq < 0:
            raise ValueError(f"第 {kln} 行: 不是 KEYWORD=VALUE 行: {line}")
        key = line[:eq].strip().upper()
        if not key:
            raise ValueError(f"第 {kln} 行: '=' 前键名为空")
        vals = line[eq + 1:].lstrip()
        # 行尾 ',' 吸收续行(逗号保留为分隔符——input.f90 trim(vals)//' '//line)
        while vals.strip() and vals.strip().endswith(","):
            if i >= len(lines):
                raise ValueError(f"键 {key}(第 {kln} 行): 续行越过文件末尾")
            vals = vals.rstrip() + " " + lines[i].lstrip()
            i += 1
        vals = vals.strip()
        if len(vals) > 512:
            raise ValueError(f"键 {key}(第 {kln} 行): 值超 512 字符")
        if key in out:
            raise ValueError(f"重复键 {key}")
        out[key] = vals
    return out


def parse_input_file(path: str | Path) -> dict[str, str]:
    return parse_input_text(Path(path).read_text(encoding="utf-8", errors="replace"))


def apply_defaults(parsed: dict[str, str]) -> dict[str, str]:
    """注入默认值(inject_defaults 移植:非必填、非条件、未写的键取默认)。"""
    out = dict(parsed)
    for spec in KEY_TABLE:
        if spec.required or spec.cond:
            continue
        if spec.name not in out:
            out[spec.name] = spec.default
    return out
