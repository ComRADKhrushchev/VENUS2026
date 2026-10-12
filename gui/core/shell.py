# -*- coding: utf-8 -*-
"""shell.py - bash/make 定位与 POSIX 路径换算(Windows 上的 WSL bash 陷阱)

设计:
  Windows 的 PATH 里 System32\\bash.exe(WSL)会排在 Git Bash 之前被
  shutil.which 命中,而 WSL 不认 C:/ 盘符路径——仓库构建(Fortran make/
  prebuild)走的是 Git Bash。find_bash() 优先探 Git Bash 的固定安装位,
  再退回 which 并排除 WSL 位。posix_path() 把 Windows 路径转正斜杠形式
  (Git Bash 接受 C:/...,且反斜杠会被 bash 当转义吃掉)。
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

_GIT_BASH_CANDIDATES = [
    r"C:\Program Files\Git\bin\bash.exe",
    r"C:\Program Files\Git\usr\bin\bash.exe",
    r"C:\Program Files (x86)\Git\bin\bash.exe",
    r"C:\Program Files (x86)\Git\usr\bin\bash.exe",
]
_WSL_MARKS = ("system32", "windowsapps")


def find_bash() -> str | None:
    """返回可用的 bash 解释器路径(Git Bash 优先);无则 None。"""
    env = os.environ.get("VENUS_GUI_BASH")
    if env and Path(env).is_file():
        return env
    for c in _GIT_BASH_CANDIDATES:
        if Path(c).is_file():
            return c
    w = shutil.which("bash")
    if w and not any(m in w.lower() for m in _WSL_MARKS):
        return w
    return None


def posix_path(p: str | Path) -> str:
    """Windows 路径 → Git Bash 可用的正斜杠形式。"""
    return str(p).replace("\\", "/")


def find_make() -> str | None:
    """make 的可用解释(Git Bash 内的 make;Windows 原生 make 亦可)。"""
    return shutil.which("make")
