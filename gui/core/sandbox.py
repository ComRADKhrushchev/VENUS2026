# -*- coding: utf-8 -*-
"""sandbox.py - 隔离树编译(临时副本里激活容器做真实 make,绝不触碰工作树)

设计:
  复制 src/(排除 build/) + scripts/ + systems/reference 到临时目录,把被检
  容器置 <tmp>/systems/<名>/ 后在 <tmp>/src 里跑 make(reference 与被检容器
  并存对编译合法——双绑 STOP 1 是运行期装配行为,drill 腿4 语义)。日志逐行
  流式回调并落盘 <tmp>/build.log(用户纪律:读数须有出处)。默认保留隔离树
  供诊断,cleanup() 显式清。
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .shell import find_bash, posix_path


@dataclass
class CompileResult:
    ok: bool
    returncode: int
    sandbox: Path      # 隔离树根(保留供诊断)
    log: Path          # build.log 落盘路径
    command: str       # 实际执行的 make 命令


def repo_root() -> Path:
    """gui/core/sandbox.py → 仓库根(向上二级到 gui,再一级)。"""
    return Path(__file__).resolve().parents[2]


def run_isolated_compile(
    container: str | Path,
    fc: str = "gfortran",
    root: Path | None = None,
    on_line=None,
    keep: bool = True,
    timeout: float = 600.0,
) -> CompileResult:
    """在隔离副本中激活 container 并 make;返回结果(日志已落盘)。

    on_line: 每行输出的回调(str);fc: 'gfortran' | 'ifx'。
    """
    bash = find_bash()
    if bash is None:
        raise RuntimeError("未找到 Git Bash(编译需要 bash+make;可用环境变量 VENUS_GUI_BASH 指定)")
    root = Path(root) if root else repo_root()
    container = Path(container).resolve()
    if not container.is_dir():
        raise FileNotFoundError(f"容器文件夹不存在: {container}")

    tmp = Path(tempfile.mkdtemp(prefix="venus_gui_sb_"))
    ignore = shutil.ignore_patterns("build", "*.o", "*.mod", "*.obj")
    shutil.copytree(root / "src", tmp / "src", ignore=ignore)
    shutil.copytree(root / "scripts", tmp / "scripts")
    (tmp / "systems").mkdir()
    shutil.copytree(root / "systems" / "reference", tmp / "systems" / "reference")
    shutil.copytree(container, tmp / "systems" / container.name)

    fc_arg = "" if fc in ("gfortran", "") else f" FC={fc}"
    cmd = f'cd "{posix_path(tmp / "src")}" && make{fc_arg}'
    log = tmp / "build.log"

    with log.open("w", encoding="utf-8", errors="replace", newline="\n") as lf:
        proc = subprocess.Popen(
            ["bash", "-c", cmd], executable=bash,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        )
        try:
            for line in proc.stdout:
                lf.write(line)
                lf.flush()
                if on_line:
                    on_line(line.rstrip("\n"))
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            lf.write(f"\n[sandbox] 超时({timeout}s),进程已终止\n")
            rc = -1
    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
    return CompileResult(ok=(rc == 0), returncode=rc, sandbox=tmp, log=log, command=cmd)


def cleanup(result: CompileResult) -> None:
    shutil.rmtree(result.sandbox, ignore_errors=True)
