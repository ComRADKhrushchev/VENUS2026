# -*- coding: utf-8 -*-
"""pytest 入口:把仓库根放入 sys.path,使 gui.* 可导入;给出 REPO_ROOT。"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
