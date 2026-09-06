# VENUS2026 完备化新程序 — src/ 构建入口（P0 骨架）
# 单次 ifx 调用完成编译+链接，产出 venus.e。
# 用法:
#   make          # 构建 venus.e
#   make clean    # 删除可执行文件与编译中间产物
# Linux: 需先加载 oneAPI 环境（source setvars.sh）使 ifx 可见。
#        新程序无 MKL 依赖，不移植旧库的 -qmkl 链接段。
# Windows: 下方路径变量均可用 make VAR=... 覆盖以适配本机 oneAPI/MSVC 版本。
# （平台段移植自 venus2026/Makefile，仅删 MKL 三库链接。）

IFX    ?= ifx
FFLAGS ?= -r8 -double-size=64 -i8 -O2 -w
EXE    := venus.e

ifeq ($(OS),Windows_NT)
INTEL_BASE ?= C:/Program Files (x86)/Intel/oneAPI
MSVC_BASE  ?= C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207
WINKIT_LIB ?= C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0
# Git Bash 的 /usr/bin/link 会遮蔽 MSVC link.exe；-Qlocation,link 直接指定链接器目录。
# 同时显式导出 LIB，使裸 Git Bash（未跑 setvars.bat）下也能链接系统库。
export PATH := $(INTEL_BASE)/compiler/2025.2/bin:$(PATH)
export LIB  := $(MSVC_BASE)/lib/x64;$(WINKIT_LIB)/um/x64;$(WINKIT_LIB)/ucrt/x64;$(INTEL_BASE)/compiler/2025.2/lib
LINKOPT := -Qlocation,link,"$(MSVC_BASE)/bin/Hostx64/x64"
else
LINKOPT :=
endif

# 源文件按模块依赖顺序排列（L5→L1→L2→L3→L0）。
# P0 骨架：SRC 暂为空列表；K3 阶段（新树独立编译）填入 L0-L5 文件，
# 填入前 make venus.e 仅验证 Makefile 解析，不产出可用可执行文件。
SRC :=

$(EXE): $(SRC)
	$(IFX) $(FFLAGS) $(LINKOPT) $(SRC) -o $(EXE)

clean:
	rm -f $(EXE) *.mod *.obj *.pdb *.ilk

.PHONY: clean
