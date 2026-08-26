# VENUS2026 — TEST 势最小构建（唯一构建入口）
# 单次 ifx 调用完成编译+链接，产出 venus_test.e。
# 用法:
#   make          # 构建 venus_test.e
#   make clean    # 删除可执行文件与编译中间产物
# Linux: 需先加载 oneAPI 环境（source setvars.sh）使 ifx 可见，MKL 经 -qmkl 链接。
# Windows: 下方路径变量均可用 make VAR=... 覆盖以适配本机 oneAPI/MSVC 版本。

IFX    ?= ifx
FFLAGS ?= -r8 -double-size=64 -i8 -O2 -w
EXE    := venus_test.e

ifeq ($(OS),Windows_NT)
INTEL_BASE ?= C:/Program Files (x86)/Intel/oneAPI
MSVC_BASE  ?= C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207
WINKIT_LIB ?= C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0
# Git Bash 的 /usr/bin/link 会遮蔽 MSVC link.exe；-Qlocation,link 直接指定链接器目录。
# 同时显式导出 LIB，使裸 Git Bash（未跑 setvars.bat）下也能链接系统库。
export PATH := $(INTEL_BASE)/compiler/2025.2/bin:$(PATH)
export LIB  := $(MSVC_BASE)/lib/x64;$(WINKIT_LIB)/um/x64;$(WINKIT_LIB)/ucrt/x64;$(INTEL_BASE)/compiler/2025.2/lib
LINKOPT := -Qlocation,link,"$(MSVC_BASE)/bin/Hostx64/x64"
MKLFLAGS ?= "$(INTEL_BASE)/mkl/2025.2/lib/mkl_intel_lp64.lib" \
            "$(INTEL_BASE)/mkl/2025.2/lib/mkl_sequential.lib" \
            "$(INTEL_BASE)/mkl/2025.2/lib/mkl_core.lib"
else
LINKOPT :=
MKLFLAGS ?= -qmkl=sequential
endif

# 源文件按模块依赖顺序排列。
SRC := \
  venus_params.f90 venus_data.f90 input_parser.f90 \
  param_mapping.f90 presets.f90 harmonic_sampling.f90 \
  src_TEST/test_potentials.f90 \
  DVDQ.f \
  venus_input.f90 \
  src_TEST/interface_TEST.f90 ENERGY.f \
  src_VENUS/ANGVEL.f src_VENUS/BAREXC.f90 src_VENUS/CENMAS.f90 \
  src_VENUS/EBOND.f src_VENUS/EIGN.f src_VENUS/EIGOUT.f src_VENUS/ENMODE.f \
  src_VENUS/VFDATE.f src_VENUS/FGMTRX.f src_VENUS/FINAL.f src_VENUS/FINLNJ.f \
  src_VENUS/FMTRX.f src_VENUS/GAMA.f src_VENUS/GFINAL.f src_VENUS/GLPAR.f \
  src_VENUS/HOMOQP.f90 src_VENUS/INITEBK.f90 src_VENUS/INITQP.f90 \
  src_VENUS/LMEXCT.f src_VENUS/LMODE.f src_VENUS/NMODE.f90 \
  src_VENUS/ORTHAN.f90 src_VENUS/POTEN.f src_VENUS/RAND0.f \
  src_VENUS/RAND1.f src_VENUS/RANDST.f src_VENUS/ROTATE.f90 \
  src_VENUS/ROTEN.f90 src_VENUS/ROTN.f90 src_VENUS/SELECT.f90 \
  src_VENUS/SURF.f src_VENUS/THRMAN.f VENUS.f90 \
  src_VENUS/WEBOND.f src_VENUS/WENMOD.f src_VENUS/WLBOND.f \
  src_VENUS/CPUSEC.f src_VENUS/RADAU.f src_VENUS/SYMPLE.f \
  src_VENUS/THERMO.f src_VENUS/THERMBATH.f src_VENUS/GASDEV.f \
  src_VENUS/POTENZ.f src_VENUS/VERLET.f90 src_VENUS/PRINFO.f \
  src_VENUS/MICROCI.f src_VENUS/QMMICRO.f src_VENUS/VOLPSCONE.f src_VENUS/DENQ.f \
  src_VENUS/JMAXCALC.f src_VENUS/PROBJ.f90 \
  src_VENUS/ROTATEX.f src_VENUS/ROTATEY.f src_VENUS/ROTATEZ.f \
  src_VENUS/ArbitraryAxisRotation.f src_VENUS/ROTATEJM.f src_VENUS/ROTATEJKM.f \
  src_VENUS/GLO.f90 \
  src_VENUS/TEST.f src_VENUS/GWRITE.f90 src_VENUS/HIST.f90

$(EXE): $(SRC)
	$(IFX) $(FFLAGS) $(LINKOPT) $(SRC) $(MKLFLAGS) -o $(EXE)

clean:
	rm -f $(EXE) *.mod *.obj *.pdb *.ilk \
	  src_VENUS/*.mod src_VENUS/*.obj src_TEST/*.mod src_TEST/*.obj

.PHONY: clean
