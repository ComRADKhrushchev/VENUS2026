# VENUS2026 completion rewrite - src/ build entry (P0 skeleton)
# A single ifx invocation compiles and links, producing venus.e.
# Usage:
#   make          # build venus.e
#   make clean    # remove the executable and build intermediates
# Linux: load the oneAPI environment first (source setvars.sh) so ifx is visible.
#        The new program has no MKL dependency; the old library's -qmkl link block is not carried over.
# Windows: the path variables below can all be overridden with make VAR=... to match local oneAPI/MSVC versions.
# (Platform block ported from venus2026/Makefile; only the three MKL library links were dropped.)

IFX    ?= ifx
FFLAGS ?= -r8 -double-size=64 -i8 -O2 -w
EXE    := venus.e

ifeq ($(OS),Windows_NT)
INTEL_BASE ?= C:/Program Files (x86)/Intel/oneAPI
MSVC_BASE  ?= C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207
WINKIT_LIB ?= C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0
# Git Bash's /usr/bin/link shadows MSVC link.exe; -Qlocation,link points ifx at the linker directory directly.
# LIB is also exported explicitly so linking against system libraries works in a bare Git Bash (without setvars.bat).
export PATH := $(INTEL_BASE)/compiler/2025.2/bin:$(PATH)
export LIB  := $(MSVC_BASE)/lib/x64;$(WINKIT_LIB)/um/x64;$(WINKIT_LIB)/ucrt/x64;$(INTEL_BASE)/compiler/2025.2/lib
LINKOPT := -Qlocation,link,"$(MSVC_BASE)/bin/Hostx64/x64"
else
LINKOPT :=
endif

# Source files are listed in module dependency order (L5->L1->L2->L3->L0).
# P0 skeleton: SRC is an empty list for now; phase K3 (independent build of the new tree)
# fills in the fixed skeleton list PLUS a wildcard pickup of systems/ (E11 "drag-in"
# ruling: each system blackbox sits in its own folder with a reg_<system>_all module;
# a prebuild step scans systems/*/reg_*.f90 and generates registry_gen.f90, so dragging
# a folder in is enough - zero edits to Makefile or driver; .mod dependencies let make
# sort the compile order). Until then `make venus.e` only checks that the Makefile parses.
SRC :=

$(EXE): $(SRC)
	$(IFX) $(FFLAGS) $(LINKOPT) $(SRC) -o $(EXE)

clean:
	rm -f $(EXE) *.mod *.obj *.pdb *.ilk

.PHONY: clean
