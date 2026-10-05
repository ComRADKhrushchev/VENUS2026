# VENUS2026 completion rewrite - src/ build entry (K3: wildcard pickup + prebuild registry + build/ isolation)
#
# Usage:
#   make              # build build/venus.e with gfortran (default)
#   make FC=ifx       # build with ifx (the project's conventional flags)
#   make clean        # remove the whole build/ tree
#
# Drag-in (E11): a new methods/ or systems/ reg_*.f90 is picked up by the wildcard on the
#   next make; ../scripts/prebuild.sh regenerates build/registry_gen.f90 (the aggregate
#   registration module - registration is a BUILD ARTIFACT, the driver carries zero
#   hand-written registration lines) and build/depend.mk (object-level use-dependencies;
#   depending on the defining unit's .o sidesteps make's two-product (.o+.mod) problem and
#   is parallel-safe). The Makefile never changes again.
#
# Layout: sources at src/*.f90, src/*/*.f90 and the repo-level systems/*/*.f90; ALL
#   artifacts (.o/.mod/registry_gen.f90/depend.mk/venus.e) land in build/ - never the
#   source tree (the K2 stale-mod shadowing lesson, E14 task-8 note).

# := not ?= - FC is a make BUILT-IN variable whose default (f77) would survive ?=
FC     := gfortran
# PER-COMPILER build trees: .mod and .o formats are mutually INCOMPATIBLE between
# compilers - sharing one build/ made a `make FC=ifx` followed by a plain `make` read
# ifx .mod files with the gfortran reader ("Unexpected EOF"). gfortran and ifx each
# get their own subtree; `make clean` removes both.
BUILD  := build/$(FC)
EXE    := $(BUILD)/venus.e
PREBUILD := ../scripts/prebuild.sh

# ---- sources (wildcard, E11) ----
# filter-out build/: the artifact directory sits inside the tree pattern (*/ *.f90 would
# pick up the generated registry_gen.f90 twice - a ghost object and a circular dep);
# build/ as a WHOLE is excluded - any source parked under build/ (probe programs of a
# concurrent session, legacy src/build content) never compiles into venus.e
SRC_TREE := $(filter-out build/%,$(wildcard *.f90 */*.f90))
SRC_SYS  := $(wildcard ../systems/*/*.f90)

OBJ_TREE := $(patsubst %.f90,$(BUILD)/%.o,$(SRC_TREE))
OBJ_SYS  := $(patsubst ../systems/%.f90,$(BUILD)/sys/%.o,$(SRC_SYS))
OBJ_GEN  := $(BUILD)/registry_gen.o
OBJS     := $(OBJ_TREE) $(OBJ_SYS) $(OBJ_GEN)

# ---- compiler fork ----
LDPOST :=
ifeq ($(FC),ifx)
IFX    ?= ifx
# Windows ifx takes the module dir as -module:DIR (colon); the space-separated
# -module DIR form is SILENTLY ignored and the .mod files land in the CWD (= src/,
# poisoning later gfortran runs) - the root cause behind every "Unexpected EOF" round
ifeq ($(OS),Windows_NT)
MODDIR := -module:$(BUILD)
else
MODDIR := -module $(BUILD)
endif
FFLAGS := -r8 -double-size=64 -i8 -O2 -w -c $(MODDIR)
LDPOST := $(LINKOPT)
else
FFLAGS := -O2 -c -J$(BUILD)
endif

# ---- platform section (the ifx-Windows link path; the P0 skeleton's -Qlocation form
#      never survived a real link - make mangles the nested quotes and ifx then finds
#      Git Bash's /usr/bin/link. The setvars-style fix: put the MSVC linker bin FIRST on
#      PATH so ifx resolves link.exe by search, and export LIB for the bare-Git-Bash case) ----
ifeq ($(OS),Windows_NT)
INTEL_BASE ?= C:/Program Files (x86)/Intel/oneAPI
MSVC_BASE  ?= C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207
WINKIT_LIB ?= C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0
export PATH := $(MSVC_BASE)/bin/Hostx64/x64:$(INTEL_BASE)/compiler/2025.2/bin:$(PATH)
export LIB  := $(MSVC_BASE)/lib/x64;$(WINKIT_LIB)/um/x64;$(WINKIT_LIB)/ucrt/x64;$(INTEL_BASE)/compiler/2025.2/lib
LINKOPT :=
endif

# ---- rules ----
# prebuild runs EVERY make (phony): make cannot sense a wildcard set SHRINKING (a
# deleted reg file leaves no stale timestamp to compare), so the aggregate and the
# dependency file are regenerated on every invocation; the script's sync-protect keeps
# the output timestamps stable when content is unchanged (no rebuild cascades)
all: prebuild $(EXE)
# NOTE: all stays FIRST - make's default goal is the first non-.PHONY target; defining
# prebuild above all would silently reduce a bare `make` to the prebuild step alone

.PHONY: prebuild
prebuild: | $(BUILD)
	bash $(PREBUILD) . $(BUILD)

$(BUILD)/%.o: %.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $< -o $@

$(BUILD)/sys/%.o: ../systems/%.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $< -o $@

# the generated unit gets an EXPLICIT rule - the wildcard pattern would resolve its
# prerequisite to a nonexistent cwd-level registry_gen.f90 (circularity). Its reg-object
# prerequisites are stated STATICALLY: depend.mk is one make-run stale whenever the reg
# set changes (a freshly dragged-in reg module must compile BEFORE registry_gen, and
# make evaluates this run's graph against last run's depend.mk - a static rule closes
# that window; the residual limitation "tree-structure changes need a second make" is
# confined to non-reg use edges, which the skeleton does not change)
SRC_REGM := $(wildcard methods/reg_*.f90)
SRC_REGS := $(wildcard ../systems/*/reg_*.f90)
OBJ_REG  := $(patsubst %.f90,$(BUILD)/%.o,$(SRC_REGM)) \
            $(patsubst ../systems/%.f90,$(BUILD)/sys/%.o,$(SRC_REGS))
$(BUILD)/registry_gen.o: $(BUILD)/registry_gen.f90 $(OBJ_REG) | $(BUILD)
	$(FC) $(FFLAGS) $< -o $@

$(EXE): $(OBJS)
	$(FC) $(OBJS) -o $(EXE) $(LDPOST)

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf build
	rm -f *.mod
# the second rm is DEFENSIVE: in gfortran's module search the CURRENT directory precedes
# the -J directory - a stale cwd-level .mod (left by an ad-hoc compile) shadows the
# build-dir one and kills downstream compiles with "Unexpected EOF" (the E14 task-8
# stale-mod lesson, recurring at K3)

.PHONY: all clean

# object-level use-dependencies. BOOTSTRAP: after `make clean` the file does not exist
# at parse time (-include would silently yield an unordered graph and the first compile
# order would be wrong) - a parse-time $(shell) generates it BEFORE the include on the
# first run only (when absent); the phony prebuild rule handles later reg-set changes
ifeq ($(wildcard $(BUILD)/depend.mk),)
_DUMMY := $(shell bash $(PREBUILD) . $(BUILD))
endif
-include $(BUILD)/depend.mk
