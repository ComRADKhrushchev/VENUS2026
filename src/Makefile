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
BUILD  := build
EXE    := $(BUILD)/venus.e
PREBUILD := ../scripts/prebuild.sh

# ---- sources (wildcard, E11) ----
# filter-out build/: the artifact directory sits inside the tree pattern (*/ *.f90 would
# pick up the generated registry_gen.f90 twice - a ghost object and a circular dep)
SRC_TREE := $(filter-out $(BUILD)/%,$(wildcard *.f90 */*.f90))
SRC_SYS  := $(wildcard ../systems/*/*.f90)

OBJ_TREE := $(patsubst %.f90,$(BUILD)/%.o,$(SRC_TREE))
OBJ_SYS  := $(patsubst ../systems/%.f90,$(BUILD)/sys/%.o,$(SRC_SYS))
OBJ_GEN  := $(BUILD)/registry_gen.o
OBJS     := $(OBJ_TREE) $(OBJ_SYS) $(OBJ_GEN)

# ---- compiler fork ----
LDPOST :=
ifeq ($(FC),ifx)
IFX    ?= ifx
FFLAGS := -r8 -double-size=64 -i8 -O2 -w -c -module $(BUILD) $(LINKOPT)
LDPOST := $(LINKOPT)
else
FFLAGS := -O2 -c -J$(BUILD)
endif

# ---- platform section (ported from the P0 skeleton; the ifx-Windows link path) ----
ifeq ($(OS),Windows_NT)
INTEL_BASE ?= C:/Program Files (x86)/Intel/oneAPI
MSVC_BASE  ?= C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.44.35207
WINKIT_LIB ?= C:/Program Files (x86)/Windows Kits/10/Lib/10.0.26100.0
# Git Bash's /usr/bin/link shadows MSVC link.exe; -Qlocation,link points ifx at the linker directory directly.
# LIB is also exported explicitly so linking against system libraries works in a bare Git Bash (without setvars.bat).
export PATH := $(INTEL_BASE)/compiler/2025.2/bin:$(PATH)
export LIB  := $(MSVC_BASE)/lib/x64;$(WINKIT_LIB)/um/x64;$(WINKIT_LIB)/ucrt/x64;$(INTEL_BASE)/compiler/2025.2/lib
LINKOPT := -Qlocation,link,"$(MSVC_BASE)/bin/Hostx64/x64"
ifeq ($(FC),ifx)
FFLAGS := -r8 -double-size=64 -i8 -O2 -w -c -module $(BUILD) $(LINKOPT)
endif
endif

# ---- rules ----
all: $(EXE)

# grouped target (make >= 4.3): one invocation emits both generated artifacts
$(BUILD)/registry_gen.f90 $(BUILD)/depend.mk &: $(PREBUILD) $(SRC_TREE) $(SRC_SYS) | $(BUILD)
	bash $(PREBUILD) . $(BUILD)

$(BUILD)/%.o: %.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $< -o $@

$(BUILD)/sys/%.o: ../systems/%.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) $< -o $@

# the generated unit gets an EXPLICIT rule - the wildcard pattern would resolve its
# prerequisite to a nonexistent cwd-level registry_gen.f90 (circularity)
$(BUILD)/registry_gen.o: $(BUILD)/registry_gen.f90 | $(BUILD)
	$(FC) $(FFLAGS) $< -o $@

$(EXE): $(OBJS)
	$(FC) $(OBJS) -o $(EXE) $(LDPOST)

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)

.PHONY: all clean

# object-level use-dependencies (regenerated with the sources; safe to be absent)
-include $(BUILD)/depend.mk
