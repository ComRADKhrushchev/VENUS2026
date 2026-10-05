# benchmark_oldvenus - clean VENUS05-lineage tree (gfortran)
# Machinery: src/*.f  (stock files, untouched, from VENUS_TSH_HO2/src_VENUS)
# System side: sys/  (PES seam, engine, channel test, writer, stubs, verlet)
# Build from this directory:  make          -> venus.e
#                             make fdcheck  -> fdcheck.exe (seam force audit)
FC = gfortran
FFLAGS = -O2 -std=legacy -w -I.

MACH := $(wildcard src/*.f)
SYSF := sys/TEST_h2ag.f sys/GWRITE_h2ag.f sys/h2ag_engine.f
SYS90 := sys/PES_h2ag.f90 sys/friction_stub.f90 sys/VERLET_czz.f90

OBJS = $(MACH:.f=.o) $(SYSF:.f=.o) $(SYS90:.f90=.o)

venus.e: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@

fdcheck: sys/fdcheck.o sys/PES_bkmp2.o sys/bkmp2_engine.o
	$(FC) $(FFLAGS) $^ -o fdcheck.exe

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f src/*.o sys/*.o venus.e fdcheck.exe

.PHONY: clean
