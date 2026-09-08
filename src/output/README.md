# output/ - L4 recording and offline post-processing

- **Online (inside the evolution loop)**: recorder.f90 - trajectory serialization
  (rec_reg_col column registry +
  rec_select/rec_frame/rec_final). One of the four legal call classes inside the evolution loop
(facade/propagator/recorder/predicate).
  - **Offline (not part of the evolution loop)**: final_state.f90 (generalized COM relative
  energies / angular momenta / final-state quantities),
  outcome.f90 (outcome-classification execution - looks up the event%outcome_tbl data table),
density.f90
(Beyer-Swinehart density of states + eigenvalue output), plot_hist.f90 (five histogram kinds,
  module name hist).

Boundary constraint: the four offline files **never touch the physical state** (they do not
modify Q/P/s/F and do not call facade/propagator),
they only do representation -> representation transforms and reductions (recording -> statistics).

Reference material (old library, read-only): GWRITE.f90 / FINAL.f / GFINAL.f / the VENUS.f90
outcome blocks /
DENQ.f / EIGOUT.f / HIST.f90.
