# VENUS2026 Project

An independent working project for streamlining and refactoring the VENUS program.

## Layout

```
venus2026/    program sources (moved out of C@Au_TDHF/venus2026_stage; 223 git-managed
              files, work products fort.*/stdout* not migrated; the original repository
              carries .git history, remote github.com/ComRADKhrushchev/VENUS2026.git)
docs/         design documents
  design_spec_completeness.md   program-completeness design spec (derivation document)
```

## Current-phase task

Streamline VENUS2026: cut redundant functionality. The architectural adjudication basis is
docs/design_spec_completeness.md — the decomposition into four classes of pure subroutines plus the
control kernel; the irreducible interfaces A (the force-contribution facade, internally a dual
registration axis) and B′ (the initial-state distribution family); and the closed selection families
(nuclear propagators, the termination-predicate grammar, outcome-classification data).

## Notes

- venus2026/ currently has no .git (only the managed files were copied during migration); when
  history is needed, go back to C@Au_TDHF/venus2026_stage, or continue working with git in the
  original repository and then sync
- Build and run conventions are in venus2026/CLAUDE.md (Makefile single entry, ifx, internal units
  Å/amu/10 fs, KEYWORD=VALUE input)
