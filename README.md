# VENUS2026 Project

An independent working project for streamlining and refactoring the VENUS program
(classical adiabatic QCT; the nonadiabatic electronic axis is contract-first scaffolding —
the first container-composed method package lands per the E27 ruling: the physical-system
container is the composition point of system-specific code).

## Layout

```
src/           the NEW PROGRAM (completion rewrite, L0-L5 layering; sole execution
               authority docs/final_completion_plan.md). Build: cd src && make [FC=ifx];
               guidance src/CLAUDE.md, extension points src/EXTENDING.md, layer map
               src/README.md. Verification: cases/ (five-case formal set, per-case
               results.txt) + scripts/check_* criteria programs.
systems/       physical-system container material (E27: the container is the composition
               point — H definition + system-specific method packages live here, one folder
               per system; reference/ = the pure-math harmonic-well container driving the
               case set)
venus2026/     the OLD minimal form — READ-ONLY reference (provenance of algorithms
               and formulas); neither migrated nor modified (its CLAUDE.md describes
               the reference tree itself)
docs/          design documents: design_spec_completeness.md (the specification, §4.2
               backfilled with the landed interfaces), final_completion_plan.md
               (execution authority), optimal_structure.md (implemented),
               completeness_audit_src.{md,json} (S-phase new-tree audit),
               completeness_audit.{md,json} (old-tree audit, historical),
               plans/ (per-phase records P0..V, S incl. the E27 ruling)
cases/         the five-case verification set (harmonic_scatter / radau_conservation /
               surface_window / b_layer / chain_frequencies)
scripts/       check_* criteria programs + probe drivers + prebuild.sh
```

## Status

Phases P0 → K1-K3 → C1-C8 → V complete; S (closing: re-audit + document wrap-up)
executed 2026-09-16 with the mid-phase E27 architecture ruling (the physical-system
container as composition point; the ws interchange is retired — historical term,
docs/plans/S.md E27). The five-case set is green and reproducible
(same-ISEED byte-identical); the old chain cross-check holds bitwise
(max|dF| = max|dE| = 0, dual compiler).
自定义 PES 接入新增 Tier 1 快通道——单文件魔数名容器零 reg 接入
（docs/plans/2026-10-12-auto-pes-container.md）。

## Notes

- Build and run conventions: src/CLAUDE.md (single-entry Makefile, wildcard drag-in
  + prebuild registration, per-compiler build trees, cases run from their folders).
- The original repository carries .git history for the legacy tree: it remains the
  provenance record; the workspace root is the active repository.
