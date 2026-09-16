# AGENTS.md — VENUS2026 workspace guide

VENUS2026 is a completion-oriented Fortran rewrite of the VENUS chemical dynamics program
(adiabatic QCT). The active codebase is `src/`; `venus2026/` is the OLD minimal form kept as a
read-only reference — never edit or migrate it. Repo docs and commit messages are largely
Chinese; match the surrounding language when writing docs/comments.

## Authority documents (read before touching sensitive areas)

- `docs/final_completion_plan.md` — sole execution authority (target end state, sequence, acceptance).
- `docs/design_spec_completeness.md` — architecture adjudication basis.
- `docs/optimal_structure.md` — target structure and material mapping.
- `src/README.md` — layer table (L0–L5), layout, and the assembly call-chain statute.
- `docs/plans/` — per-phase task records (C1–C8, V, K1–K3...); history lives here, not in code comments.

## Build & run

All builds go through `src/Makefile` (single entry):

```
cd src
make              # gfortran (default) -> build/gfortran/venus.e
make FC=ifx       # ifx (project's conventional flags: -r8 -double-size=64 -i8 -O2 -w)
make clean        # removes the whole build/ tree (both compiler subtrees)
```

- Registration is a BUILD ARTIFACT: `scripts/prebuild.sh` regenerates `build/<FC>/registry_gen.f90`
  and `depend.mk` on every make. A new method = one file in `src/methods/` (or `systems/<name>/`)
  plus a thin `reg_*.f90`; the wildcard picks it up — never hand-edit the Makefile or registry.
- Artifacts land only in `build/$(FC)/` (per-compiler trees; .mod formats are incompatible between
  compilers). Stale `src/*.mod` or `*.obj` in the source tree cause "Unexpected EOF" — never compile
  ad hoc into `src/`.
- Windows ifx quirks are already handled in the Makefile (MSVC link.exe path shadowing, `-module:`).
- Run: `cd cases/<name> && ../../src/build/<FC>/venus.e` — input is hardcoded `input_qct.txt`
  in the run directory; outputs are `fort.*` / `stdout*` (gitignored).

## Architecture rules that matter for edits

- Layers L0–L5 depend only downward (see `src/README.md`). The control kernel (`driver.f90`)
  sees physics only through `facade/` (E14 four-entry contract: force_eval / elec_prop /
  mqc_statistic; termination via the container_term slot).
- **E27 composition root (2026-09-16)**: the container (`systems/<sys>/`) delivers the
  composed active force kernel (`container_pes_i(q, mass, g, v)` in eV, fully assigned);
  the former ws interchange and the h_form/h_deliver matchmaking are retired. System-specific
  electronic method packages live in the system folder and register from there
  (`method_reg(name, prop, match [, sync])`); statistics write-back via `occ_set`/`n_hop_set`.
- Registry-family members (methods/, systems/) never depend on each other; they may use utils/.
- Unit conversion happens only at the force_facade boundary; internal units Å/amu/10 fs.
- One single random stream for the whole library (`utils/rng.f90`) — never add another RNG.
- "Legislated constants" (rounding conventions, e.g. conversion factors) are fixed by decision
  recorded in docs/plans; use them verbatim in analytic checks, not CODATA re-derivations.
- Input grammar is KEYWORD=VALUE (`src/input.f90`); `GWRITE_LEVEL` bare default is 2.

## Conventions & verification

- File headers: one-line purpose + short Design note (≤12 lines); subroutine headers 1–2 lines
  wrapped in `---`; numbered `! 1.` `! 2.` blocks for main logic. No historical E-number
  adjudication narratives, old-repo line references, or task-phase tags in source comments —
  those live in docs/plans.
- `scripts/check_*.f90` are criterion/verification drivers (gfortran, self-contained). Verify
  changes with the relevant check program plus a clean dual build (`make` and `make FC=ifx`).
- Commit messages in this repo are detailed Chinese narratives describing evidence and
  verification legs; follow that style.
