# Julia Replication of QSP Analyses

This directory contains a self-contained Julia implementation of the quantitative systems pharmacology (QSP) analyses and figure generation that were originally implemented in R.  The Julia code mirrors the structure of the existing workflow and can be used as a drop-in replacement for running model fitting, dose exploration, figure generation and table creation.

## Contents

- `Project.toml` – Julia project file listing dependencies.
- `analysis/` – model definition, parameter fitting and dose prediction simulations.
- `plotting/` – figure-generation utilities mirroring the original R graphics.
- `scripts/` – ready-to-run command-line entry points orchestrating the workflow.
- `docs/` – additional task-focused documentation.

## Getting started

1. Install Julia (version 1.9 or later is recommended).
2. From the repository root, activate the Julia project:
   ```julia
   import Pkg
   Pkg.activate("code/julia")
   Pkg.instantiate()
   ```
3. Run the scripted workflows as described in the readme files within `analysis/`, `plotting/`, and `scripts/`.

The Julia implementation reads the same input files as the R code (`data/processed-data/*`) and writes outputs to the same locations in `results/`.  All file paths are resolved relative to the repository root.
