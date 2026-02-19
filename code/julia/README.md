# Julia Workflow for QSP Probenecid Project

This folder contains a Julia re-implementation of the current **analysis workflow** in R for:

- model fitting (`run-fit.R` equivalent)
- dose-response predictions (`run-dose-predictions.R` equivalent)
- plotting and parameter tables (`run-plotting.R` equivalent)

This Julia workflow does **not** re-implement:

- anything in `processing-code`
- the Shiny app (`app.R`)

## Why this Julia version exists

The R fitting workflow is effective but uses local NLopt methods with multistart. In broad and rugged objective landscapes, local methods can spend too much compute in locally good basins. The Julia implementation uses a **hybrid global-local strategy**:

1. **Global exploration** with Differential Evolution (`BlackBoxOptim`)
2. **Local polishing** under bounds with `Optim` (`Fminbox + NelderMead`)
3. Optional **local restarts** around the global best point

This typically gives stronger global coverage than local-only multistart while retaining a robust final refinement.

Sampling-stage behavior:
- Stage 1 (baseline) uses global DE + local refinement.
- Stage 2 (fixed-parameter samples) uses local-only refinement, seeded from the stage-1 best-fit parameters.

## Folder structure

- `Project.toml`: Julia dependencies
- `src/handel-qsp.jl`: top-level module file that includes all subfiles
- `src/*.jl`: workflow implementation modules
- `scripts/run-fit.jl`: fitting entry point
- `scripts/run-dose-predictions.jl`: prediction entry point
- `scripts/run-plotting.jl`: plotting entry point
- `docs/optimization-strategy.md`: deeper explanation of optimizer design and trade-offs

## Installation and first run

From repository root:

```bash
julia --project=code/julia -e 'import Pkg; Pkg.instantiate()'
```

Then run scripts:

```bash
julia --project=code/julia code/julia/scripts/run-fit.jl
julia --project=code/julia code/julia/scripts/run-dose-predictions.jl
julia --project=code/julia code/julia/scripts/run-plotting.jl
```

## Re-running after code changes (without restarting Julia)

If you work in the Julia REPL, use:

```julia
using Revise
include("code/julia/scripts/run-fit.jl")
```

After edits, run `include("code/julia/scripts/run-fit.jl")` again.
The scripts now call functions with `HandelQSP.<name>` qualification, which avoids
name ambiguity issues in `Main` during repeated includes.

## Main outputs

The Julia scripts write files to the same high-level output locations already used by the project.

### Fitting

- `results/output/<model>-bestfit-multistart.jld2`
- `results/output/<model>-bestfit-sample.jld2` (when sampling stage is enabled)

Each saved bestfit object contains only fields needed downstream:

- `fitpars`
- `fitparnames`
- `fixedpars`
- `Y0`
- `fitdata`
- `parlabels`
- `objective`

### Dose predictions

- `results/output/<model>-dose-response-results.jld2`
- `results/output/<model>-dose-response-messages.csv`

### Plotting

- `results/figures/*-julia.png`
- `results/tables/*-julia.csv`

Time-series plots follow the R layout with 9 state panels:
- `U`, `I`, `V`
- `F`, `A`, `S`
- `Ad`, `Ac`, `At`

## Configuration model

Each script keeps a **single user settings block** near the top. For fitting, `scripts/run-fit.jl` now exposes all main knobs (optimizer budgets, solver settings, sigma handling, sampling-stage settings, output path, and worker count) in one place, similar to `run-fit.R`.

### Key fit settings to tune first

- `model_choice`
- `run_sampling_stage`
- `global_maxeval`
- `global_trace_mode`
- `verbose_fit_log`
- `reuse_previous_bestfit` / `previous_bestfit_file` / `previous_bestfit_n`
- `local_maxiters`
- `n_local_restarts`
- `sampling_local_optimizer`
- `sampling_local_maxiters`
- `solver_settings` (default uses stiff solver `rodas5p`)
- `verbose_fit_log` and `log_top_n` for runtime diagnostics

Parallel note:
- `n_workers` is honored for threaded sampling-stage fitting (capped by available Julia threads).
- Julia thread count is set at launch; start Julia with e.g. `JULIA_NUM_THREADS=8`.

If runtime is too long, reduce `global_maxeval` first.
If solutions look under-explored, increase `global_maxeval` and possibly `global_population` in `default_fit_settings`.

### Minimal smoke-test settings

For a fast end-to-end test of fitting, use the `quick_smoke_test` switch in `scripts/run-fit.jl`:

- set `quick_smoke_test = true`
- keep `run_sampling_stage = false`

This profile uses:

- `global_maxeval = 1500`
- `global_population = 30`
- `local_maxiters = 200`
- `n_local_restarts = 3`
- `solver_settings = (solvertype = "rodas5p", tols = 1e-8, tfinal = 7.0, dt = 0.2)`

These settings are intended for workflow validation only, not final inference.

## Notes for users new to Julia

- Julia compiles code on first run; first execution is slower than subsequent runs.
- Most workflow functions return plain containers (`NamedTuple`, `Dict`, `DataFrame`) to keep debugging straightforward.
- The project uses `.jld2` files instead of `.Rds` because they preserve nested Julia objects naturally.

## Mapping to R files

- `scripts/run-fit.jl` -> `code/analysis-code/run-fit.R`
- `scripts/run-dose-predictions.jl` -> `code/analysis-code/run-dose-predictions.R`
- `scripts/run-plotting.jl` -> `code/plotting-code/run-plotting.R`

## Known differences from R implementation

- Solver backend is DifferentialEquations.jl (not deSolve), but equations and dosing logic are matched.
- Fitting uses hybrid DE + local refinement instead of the original staged NLopt multistart.
- Output format is `.jld2` instead of `.Rds`.

These differences are intentional and are aimed at more robust global search behavior.
