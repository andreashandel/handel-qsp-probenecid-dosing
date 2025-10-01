# Julia workflow quick reference

This document summarises the end-to-end execution order for reproducing the complete analysis pipeline in Julia.

1. **Model fitting** – estimates dynamic parameters by minimising the squared error objective.
   ```bash
   julia --project=code/julia code/julia/scripts/run_fit_model.jl
   ```
2. **Dose exploration** – reuses the serialized fit results to run all dosing regimens.
   ```bash
   julia --project=code/julia code/julia/scripts/run_dose_predictions.jl
   ```
3. **Best-fit visualisations** – generates the baseline fit plot and residual diagnostics.
   ```bash
   julia --project=code/julia code/julia/scripts/make_fit_figures.jl
   ```
4. **Dose response figures** – creates the panels comparing alternative schedules.
   ```bash
   julia --project=code/julia code/julia/scripts/make_dose_response_figures.jl
   ```
5. **Time-series figures** – renders supplemental trajectory plots for each schedule.
   ```bash
   julia --project=code/julia code/julia/scripts/make_timeseries_figures.jl
   ```
6. **Tables** – exports the parameter summary tables.
   ```bash
   julia --project=code/julia code/julia/scripts/make_tables.jl
   ```

All scripts accept the `--help` flag (for example `julia ... run_fit_model.jl --help`) to display available keyword overrides.
