# Analysis workflow (Julia)

The files in this folder reproduce the ODE model simulation, parameter fitting and dose exploration steps that were previously implemented in `code/analysis-code` (R).  All functions are documented and can be imported from other Julia code, while the companion scripts in `code/julia/scripts` provide ready-to-run entry points.

## Files

- `SimulateModel.jl` – defines the pharmacodynamic/pharmacokinetic model and exposes [`simulate_model`] for integrating the ODEs with dosing callbacks.
- `FitModel.jl` – implements data ingestion, objective construction and NLopt-based calibration.
- `DosePredictions.jl` – runs counterfactual dosing scenarios for calibrated parameter samples.

## Running the analysis

With the Julia project activated (see `code/julia/README.md`), you can run the full calibration pipeline from the command line:

```bash
julia --project=code/julia code/julia/scripts/run_fit_model.jl
```

This command writes a serialized result object (`results/output/bestfit.jls`) mirroring the structure of the original RDS file.  The subsequent dose exploration can then be launched with:

```bash
julia --project=code/julia code/julia/scripts/run_dose_predictions.jl
```

Both scripts expose keyword arguments (documented inside the files) for adjusting solver tolerances, Latin hypercube sampling of fixed parameters, and the choice of NLopt algorithm.
