"""
`HandelQSP`
============

This module is a Julia re-implementation of the core analysis workflow currently
implemented in R for the probenecid QSP project. The design goal is functional
parity with three R entry points:

- `run-fit.R`
- `run-dose-predictions.R`
- `run-plotting.R`

The module is intentionally split into multiple files so that users unfamiliar
with Julia can navigate one concern at a time (data loading, simulation,
objective evaluation, fitting, predictions, plotting).

Usage pattern
-------------
The scripts in `code/julia/scripts/` load this module and call high-level helper
functions. You can also work interactively from the Julia REPL by including
this file and then calling exported functions.

```julia
include("code/julia/src/handel-qsp.jl")
using .HandelQSP
cfg = default_fit_settings("model1")
```
"""
module HandelQSP

using CSV
using DataFrames
using Dates
using DifferentialEquations
using BlackBoxOptim
using NLopt
using QuasiMonteCarlo
using Statistics
using Random
using Printf
using LinearAlgebra
using Distributions
using JLD2
using CairoMakie

include("paths-function.jl")
include("virus-transform-function.jl")
include("data-and-config-function.jl")
include("simulators-function.jl")
include("objective-function.jl")
include("fit-workflow-function.jl")
include("dose-predictions-function.jl")
include("plot-workflow-function.jl")

export project_root
export virus_quantity_name, transform_virus, inverse_transform_virus
export load_fit_data, load_fixed_parameters, compute_sigma_settings, build_model_config
export default_fit_settings, run_fit_workflow
export default_dose_prediction_settings, run_dose_prediction_workflow
export default_plot_settings, run_plot_workflow
export save_julia_object, load_julia_object

end
