#!/usr/bin/env julia

"""
Julia workflow entry point: dose predictions
===========================================

This script is the Julia counterpart of `code/analysis-code/run-dose-predictions.R`.

Expected input
--------------
- `results/output/<model>-bestfit-sample.jld2`
  or fallback `results/output/<model>-bestfit-multistart.jld2`

Run command
-----------
`julia --project=code/julia code/julia/scripts/run-dose-predictions.jl`
or run with: include("code/julia/scripts/run-dose-predictions.jl")
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "handel-qsp.jl"))
import .HandelQSP

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choice = "model2"

settings = HandelQSP.default_dose_prediction_settings(model_choice)
HandelQSP.run_dose_prediction_workflow(settings)

