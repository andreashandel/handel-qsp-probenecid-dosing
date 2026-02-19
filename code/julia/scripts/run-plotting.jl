#!/usr/bin/env julia

"""
Julia workflow entry point: plotting and tables
===============================================

This script is the Julia counterpart of `code/plotting-code/run-plotting.R`.

Run command
-----------
`julia --project=code/julia code/julia/scripts/run-plotting.jl`
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "handel-qsp.jl"))
import .HandelQSP

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choice = "model1"
nsamp = 1

settings = merge(HandelQSP.default_plot_settings(model_choice), (nsamp = nsamp,))
HandelQSP.run_plot_workflow(settings)
