#!/usr/bin/env julia

"""
Julia workflow entry point: structural identifiability
======================================================

This script evaluates structural identifiability for the parameters currently
fitted in model1 and model2 using StructuralIdentifiability.jl.

Run command
-----------
`julia --project=code/julia code/julia/scripts/run-structural-identifiability.jl`
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "handel-qsp.jl"))
import .HandelQSP

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choices = ["model1", "model2"]
prob_threshold = 0.99

settings = merge(
    HandelQSP.default_structural_identifiability_settings(),
    (
        model_choices = model_choices,
        prob_threshold = prob_threshold,
    ),
)

HandelQSP.run_structural_identifiability_workflow(settings)
