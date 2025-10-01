#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "analysis", "FitModel.jl"))
include(joinpath(@__DIR__, "..", "plotting", "TableExporter.jl"))
using .ModelFitting: FitResult
using .TableExporter
using Serialization

bestfit_path = length(ARGS) >= 1 ? ARGS[1] : "results/output/bestfit.jls"
output_dir = length(ARGS) >= 2 ? ARGS[2] : "results/tables"

open(bestfit_path, "r") do io
    bestfit_list = Serialization.deserialize(io)
    export_parameter_tables(bestfit_list; output_dir = output_dir)
end
