#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "analysis", "FitModel.jl"))
include(joinpath(@__DIR__, "..", "plotting", "DiagnosticPlots.jl"))
using .ModelFitting: FitResult
using .DiagnosticPlots
using Serialization

function load_bestfit(path)
    open(path, "r") do io
        Serialization.deserialize(io)
    end
end

bestfit_path = length(ARGS) >= 1 ? ARGS[1] : "results/output/bestfit.jls"
bestfit_list = load_bestfit(bestfit_path)

mkpath("results/figures")

for (idx, bestfit) in enumerate(bestfit_list)
    fig_grid = residual_grid_plot(bestfit)
    save(joinpath("results/figures", "residuals$(idx).png"), fig_grid)

    fig_combined = residual_combined_plot(bestfit)
    save(joinpath("results/figures", "residuals-combined$(idx).png"), fig_combined)
end
