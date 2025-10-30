#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "plotting", "DoseResponsePlots.jl"))
using .DoseResponsePlots
using Serialization

function load_results(path)
    open(path, "r") do io
        Serialization.deserialize(io)
    end
end

results_path = "results/output/dose-response-results.jls"
if length(ARGS) >= 1
    results_path = ARGS[1]
end

sim_list = load_results(results_path)

save_dose_response_figures(sim_list, "results/figures/dose-response-baseline.png";
    scenarios = ["baseline"])
save_dose_response_figures(sim_list, "results/figures/dose-response-txstart.png";
    scenarios = ["baseline", "d2 start", "d3 start"])
save_dose_response_figures(sim_list, "results/figures/dose-response-txinterval.png";
    scenarios = ["baseline", "daily tx", "single tx"])
