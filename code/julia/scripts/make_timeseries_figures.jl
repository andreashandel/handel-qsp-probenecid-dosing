#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "analysis", "FitModel.jl"))
include(joinpath(@__DIR__, "..", "plotting", "TimeSeriesPlots.jl"))
include(joinpath(@__DIR__, "..", "analysis", "DosePredictions.jl"))
using .ModelFitting: FitResult
using .TimeSeriesPlots
using .DosePredictions
using Serialization

function load_serialized(path)
    open(path, "r") do io
        Serialization.deserialize(io)
    end
end

bestfit_path = "results/output/bestfit.jls"
dose_path = "results/output/dose-response-results.jls"
if length(ARGS) >= 1
    bestfit_path = ARGS[1]
end
if length(ARGS) >= 2
    dose_path = ARGS[2]
end

bestfit_list = load_serialized(bestfit_path)
sim_list = load_serialized(dose_path)

bestfit = bestfit_list[1]
fig = plot_timeseries(data = bestfit.fitdata, modelfit = bestfit.odeout_df,
    dose_levels = ["no drug", "10 mg/kg", "100 mg/kg"])
mkpath("results/figures")
save("results/figures/bestfit1.png", fig)

baseline_ts = filter(:Schedule => ==("s1"), sim_list[1][:timeseries_df])
fig_s1 = plot_timeseries(data = nothing, modelfit = baseline_ts,
    dose_levels = ["1 mg/kg", "10 mg/kg", "100 mg/kg"])
save("results/figures/timeseries-baseline.png", fig_s1)

for (sched, filename) in [("s2", "timeseries-d2tx.png"),
        ("s3", "timeseries-d3tx.png"),
        ("s4", "timeseries-dailytx.png"),
        ("s5", "timeseries-singletx.png")]
    ts = filter(:Schedule => ==(sched), sim_list[1][:timeseries_df])
    fig_sched = plot_timeseries(data = nothing, modelfit = ts,
        dose_levels = ["1 mg/kg", "10 mg/kg", "100 mg/kg"])
    save(joinpath("results/figures", filename), fig_sched)
end
