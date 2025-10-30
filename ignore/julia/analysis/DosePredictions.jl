module DosePredictions

using DataFrames
using Serialization

include("SimulateModel.jl")
include("FitModel.jl")
using .QSPModel: simulate_model
using .ModelFitting: FitResult, namedtuple_from

export simulate_dose_predictions, run_dose_predictions, load_fit_results

trapz(x::AbstractVector, y::AbstractVector) = sum(diff(x) .* (y[1:end-1] .+ y[2:end]) .* 0.5)

"""
    load_fit_results(path) -> Vector{FitResult}

Convenience loader for the serialized fit results produced by
[`ModelFitting.fit_model`].
"""
function load_fit_results(path::AbstractString)
    open(path, "r") do io
        return Serialization.deserialize(io)
    end
end

function compute_percent_reduction(df::DataFrame)
    grouped = groupby(df, :Schedule)
    out = DataFrame()
    for sub in grouped
        base_v = first(sub.AUCV)
        base_f = first(sub.AUCF)
        base_s = first(sub.AUCS)
        sub.perc_AUCV = (base_v .- sub.AUCV) ./ base_v .* 100
        sub.perc_AUCF = (base_f .- sub.AUCF) ./ base_f .* 100
        sub.perc_AUCS = (base_s .- sub.AUCS) ./ base_s .* 100
        append!(out, sub)
    end
    return out
end

"""
    simulate_dose_predictions(bestfit; tfinal, dt, ts_doses)

Mirror of `code/analysis-code/dose-predictions-simulator-function.R`.  Returns
a dictionary containing the summary area-under-the-curve metrics, percentage
reductions relative to the no-drug baseline and the stored time-series for
selected doses.
"""
function simulate_dose_predictions(bestfit::FitResult;
    tfinal::Float64 = 7.0, dt::Float64 = 0.01,
    ts_doses::Vector{Float64} = [1.0, 10.0, 100.0])

    params_nt = NamedTuple{Tuple(bestfit.fitparnames)}(Tuple(bestfit.solution))
    fixed_names = collect(keys(bestfit.fixedpars))
    fixed_values = [bestfit.fixedpars[name] for name in fixed_names]
    fixed_nt = namedtuple_from(fixed_names, fixed_values)
    Y0 = bestfit.Y0

    base_doses = sort(unique(vcat(0.0, ts_doses, 10 .^ range(-3, 3; length = 20))))

    function simulate_dose_response(doses::Vector{Float64}, txstart::Float64,
        txend::Float64, txinterval::Float64, schedule_name::String)
        summary_df = DataFrame(Dose = doses, AUCV = fill(NaN, length(doses)),
            AUCF = fill(NaN, length(doses)), AUCS = fill(NaN, length(doses)),
            Schedule = fill(schedule_name, length(doses)))
        ts_store = DataFrame[]

        for (idx, dose) in enumerate(doses)
            kwargs = merge(Y0, params_nt, fixed_nt,
                (Ad0 = dose / 50.0, txstart = txstart, txinterval = txinterval,
                    txend = txend, tstart = 0.0, tfinal = tfinal, dt = dt))
            odeout = simulate_model(; kwargs...)

            summary_df.AUCV[idx] = trapz(odeout.time, log10.(max.(1.0, odeout.V)))
            summary_df.AUCF[idx] = trapz(odeout.time, odeout.F)
            summary_df.AUCS[idx] = trapz(odeout.time, odeout.S)

            if dose in ts_doses
                tmp = copy(odeout)
                tmp.Dose = fill(dose, nrow(tmp))
                tmp.Schedule = fill(schedule_name, nrow(tmp))
                push!(ts_store, tmp)
            end
        end

        timeseries_df = isempty(ts_store) ? DataFrame() : reduce(vcat, ts_store)
        return (summary = summary_df, timeseries = timeseries_df)
    end

    schedule_defs = [
        (txstart = 1.0, txend = 4.0, txinterval = 0.5, name = "s1"),
        (txstart = 2.0, txend = 5.0, txinterval = 0.5, name = "s2"),
        (txstart = 3.0, txend = 6.0, txinterval = 0.5, name = "s3"),
        (txstart = 1.0, txend = 4.0, txinterval = 1.0, name = "s4"),
        (txstart = 1.0, txend = 1.0, txinterval = 1.0, name = "s5"),
    ]

    results = map(schedule_defs) do sched
        simulate_dose_response(base_doses, sched.txstart, sched.txend,
            sched.txinterval, sched.name)
    end

    all_results_df = reduce(vcat, getindex.(results, :summary))
    timeseries_df = reduce(vcat, getindex.(results, :timeseries))

    label_map = Dict(
        "s1" => "baseline",
        "s2" => "d2 start",
        "s3" => "d3 start",
        "s4" => "daily tx",
        "s5" => "single tx",
    )

    reduction_df = compute_percent_reduction(all_results_df)
    reduction_df.Scenario = [label_map[s] for s in reduction_df.Schedule]

    return Dict(
        :all_results_df => all_results_df,
        :reduction_df => reduction_df,
        :timeseries_df => timeseries_df,
    )
end

"""
    run_dose_predictions(bestfit_list; result_path)

Iterate `simulate_dose_predictions` across all `FitResult` objects.  When
`result_path` is non-empty the returned list is serialized to disk.
"""
function run_dose_predictions(bestfit_list::Vector{FitResult};
    result_path::AbstractString = "results/output/dose-response-results.jls")
    sims = [simulate_dose_predictions(fit) for fit in bestfit_list]
    if !isempty(result_path)
        mkpath(dirname(result_path))
        open(result_path, "w") do io
            Serialization.serialize(io, sims)
        end
    end
    return sims
end

end # module
