module ModelFitting

using CSV
using CategoricalArrays
using DataFrames
using Dates
using Statistics
using Optimization
using OptimizationNLopt
using QuasiMonteCarlo
using Random
using Serialization

include("SimulateModel.jl")
using .QSPModel: simulate_model

export fit_model, FitResult, load_fit_inputs, namedtuple_from

const DEFAULT_SCENARIO_ORDER = [
    "NoTreatment",
    "PanCytoVir10mg",
    "PanCytoVir100mg",
]

struct FitResult
    solution::Vector{Float64}
    objective::Float64
    x0::Vector{Float64}
    fitparnames::Vector{Symbol}
    fixedpars::Dict{Symbol, Float64}
    Y0::NamedTuple
    fitdata::DataFrame
    parlabels::Vector{String}
    simresult::DataFrame
    odeout_df::DataFrame
    algorithm::Symbol
end

"""
    load_fit_inputs(; data_path, fixed_path) -> (fitdata, fixedpars_df, dose_map)

Load the processed experimental data and fixed parameter definitions used by
the QSP calibration workflow.  The returned data frame mirrors the structure
of the original R implementation and includes additional helper columns used by
Julia scripts.
"""
function load_fit_inputs(; data_path::AbstractString,
    fixed_path::AbstractString)
    fitdata = CSV.read(data_path, DataFrame)
    fitdata.Scenario = categorical(fitdata.Scenario,
        ordered = true, levels = DEFAULT_SCENARIO_ORDER)
    dose_map = Dict(
        DEFAULT_SCENARIO_ORDER[1] => 0.0,
        DEFAULT_SCENARIO_ORDER[2] => 10.0,
        DEFAULT_SCENARIO_ORDER[3] => 100.0,
    )
    fitdata.Dose = [dose_map[string(s)] for s in fitdata.Scenario]
    fitdata.xvals = fitdata.Day

    fixedpars_df = CSV.read(fixed_path, DataFrame)
    return fitdata, fixedpars_df, dose_map
end

function _clamp!(values::Vector{Float64}, lb::Vector{Float64},
    ub::Vector{Float64})
    for i in eachindex(values)
        values[i] = clamp(values[i], lb[i], ub[i])
    end
    return values
end

function _nlopt_solver(alg::Symbol)
    mapping = Dict(
        :LN_COBYLA => OptimizationNLopt.LN_COBYLA(),
        :LN_NELDERMEAD => OptimizationNLopt.LN_NELDERMEAD(),
        :LN_SBPLX => OptimizationNLopt.LN_SBPLX(),
    )
    return get(mapping, alg) do
        error("Unsupported NLopt algorithm: $(alg)")
    end
end

namedtuple_from(names::Vector{Symbol}, values::Vector{Float64}) =
    NamedTuple{Tuple(names)}(Tuple(values))

dict_from(names::Vector{Symbol}, values::Vector{Float64}) =
    Dict(zip(names, values))

function _align_predictions(model_df::DataFrame, data_df::DataFrame,
    variable::Symbol; logvirus::Bool = false)
    if nrow(data_df) == 0
        return DataFrame()
    end
    subset = select(model_df, :time, variable)
    rename!(subset, :time => :xvals, variable => :Predicted)
    if logvirus
        subset.Predicted .= log10.(max.(1.0, subset.Predicted))
    end
    merged = leftjoin(data_df, subset, on = :xvals)
    return merged
end

function _objective_factory(fitdata::DataFrame, Y0::NamedTuple,
    fixed_names::Vector{Symbol}, fixed_values::Vector{Float64},
    param_names::Vector{Symbol}, tfinal::Float64, dt::Float64,
    doses::Vector{Float64}, scenarios::Vector{String},
    Vmax::Float64, Innmax::Float64, Symmax::Float64)

    fixed_nt = namedtuple_from(fixed_names, fixed_values)

    function objective(x)
        params_nt = NamedTuple{Tuple(param_names)}(Tuple(x))
        total = 0.0

        for (dose, scen) in zip(doses, scenarios)
            kwargs = merge(Y0, params_nt, fixed_nt,
                (Ad0 = dose, txstart = 1.0, txinterval = 0.5,
                    txend = 4.0, tstart = 0.0, tfinal = tfinal, dt = dt))
            odeout = try
                simulate_model(; kwargs...)
            catch err
                return 1.0e10
            end

            scen_mask = String.(fitdata.Scenario) .== scen
            scen_data = fitdata[scen_mask, :]

            Vvals = scen_data[scen_data.Quantity .== "LogVirusLoad", :]
            Innvals = scen_data[scen_data.Quantity .== "IL6", :]
            Symvals = scen_data[scen_data.Quantity .== "Weight", :]

            Vpred = _align_predictions(odeout, Vvals; variable = :V,
                logvirus = true)
            Innpred = _align_predictions(odeout, Innvals; variable = :F)
            Sympred = _align_predictions(odeout, Symvals; variable = :S)

            if (nrow(Vpred) == 0) || any(ismissing, Vpred.Predicted)
                return 1.0e10
            end
            if (nrow(Innpred) == 0) || any(ismissing, Innpred.Predicted)
                return 1.0e10
            end
            if (nrow(Sympred) == 0) || any(ismissing, Sympred.Predicted)
                return 1.0e10
            end

            Virtot = mean(((Vpred.Predicted .- Vpred.Value) ./ Vmax) .^ 2)
            Inntot = mean(((Innpred.Predicted .- Innpred.Value) ./ Innmax) .^ 2)
            Symtot = mean(((Sympred.Predicted .- Sympred.Value) ./ Symmax) .^ 2)
            total += Virtot + Inntot + Symtot
        end
        return total
    end

    return objective
end

"""
    fit_model(; kwargs...) -> Vector{FitResult}

Calibrate the QSP model using the same logic as `code/analysis-code/fit-model.R`.
The routine constructs the least-squares objective across all three experimental
scenarios, performs optional Latin hypercube sampling of fixed parameters and
stores rich diagnostic metadata for downstream plotting.

## Keyword arguments
- `data_path`, `fixed_path`: input CSV locations.
- `result_path`: output path for the serialized Julia results (`.jls`).
- `nsamp`: number of Latin hypercube samples for fixed parameters (0 → only
  baseline).
- `rng`: random number generator used for sampling.
- `algorithm`: NLopt algorithm symbol (`:LN_COBYLA`, `:LN_NELDERMEAD`,
  `:LN_SBPLX`).
- `maxeval`, `ftol_rel`: optimizer settings mirroring the R script.
- `tfinal`, `dt`: simulation horizon and output cadence for the objective.
- `save_results`: toggle writing to disk.
- `verbose`: whether to log progress to stdout.
"""
function fit_model(; data_path::AbstractString = "data/processed-data/processeddata.csv",
    fixed_path::AbstractString = "data/processed-data/fixed-parameters.csv",
    result_path::AbstractString = "results/output/bestfit.jls",
    nsamp::Int = 0, rng::AbstractRNG = Random.MersenneTwister(123),
    algorithm::Symbol = :LN_COBYLA, maxeval::Int = 1_000,
    ftol_rel::Float64 = 1.0e-10, tfinal::Float64 = 7.0,
    dt::Float64 = 0.02, save_results::Bool = true,
    verbose::Bool = true)

    fitdata, fixedpars_df, dose_map = load_fit_inputs(; data_path = data_path,
        fixed_path = fixed_path)

    fixed_names = Symbol.(strip.(String.(fixedpars_df.parname)))
    fixed_values = parse.(Float64, strip.(String.(fixedpars_df.value)))

    idx_emax_v = findfirst(==(Symbol("Emax_V")), fixed_names)
    if idx_emax_v === nothing
        error("Fixed parameter table must contain Emax_V")
    end

    Vvals_all = filter(:Quantity => ==("LogVirusLoad"), fitdata)
    Innvals_all = filter(:Quantity => ==("IL6"), fitdata)
    Symvals_all = filter(:Quantity => ==("Weight"), fitdata)
    Vmax = maximum(Vvals_all.Value)
    Innmax = maximum(Innvals_all.Value)
    Symmax = maximum(Symvals_all.Value)

    Y0 = (Ad = 0.0, Ac = 0.0, At = 0.0, U = 1.0e7, I = 0.0, V = 1.0,
        F = 0.0, A = 0.0, S = 0.0)

    par_ini_full = (
        b = 1.0e-8,
        k = 1.0e-4,
        p = 2.0e3,
        kF = 1.0,
        cV = 10.0,
        gF = 0.1,
        hV = 1.0e4,
        Fmax = 5.0,
        hF = 1.0,
        gS = 10.0,
        cS = 1.0,
        Emax_F = 0.5,
        C50_F = 1.0,
        C50_V = 1.0,
    )

    lb = [1.0e-12, 1.0e-8, 1.0e-3, 1.0e-10, 1.0e-2, 1.0e-3, 1.0e-3,
        1.0, 1.0e-4, 1.0e-4, 1.0e-2, 1.0e-3, 1.0e-7, 1.0e-7]
    ub = [1.0e-2, 1.0e5, 1.0e10, 1.0e3, 1.0e3, 1.0e3, 1.0e8,
        1.0e2, 1.0e5, 1.0e4, 1.0e4, 1.0, 1.0e6, 1.0e6]

    par_ini = collect(values(par_ini_full))
    _clamp!(par_ini, lb, ub)
    param_names = collect(keys(par_ini_full))

    parlabels = [
        "Virus infection rate",
        "Adaptive response clearance rate",
        "Virus production rate",
        "Innate response supression strength",
        "Virus removal rate",
        "Maximum innate response induction",
        "Half maximum of innate response induction",
        "Maximum innate response",
        "Innate response decay rate",
        "Symptom induction rate",
        "Symptom decay rate",
        "Maximum innate response supression",
        "Half maximum of innate response effect",
        "Half maximum of virus suppression effect",
    ]

    samples = Vector{Vector{Float64}}()
    push!(samples, copy(fixed_values))

    if nsamp > 0
        lower = 0.5 .* fixed_values
        upper = 2.0 .* fixed_values
        sample_mat = QuasiMonteCarlo.sample(
            QuasiMonteCarlo.LatinHypercubeSample(), lower, upper, nsamp; rng = rng)
        for j in 1:nsamp
            push!(samples, collect(sample_mat[:, j]))
        end
    end

    for vals in samples
        vals[idx_emax_v] = 1.0
    end

    solver = _nlopt_solver(algorithm)
    scenario_levels = String.(levels(fitdata.Scenario))
    doses = [dose_map[scen] / 50.0 for scen in scenario_levels]

    results = Vector{FitResult}(undef, length(samples))
    start_time = Dates.now()

    for (idx, fixed_vals) in enumerate(samples)
        if verbose
            println("processing sample $(idx) at $(Dates.format(Dates.now(), "HH:MM:SS"))")
        end

        objective = _objective_factory(fitdata, Y0, fixed_names, fixed_vals,
            param_names, tfinal, dt, doses, scenario_levels,
            Vmax, Innmax, Symmax)
        optf = OptimizationFunction((x, _) -> objective(x), Optimization.NoGradient())
        prob = OptimizationProblem(optf, par_ini; lb = lb, ub = ub)
        sol = Optimization.solve(prob, solver; maxeval = maxeval, ftol_rel = ftol_rel)

        params = collect(sol.minimizer)
        params_nt = NamedTuple{Tuple(param_names)}(Tuple(params))
        fixed_nt = namedtuple_from(fixed_names, fixed_vals)

        tfinal_sim = 7.5
        dt_sim = 0.01
        odeout_list = Vector{DataFrame}()
        for (dose, scen) in zip(doses, scenario_levels)
            kwargs = merge(Y0, params_nt, fixed_nt,
                (Ad0 = dose, txstart = 1.0, txinterval = 0.5, txend = 4.0,
                    tstart = 0.0, tfinal = tfinal_sim, dt = dt_sim))
            push!(odeout_list, simulate_model(; kwargs...))
        end

        annotated = DataFrame()
        dose_values_mgkg = [dose_map[scen] for scen in scenario_levels]
        for (df, scen, dose_mgkg) in zip(odeout_list, scenario_levels, dose_values_mgkg)
            tmp = copy(df)
            tmp.Dose = fill(dose_mgkg, nrow(tmp))
            tmp.Scenario = fill(scen, nrow(tmp))
            annotated = vcat(annotated, tmp)
        end

        results[idx] = FitResult(params, sol.minimum, copy(par_ini),
            copy(param_names), dict_from(fixed_names, fixed_vals), Y0,
            copy(fitdata), parlabels, annotated, annotated, algorithm)
    end

    if save_results
        mkpath(dirname(result_path))
        open(result_path, "w") do io
            Serialization.serialize(io, results)
        end
    end

    if verbose
        runtime = Dates.now() - start_time
        println("Model fitting completed in $(Dates.value(runtime) / 60_000) minutes.")
    end

    return results
end

end # module
