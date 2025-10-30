#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "analysis", "FitModel.jl"))
using .ModelFitting

function print_help()
    println("Run the QSP model calibration using NLopt.")
    println("\nOptions:")
    println("  --nsamp <int>         Number of Latin hypercube samples for fixed parameters (default: 0)")
    println("  --algorithm <name>    NLopt algorithm (LN_COBYLA, LN_NELDERMEAD, LN_SBPLX)")
    println("  --maxeval <int>       Maximum optimizer evaluations (default: 1000)")
    println("  --ftol <float>        Relative function tolerance (default: 1e-10)")
    println("  --data <path>         Path to processed data CSV")
    println("  --fixed <path>        Path to fixed parameter CSV")
    println("  --output <path>       Output file for serialized results (default: results/output/bestfit.jls)")
    println("  --no-save             Do not write results to disk")
    println("  --help                Display this message")
end

function parse_args()
    kwargs = Dict{String, String}()
    flags = Set{String}()
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--help"
            print_help()
            exit(0)
        elseif arg == "--no-save"
            push!(flags, arg)
            i += 1
        elseif startswith(arg, "--")
            if i == length(ARGS)
                error("Missing value for argument $(arg)")
            end
            kwargs[arg] = ARGS[i + 1]
            i += 2
        else
            error("Unknown argument: $(arg)")
        end
    end
    return kwargs, flags
end

kwargs, flags = parse_args()

nsamp = parse(Int, get(kwargs, "--nsamp", "0"))
alg = Symbol(get(kwargs, "--algorithm", "LN_COBYLA"))
maxeval = parse(Int, get(kwargs, "--maxeval", "1000"))
ftol = parse(Float64, get(kwargs, "--ftol", "1e-10"))
data_path = get(kwargs, "--data", "data/processed-data/processeddata.csv")
fixed_path = get(kwargs, "--fixed", "data/processed-data/fixed-parameters.csv")
output_path = get(kwargs, "--output", "results/output/bestfit.jls")
save_results = !("--no-save" in flags)

fit_model(; data_path = data_path, fixed_path = fixed_path,
    result_path = output_path, nsamp = nsamp, algorithm = alg,
    maxeval = maxeval, ftol_rel = ftol, save_results = save_results)
