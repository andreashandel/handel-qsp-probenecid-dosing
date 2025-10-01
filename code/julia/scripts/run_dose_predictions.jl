#!/usr/bin/env julia
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "analysis", "FitModel.jl"))
include(joinpath(@__DIR__, "..", "analysis", "DosePredictions.jl"))
using .ModelFitting: FitResult
using .DosePredictions
using Serialization

function print_help()
    println("Simulate dose-response scenarios for fitted parameter sets.")
    println("\nOptions:")
    println("  --bestfit <path>   Path to serialized fit results (default: results/output/bestfit.jls)")
    println("  --output <path>    Output file for dose predictions (default: results/output/dose-response-results.jls)")
    println("  --help             Show this message")
end

function parse_args()
    kwargs = Dict{String, String}()
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--help"
            print_help()
            exit(0)
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
    return kwargs
end

kwargs = parse_args()

bestfit_path = get(kwargs, "--bestfit", "results/output/bestfit.jls")
output_path = get(kwargs, "--output", "results/output/dose-response-results.jls")

open(bestfit_path, "r") do io
    bestfit_list = Serialization.deserialize(io)
    run_dose_predictions(bestfit_list; result_path = output_path)
end
