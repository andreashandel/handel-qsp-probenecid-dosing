#!/usr/bin/env julia

# This script is designed to be used as a teaching aid.  Every section of the
# code contains detailed commentary explaining not only *what* each command
# does, but *why* it is needed when working with the quantitative systems
# pharmacology (QSP) model implemented in this repository.  A beginner should
# be able to follow along and experiment with the simulation without needing to
# understand the separate parameter fitting workflow.

# Import Julia's package manager so we can activate the local project
# environment.  Activating ensures all dependencies specified in `Project.toml`
# are available when the script is executed from the command line.
using Pkg

# The simulation utilities live inside the `code/julia` project directory.  We
# activate that environment (the `..` path relative to this file) so that all
# required packages—such as DifferentialEquations.jl and DataFrames.jl—are
# loaded from the correct manifest.
Pkg.activate(joinpath(@__DIR__, ".."))

# The actual differential equation model is defined in
# `analysis/SimulateModel.jl`.  `include` loads that file so we can access the
# `simulate_model` function that integrates the ODE system.
include(joinpath(@__DIR__, "..", "analysis", "SimulateModel.jl"))

# Bring the exported `simulate_model` routine into scope.  The leading dot
# tells Julia to look for the module relative to the script's module (this
# file).  After this line we can call `simulate_model` directly.
using .QSPModel: simulate_model

# These packages are used for reading the fixed-parameter CSV file and for
# storing the simulation results in a tabular format that is convenient to
# inspect or save to disk.
using CSV
using DataFrames

# -----------------------------------------------------------
# Command-line interface helpers
# -----------------------------------------------------------

# Explain to the user how to run the script.  This function prints a friendly
# usage message and is triggered when the `--help` flag is supplied.
function print_help()
    println("Explore the QSP simulation without running the parameter fitting pipeline.")
    println("\nOptions:")
    println("  --dose <float>         Drug amount added to the depot compartment per dose (mg, default: 100.0)")
    println("  --txstart <float>      Time of the first dose in days (default: 1.0)")
    println("  --txinterval <float>   Time between doses in days (default: 0.5)")
    println("  --txend <float>        Time of the last dose in days (default: 4.0)")
    println("  --duration <float>     Length of the simulation window in days (default: 7.5)")
    println("  --dt <float>           Sampling interval for saved results in days (default: 0.01)")
    println("  --fixed <path>         Path to the fixed parameter CSV (default: data/processed-data/fixed-parameters.csv)")
    println("  --output <path>        Optional path to save the simulated time series as CSV")
    println("  --help                 Show this message and exit")
end

# Parse command-line arguments into a dictionary of string key/value pairs.
# The parser is intentionally minimalistic so that beginners can see exactly
# how arguments flow into the simulation settings.
function parse_args()
    # `kwargs` stores options that expect a value (e.g., `--dose 50`).
    kwargs = Dict{String, String}()
    # Iterate over the `ARGS` vector supplied by Julia when the script is run.
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--help"
            # If the user requests help we display the usage information and
            # exit immediately.  `exit(0)` signals successful termination.
            print_help()
            exit(0)
        elseif startswith(arg, "--")
            # All other arguments are expected to come in `--flag value` pairs.
            if i == length(ARGS)
                error("Missing value for argument $(arg)")
            end
            kwargs[arg] = ARGS[i + 1]
            i += 2
        else
            # Encountering a token that does not start with `--` likely means
            # the user made a typo, so we stop with a descriptive error.
            error("Unknown argument: $(arg)")
        end
    end
    return kwargs
end

# -----------------------------------------------------------
# Data loading utilities
# -----------------------------------------------------------

# Load the pharmacokinetic constants that remain fixed during the simulation
# from the CSV file used by the fitting workflow.  Returning a `NamedTuple`
# lets us merge these parameters directly into the `simulate_model` call.
function load_fixed_parameters(path::AbstractString)
    # Read the CSV into a `DataFrame`.  Each row contains the parameter name and
    # its numerical value.
    df = CSV.read(path, DataFrame)
    # The simulation function expects parameters as symbols, so we convert the
    # string names into `Symbol` objects and collect the associated numeric
    # values.  We also strip whitespace to guard against formatting issues.
    names = Symbol.(strip.(String.(df.parname)))
    values = parse.(Float64, strip.(String.(df.value)))
    # Construct the `NamedTuple` using the helper function from `Base`.  The
    # order of the names and values vectors must align so that each parameter
    # receives the correct value.
    return NamedTuple{Tuple(names)}(Tuple(values))
end

# A tiny helper that computes the trapezoidal rule.  We use it to mimic the
# area-under-the-curve summaries reported elsewhere in the project.  The
# implementation is written out explicitly so beginners can follow the math.
function trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    # Ensure the input vectors have matching lengths before performing arithmetic.
    length(x) == length(y) || error("x and y must have the same length")
    # The trapezoidal rule approximates the integral by summing the area of
    # trapezoids formed between consecutive points.  We loop over each interval
    # and accumulate the contribution.
    total = 0.0
    for i in 1:(length(x) - 1)
        Δx = x[i + 1] - x[i]
        avg_height = (y[i] + y[i + 1]) / 2
        total += Δx * avg_height
    end
    return total
end

# -----------------------------------------------------------
# Main simulation workflow
# -----------------------------------------------------------

function main()
    # Parse the command-line options once at the start of the program.
    kwargs = parse_args()

    # Resolve the repository root so that default file paths work regardless of
    # where the script is executed from.
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))

    # Pull the simulation settings from the parsed arguments, falling back to
    # sensible defaults when the user does not provide a value.  Because all
    # entries are stored as strings we explicitly convert them to `Float64`.
    dose = parse(Float64, get(kwargs, "--dose", "100.0"))
    txstart = parse(Float64, get(kwargs, "--txstart", "1.0"))
    txinterval = parse(Float64, get(kwargs, "--txinterval", "0.5"))
    txend = parse(Float64, get(kwargs, "--txend", "4.0"))
    duration = parse(Float64, get(kwargs, "--duration", "7.5"))
    dt = parse(Float64, get(kwargs, "--dt", "0.01"))

    # Build the path to the fixed-parameter CSV.  Users can override the default
    # to experiment with alternative parameter sets.
    fixed_path = get(kwargs, "--fixed", joinpath(repo_root, "data", "processed-data", "fixed-parameters.csv"))

    # Optionally the user can request that the simulated time series be written
    # to a CSV file for further analysis in spreadsheets or other tools.
    output_path = get(kwargs, "--output", "")

    # Load the fixed pharmacokinetic parameters from disk.
    fixed_nt = load_fixed_parameters(fixed_path)

    # Define the initial state of the system.  The values mirror the defaults
    # used throughout the project and represent the starting viral load,
    # uninfected cells and immune response levels before treatment.
    Y0 = (Ad = 0.0, Ac = 0.0, At = 0.0, U = 1.0e7, I = 0.0, V = 1.0,
        F = 0.0, A = 0.0, S = 0.0)

    # Choose a representative set of dynamic parameters.  These numbers match
    # the initial guesses used in the fitting routine, providing a reasonable
    # qualitative trajectory without requiring optimisation.
    dyn_params = (
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

    # Combine all simulation inputs into a single named tuple using `merge`.
    # This mirrors how the fitting code constructs keyword arguments for
    # `simulate_model`, ensuring both workflows stay consistent.
    sim_kwargs = merge(Y0, dyn_params, fixed_nt, (
        Ad0 = dose,
        txstart = txstart,
        txinterval = txinterval,
        txend = txend,
        tstart = 0.0,
        tfinal = duration,
        dt = dt,
    ))

    # Run the ODE solver.  The result is returned as a `DataFrame` with one row
    # per saved time point and columns for each state variable.
    odeout = simulate_model(; sim_kwargs...)

    # Provide immediate feedback by displaying the first few rows so users can
    # see the structure of the result table.
    println("\nPreview of simulated state variables:")
    show(first(odeout, 5))
    println("\n")

    # Compute simple summary metrics to demonstrate how the raw output can be
    # interpreted.  We calculate peak viral load, the time of that peak and the
    # area under the viral-load curve (AUC) after applying a log10 transform.
    peak_index = argmax(odeout.V)
    peak_time = odeout.time[peak_index]
    peak_virus = odeout.V[peak_index]
    auc_v = trapz(odeout.time, log10.(max.(1.0, odeout.V)))

    println("Peak viral load: $(peak_virus) copies/mL at day $(round(peak_time; digits = 2))")
    println("Log10 viral load AUC: $(round(auc_v; digits = 3))")

    # When the user requests an output file we save the entire trajectory.  The
    # surrounding `if` guard prevents unnecessary disk writes during quick tests.
    if !isempty(output_path)
        CSV.write(output_path, odeout)
        println("\nSaved time series to $(output_path)")
    end

    println("\nSimulation complete.  Modify the command-line options to explore different dosing schedules!")
end

# Only run `main` when this script is executed directly.  This allows the file
# to be `include`d from other Julia code without automatically starting a
# simulation, which can be useful for interactive notebook sessions.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
