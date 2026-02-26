#!/usr/bin/env julia

"""
Julia workflow entry point: fit model parameters
================================================

This script is the Julia counterpart of `code/analysis-code/run-fit.R`.

How to use
----------
1. Start Julia from repository root.
2. Run:
   `julia --project=code/julia code/julia/scripts/run-fit.jl`
   or activate julia from Positron by rnning first line, then run script from julia console: 
   include("code/julia/scripts/run-fit.jl")


Configuration strategy
----------------------
The script intentionally keeps all user-editable settings in one block below.
If you are migrating from R, start by changing only:
- `model_choice`
- `run_multistart_stage`
- `run_sampling_stage`
- optimizer budgets (`global_maxeval`, `local_maxiters`)
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..")) #one up from current file
#Pkg.activate(joinpath(@__DIR__, ".")) #same folder as current file
Pkg.instantiate() #install packages as needed

include(joinpath(@__DIR__, "..", "src", "handel-qsp.jl"))
import .HandelQSP

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choice = "model1"                 # "model1" or "model2"

# Stage toggles:
# - run_multistart_stage = true,  run_sampling_stage = false -> multistart only
# - run_multistart_stage = false, run_sampling_stage = true  -> sampling only
# - run_multistart_stage = true,  run_sampling_stage = true  -> both
run_multistart_stage = false
run_sampling_stage = true

# Parallelism / compute budget.
# Note: Julia thread count is set at process launch via JULIA_NUM_THREADS.
# `n_workers` is consumed by the workflow where parallel sections exist.
n_workers = 25
verbose_fit_log = true
log_top_n = 5
reuse_previous_bestfit = true
previous_bestfit_file = joinpath(HandelQSP.project_root(), "results", "output", "$(model_choice)-bestfit-multistart.jld2")
previous_bestfit_n = 5

# Objective + parameter handling.
sigma_to_fit = String[]                 # e.g. ["sigma_add_VirusLoad"]
weight_mode = "equal_quantity"          # "equal_quantity" or "per_block"
user_fixed_params = Dict{String,Float64}() # e.g. Dict("Emax_F" => 1.0)
use_log_space = true                    # fit positive parameters in log space

# ODE solver settings.
solver_settings = (
    solvertype = "rodas5p",             # stiff: "rodas5p", "trbdf2", "rosenbrock23"
    tols = 1e-10,
    tfinal = 7.0,
    dt = 0.1,
)

# Global optimizer settings
# Backends:
#   :blackboxoptim -> Differential Evolution (:adaptive_de_rand_1_bin_radiuslimited)
#   :nlopt        -> NLopt global algorithms (set `global_nlopt_algorithm`)
global_optimizer = :blackboxoptim
global_maxeval = 100_000                 # max objective evaluations in stage-1 global search
global_population = 200                  # BlackBoxOptim DE population size
global_nthreads = max(1, min(n_workers, Threads.nthreads())) # used by BlackBoxOptim DE
global_trace_mode = :compact            # :silent, :compact, or :verbose
global_trace_interval = 120.0           # seconds between global progress lines
# NLopt global algorithm options:
# use NLopt enum symbols directly (e.g. :GN_ESCH, :GN_MLSL_LDS, :GN_ISRES, :GN_DIRECT)
global_nlopt_algorithm = :GN_ESCH       # used only when global_optimizer = :nlopt
global_nlopt_population = 200             # 0 = NLopt default; >0 sets population when supported

# Local settings
n_local_restarts = 10                   # global-seeded local starts: 1 global best + (n_local_restarts-1) jittered
local_jitter_scale = 0.1                # jitter multiplier on parameter span for DE-neighbor starts
local_maxiters = 2000                  # per-local-refinement max iterations
# Available NLopt local optimizer symbols:
# use NLopt enum symbols directly (e.g. :LN_BOBYQA, :LN_SBPLX, :LN_NELDERMEAD, :LN_COBYLA)
local_optimizer = :LN_COBYLA            # baseline local optimizer
local_show_trace = false

# Sampling-stage settings (used only if run_sampling_stage=true).
# This stage runs LOCAL-ONLY fitting per fixed-parameter sample, seeded from
# the best multistart fit (from this run or from saved multistart output).
nsamp = 40                              # number of LHS fixed-parameter samples
fixed_overrides = Dict("Emax_V" => 1.0) # forced fixed values for all samples
sample_lower_factor = 0.5               # LHS lower multiplier for fixed-parameter sampling
sample_upper_factor = 2.0               # LHS upper multiplier for fixed-parameter sampling
sample_seed = 1234                      # reproducible sampling seed
sampling_use_log_space = true           # local-only sampling fit in log space
sampling_local_maxiters = 2000         # sampling stage local optimizer iterations
# Available NLopt local optimizer symbols:
# use NLopt enum symbols directly (e.g. :LN_BOBYQA, :LN_SBPLX, :LN_NELDERMEAD, :LN_COBYLA)
sampling_local_optimizer = :LN_COBYLA   # sampling stage optimizer
sampling_local_show_trace = false

# Output location.
output_dir = joinpath(HandelQSP.project_root(), "results", "output")

# Quick smoke-test profile for fast end-to-end validation.
# Set to `false` for full local-restart counts (e.g., 10 DE-seeded starts).
quick_smoke_test = false

settings = merge(
    HandelQSP.default_fit_settings(model_choice),
    (
        model_choice = model_choice,
        run_multistart_stage = run_multistart_stage,
        run_sampling_stage = run_sampling_stage,
        n_workers = n_workers,
        verbose_fit_log = verbose_fit_log,
        log_top_n = log_top_n,
        reuse_previous_bestfit = reuse_previous_bestfit,
        previous_bestfit_file = previous_bestfit_file,
        previous_bestfit_n = previous_bestfit_n,
        sigma_to_fit = sigma_to_fit,
        weight_mode = weight_mode,
        user_fixed_params = user_fixed_params,
        solver_settings = solver_settings,
        use_log_space = use_log_space,
        global_optimizer = global_optimizer,
        global_maxeval = global_maxeval,
        global_population = global_population,
        global_nthreads = global_nthreads,
        global_trace_mode = global_trace_mode,
        global_trace_interval = global_trace_interval,
        global_nlopt_algorithm = global_nlopt_algorithm,
        global_nlopt_population = global_nlopt_population,
        n_local_restarts = n_local_restarts,
        local_jitter_scale = local_jitter_scale,
        local_maxiters = local_maxiters,
        local_optimizer = local_optimizer,
        local_show_trace = local_show_trace,
        nsamp = nsamp,
        fixed_overrides = fixed_overrides,
        sample_lower_factor = sample_lower_factor,
        sample_upper_factor = sample_upper_factor,
        sample_seed = sample_seed,
        sampling_use_log_space = sampling_use_log_space,
        sampling_local_maxiters = sampling_local_maxiters,
        sampling_local_optimizer = sampling_local_optimizer,
        sampling_local_show_trace = sampling_local_show_trace,
        output_dir = output_dir,
    ),
)

if quick_smoke_test
    settings = merge(
        settings,
        (
            # Smaller budgets to keep runtime short while still exercising
            # global + local optimizer paths.
            global_maxeval = 1000,
            global_population = 20,
            local_maxiters = 200,
            sampling_local_maxiters = 200,
            # Slightly looser tolerances and coarser output grid for speed.
            solver_settings = (solvertype = "rodas5p", tols = 1e-8, tfinal = 7.0, dt = 0.2),
        ),
    )
end

if settings.n_workers > Threads.nthreads()
    @warn "Requested n_workers exceeds available Julia threads; launch Julia with larger JULIA_NUM_THREADS if needed." requested = settings.n_workers available = Threads.nthreads()
end

HandelQSP.run_fit_workflow(settings)
