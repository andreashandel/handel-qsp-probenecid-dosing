"""
ODE simulators for model1 and model2
====================================

This file ports the two R simulator functions to DifferentialEquations.jl.

Important implementation notes
------------------------------
1. State ordering is explicit (`state_symbols`) so we can safely read/write by
   index while keeping equations readable.
2. Dosing is implemented through a `PresetTimeCallback`, which matches the R
   event logic that adds drug mass to the depot compartment at scheduled times.
3. We enforce `maxiters = 20_000` as requested to avoid very long integrator
   runs in pathological parameter regions.
"""

"""
Simple trapezoidal integration utility used in dose-prediction summaries.
"""
function trapz_area(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    n == length(y) || error("x and y must have same length")
    n < 2 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        dx = Float64(x[i + 1] - x[i])
        acc += dx * (Float64(y[i]) + Float64(y[i + 1])) / 2.0
    end
    return acc
end

"""
Build a merged parameter dictionary from fitted and fixed parameter containers.

Both inputs are expected to be `Dict{String,Float64}` and may contain sigma
parameters. Sigma values are kept in the merged dictionary because objective
computation needs them; ODE RHS simply ignores sigma keys.
"""
function merged_parameters(fitpars::Dict{String,Float64}, fixedpars::Dict{String,Float64})
    out = Dict{String,Float64}()
    for (k, v) in fixedpars
        out[k] = v
    end
    for (k, v) in fitpars
        out[k] = v
    end
    return out
end

"""
Generate dosing times from schedule settings.
"""
function dosing_times(txstart::Real, txend::Real, txinterval::Real)
    if txend < txstart
        return Float64[]
    end
    return collect(Float64(txstart):Float64(txinterval):Float64(txend))
end

"""
Return the solver object corresponding to `solvertype`.

Default behavior is intentionally stiff-first (`Rodas5P`) because model fitting
often explores parameter regions that are numerically stiff. This reduces
`maxiters` exits compared with non-stiff defaults.
"""
function choose_solver(solvertype::String)
    s = lowercase(strip(solvertype))
    if s in ("vode", "lsoda", "stiff", "bdf", "rodas5", "rodas5p")
        return Rodas5P()
    elseif s == "trbdf2"
        return TRBDF2()
    elseif s == "rosenbrock23"
        return Rosenbrock23()
    elseif s == "kencarp4"
        return KenCarp4()
    elseif s == "tsit5"
        return Tsit5()
    end
    @warn "Unknown solvertype '$solvertype'. Falling back to Rodas5P()."
    return Rodas5P()
end

"""
Simulate model1 for one scenario/dose.

Inputs
------
- `allpars`: merged dictionary containing model parameters and dosing schedule.
- `times`: output times (must include tstart).

Returns
-------
Named tuple with:
- `df`: simulation DataFrame.
- `ok`: boolean flag indicating a successful solve.
- `retcode`: solver return code for diagnostics.
"""
function simulate_model1(allpars::Dict{String,Float64}; times::Vector{Float64}, solvertype::String = "vode")
    state_symbols = [:Ad, :Ac, :At, :U, :I, :V, :F, :A, :S]
    y0 = [allpars[string(s)] for s in state_symbols]

    function rhs!(du, u, p, t)
        Ad, Ac, At, U, I, V, F, A, S = u
        b = p["b"]; cI = p["cI"]; k = p["k"]; pV = p["p"]
        kF = p["kF"]; cV = p["cV"]; gF = p["gF"]; hV = p["hV"]
        Fmax = p["Fmax"]; cF = p["cF"]; hF = p["hF"]; gA = p["gA"]
        gS = p["gS"]; cS = p["cS"]
        Emax_V = p["Emax_V"]; C50_V = p["C50_V"]; Emax_F = p["Emax_F"]; C50_F = p["C50_F"]
        ka = p["ka"]; Vc = p["Vc"]; Vt = p["Vt"]; Q = p["Q"]
        Vmax = p["Vmax"]; Km = p["Km"]; fmax = p["fmax"]; f50 = p["f50"]

        Ct = At / Vt
        fu = fmax * Ct / (f50 + Ct)
        Cu = fu * Ct
        fV = Emax_V * Cu / (C50_V + Cu)
        fF = Emax_F * Cu / (C50_F + Cu)

        du[1] = -ka * Ad
        du[2] = ka * Ad - Q / Vc * Ac + Q / Vt * At - Vmax * Ac / (Km * Vc + Ac)
        du[3] = Q / Vc * Ac - Q / Vt * At
        du[4] = -b * U * V
        du[5] = b * U * V - cI * I - k * I * A
        du[6] = (1.0 - fV) * pV * I / (1.0 + kF * F) - cV * V

        # Matches current R model1 implementation.
        du[7] = (1.0 - fF) * gF * V / (V + hV) * (Fmax - F) - cF * F
        du[8] = F / (F + hF) + gA * A
        du[9] = gS * F - cS * S
    end

    tspan = (allpars["tstart"], allpars["tfinal"])
    prob = ODEProblem(rhs!, y0, tspan, allpars)

    dts = dosing_times(allpars["txstart"], allpars["txend"], allpars["txinterval"])
    ad_index = findfirst(==(Symbol("Ad")), state_symbols)
    dose_affect! = function(integrator)
        integrator.u[ad_index] += allpars["Ad0"] / 50.0
    end
    cb = isempty(dts) ? CallbackSet() : PresetTimeCallback(dts, dose_affect!)

    sol = solve(
        prob,
        choose_solver(solvertype),
        callback = cb,
        abstol = allpars["tols"],
        reltol = allpars["tols"],
        saveat = times,
        maxiters = 20_000,
        verbose = false,
    )

    ok = SciMLBase.successful_retcode(sol.retcode)
    df = DataFrame(time = Array(sol.t))
    for (i, s) in enumerate(state_symbols)
        df[!, String(s)] = Array(sol[i, :])
    end
    return (df = df, ok = ok, retcode = sol.retcode)
end

"""
Simulate model2 for one scenario/dose.

The structure mirrors `simulate_model1`; only state definition and ODE terms
that differ in model2 are changed.
"""
function simulate_model2(allpars::Dict{String,Float64}; times::Vector{Float64}, solvertype::String = "vode")
    state_symbols = [:Ad, :Ac, :At, :U, :E, :I, :V, :F, :A, :S]
    y0 = [allpars[string(s)] for s in state_symbols]

    function rhs!(du, u, p, t)
        Ad, Ac, At, U, E, I, V, F, A, S = u
        b = p["b"]; cE = p["cE"]; cI = p["cI"]; k = p["k"]; pV = p["p"]
        kF = p["kF"]; cV = p["cV"]; gF = p["gF"]; hV = p["hV"]
        Fmax = p["Fmax"]; cF = p["cF"]; hF = p["hF"]; gA = p["gA"]
        gS = p["gS"]; cS = p["cS"]
        Emax_V = p["Emax_V"]; C50_V = p["C50_V"]; Emax_F = p["Emax_F"]; C50_F = p["C50_F"]
        ka = p["ka"]; Vc = p["Vc"]; Vt = p["Vt"]; Q = p["Q"]
        Vmax = p["Vmax"]; Km = p["Km"]; fmax = p["fmax"]; f50 = p["f50"]

        Ct = At / Vt
        fu = fmax * Ct / (f50 + Ct)
        Cu = fu * Ct
        fV = Emax_V * Cu / (C50_V + Cu)
        fF = Emax_F * Cu / (C50_F + Cu)
        logV = transform_virus(V)

        du[1] = -ka * Ad
        du[2] = ka * Ad - Q / Vc * Ac + Q / Vt * At - Vmax * Ac / (Km * Vc + Ac)
        du[3] = Q / Vc * Ac - Q / Vt * At
        du[4] = -b * U * V
        du[5] = b * U * V - cE * E
        du[6] = cE * E - cI * I - k * I * A
        du[7] = (1.0 - fV) * pV * I / (1.0 + kF * F) - cV * V
        du[8] = (1.0 - fF) * gF * logV / (logV + hV) * (Fmax - F) - cF * F
        du[9] = logV * F / (logV * F + hF) + gA * A
        du[10] = gS * F * logV - cS * S
    end

    tspan = (allpars["tstart"], allpars["tfinal"])
    prob = ODEProblem(rhs!, y0, tspan, allpars)

    dts = dosing_times(allpars["txstart"], allpars["txend"], allpars["txinterval"])
    ad_index = findfirst(==(Symbol("Ad")), state_symbols)
    dose_affect! = function(integrator)
        integrator.u[ad_index] += allpars["Ad0"] / 50.0
    end
    cb = isempty(dts) ? CallbackSet() : PresetTimeCallback(dts, dose_affect!)

    sol = solve(
        prob,
        choose_solver(solvertype),
        callback = cb,
        abstol = allpars["tols"],
        reltol = allpars["tols"],
        saveat = times,
        maxiters = 20_000,
        verbose = false,
    )

    ok = SciMLBase.successful_retcode(sol.retcode)
    df = DataFrame(time = Array(sol.t))
    for (i, s) in enumerate(state_symbols)
        df[!, String(s)] = Array(sol[i, :])
    end
    return (df = df, ok = ok, retcode = sol.retcode)
end

"""
Single dispatch wrapper that picks model1 or model2 simulator.
"""
function simulate_model(model_choice::String, allpars::Dict{String,Float64}; times::Vector{Float64}, solvertype::String = "vode")
    if model_choice == "model1"
        return simulate_model1(allpars; times = times, solvertype = solvertype)
    elseif model_choice == "model2"
        return simulate_model2(allpars; times = times, solvertype = solvertype)
    else
        error("Unknown model choice: $model_choice")
    end
end
