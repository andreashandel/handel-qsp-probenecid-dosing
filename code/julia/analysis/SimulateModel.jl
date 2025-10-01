module QSPModel

using DifferentialEquations
using DiffEqCallbacks
using DataFrames

export simulate_model

"""
    simulate_model(; kwargs...) -> DataFrame

Integrate the QSP ordinary differential equation model describing SARS-CoV-2
infection dynamics and the effect of the investigational compound.  The
function mirrors the behaviour of the original R implementation in
`model-simulator-function.R`.

# Keyword arguments
- `Ad`, `Ac`, `At`, `U`, `I`, `V`, `F`, `A`, `S`: initial conditions for the
  depot, central, target drug compartments and host response variables.
- `b`, `cI`, `k`, `p`, `kF`, `cV`, `gF`, `hV`, `Fmax`, `cF`, `hF`, `gA`, `gS`,
  `cS`, `Emax_V`, `C50_V`, `Emax_F`, `C50_F`: dynamical parameters to be
  estimated or fixed during fitting.
- `ka`, `Vc`, `Vt`, `Q`, `Vmax`, `Km`: pharmacokinetic parameters (typically
  fixed) controlling drug absorption and distribution.
- `txstart`, `txinterval`, `txend`: treatment schedule (days).  When
  `txstart > txend` no doses are administered.
- `Ad0`: per-dose amount (mg) that is added to the depot compartment.
- `tstart`, `tfinal`: simulation window (days).
- `dt`: output sampling interval (days).  The solver uses adaptive steps and
  is forced to save results on this grid to match the R workflow.
- `solver`: any `OrdinaryDiffEq.jl` compatible integrator (defaults to
  `Tsit5()`).
- `abstol`, `reltol`: integration tolerances.

# Returns
A `DataFrame` with the same nine state variables and time column as the R
function.  The columns are ordered as `[:time, :Ad, :Ac, :At, :U, :I, :V, :F,
:A, :S]`.
"""
function simulate_model(; Ad::Float64 = 0.0, Ac::Float64 = 0.0, At::Float64 = 0.0,
    U::Float64 = 1.0e7, I::Float64 = 0.0, V::Float64 = 1.0, F::Float64 = 0.0,
    A::Float64 = 1.0, S::Float64 = 0.0, b::Float64 = 1.0e-8,
    cI::Float64 = 1.0, k::Float64 = 1.0e-6, p::Float64 = 1.0e3,
    kF::Float64 = 1.0e-3, cV::Float64 = 10.0, gF::Float64 = 1.0,
    hV::Float64 = 1.0e4, Fmax::Float64 = 10.0, cF::Float64 = 5.0,
    hF::Float64 = 10.0, gA::Float64 = 1.0, gS::Float64 = 1.0,
    cS::Float64 = 2.0, Emax_V::Float64 = 1.0, C50_V::Float64 = 10.0,
    Emax_F::Float64 = 0.5, C50_F::Float64 = 10.0, ka::Float64 = 2.0,
    Vc::Float64 = 4.0, Vt::Float64 = 20.0, Q::Float64 = 5.0,
    Vmax::Float64 = 1.6, Km::Float64 = 60.0, txstart::Float64 = 1.0,
    txinterval::Float64 = 0.5, txend::Float64 = 4.0, Ad0::Float64 = 100.0,
    tstart::Float64 = 0.0, tfinal::Float64 = 10.0, dt::Float64 = 0.01,
    solver = DifferentialEquations.Tsit5(), abstol::Float64 = 1.0e-8,
    reltol::Float64 = 1.0e-8)

    state0 = [Ad, Ac, At, U, I, V, F, A, S]
    pars = (; b, cI, k, p, kF, cV, gF, hV, Fmax, cF, hF, gA, gS, cS,
        Emax_V, C50_V, Emax_F, C50_F, ka, Vc, Vt, Q, Vmax, Km)

    function qsp_rhs!(du, u, p, t)
        Ad, Ac, At, U, I, V, F, A, S = u
        (; b, cI, k, p, kF, cV, gF, hV, Fmax, cF, hF, gA, gS, cS,
            Emax_V, C50_V, Emax_F, C50_F, ka, Vc, Vt, Q, Vmax, Km) = p

        du[1] = -ka * Ad
        du[2] = ka * Ad - Q / Vc * Ac + Q / Vt * At - Vmax * Ac / (Km * Vc + Ac)
        du[3] = Q / Vc * Ac - Q / Vt * At

        Ct = At / Vt
        f_V = Emax_V * Ct / (C50_V + Ct)
        f_F = Emax_F * Ct / (C50_F + Ct)

        du[4] = -b * U * V
        du[5] = b * U * V - cI * I - k * I * A
        du[6] = (1 - f_V) * p * I / (1 + kF * F) - cV * V
        du[7] = (1 - f_F) * gF * (V / (V + hV)) * (Fmax - F) - cF * F
        du[8] = V * F / (V * F + hF) + gA * A
        du[9] = gS * F - cS * S

        return nothing
    end

    affect!(integrator) = (integrator.u[1] += Ad0)
    dosing_times = if txstart <= txend
        collect(txstart:txinterval:txend)
    else
        Float64[]
    end
    callback = isempty(dosing_times) ? nothing : DiffEqCallbacks.PresetTimeCallback(dosing_times, affect!)

    tspan = (tstart, tfinal)
    prob = DifferentialEquations.ODEProblem(qsp_rhs!, state0, tspan, pars)
    save_times = collect(range(tstart, tfinal; step = dt))

    sol = DifferentialEquations.solve(prob, solver;
        saveat = save_times, abstol = abstol, reltol = reltol,
        callback = callback)

    if sol.retcode != DifferentialEquations.ReturnCode.Success
        error("ODE solver failed with return code $(sol.retcode)")
    end

    arr = reduce(hcat, sol.u)'
    names = [:Ad, :Ac, :At, :U, :I, :V, :F, :A, :S]
    df = DataFrames.DataFrame(time = sol.t)
    for (idx, name) in enumerate(names)
        df[!, name] = arr[:, idx]
    end
    return df
end

end # module
