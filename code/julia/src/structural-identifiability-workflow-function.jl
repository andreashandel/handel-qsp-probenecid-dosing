"""
Structural identifiability workflow for fitted ODE parameters
==============================================================

This file adds a reproducible wrapper around StructuralIdentifiability.jl for
the two project models. It focuses on the parameters currently fitted in
`run-fit.jl`, while substituting all fixed parameters from the existing CSV
files.

Important modeling notes
------------------------
1. The fitting workflow transforms viral load for the objective. For
   identifiability, model1 uses `V` directly as the virus output because the
   transform is invertible.
2. Model2 uses `log10(V + 1)` inside the ODE right-hand side, which is not
   directly rational. We represent this with an auxiliary state `LV` and
   `LV' = V' / (log10_base * (V + 1))`, then keep `log10_base` fixed. This
   yields a rational ODE form compatible with StructuralIdentifiability.jl.
"""

using Logging
using StructuralIdentifiability

const MODEL1_SI = StructuralIdentifiability.@ODEmodel(
    Ad'(t) = -ka * Ad(t) + u(t),
    Ac'(t) = ka * Ad(t) - Q / Vc * Ac(t) + Q / Vt * At(t) - Vmax * Ac(t) / (Km * Vc + Ac(t)),
    At'(t) = Q / Vc * Ac(t) - Q / Vt * At(t),
    U'(t) = -b * U(t) * V(t),
    I'(t) = b * U(t) * V(t) - cI * I(t) - k * I(t) * A(t),
    V'(t) = (1 - Emax_V * (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))) / (C50_V + (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))))) * p * I(t) / (1 + kF * F(t)) - cV * V(t),
    F'(t) = (1 - Emax_F * (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))) / (C50_F + (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))))) * gF * V(t) / (V(t) + hV) * (Fmax - F(t)) - cF * F(t),
    A'(t) = F(t) / (F(t) + hF) + gA * A(t),
    S'(t) = gS * F(t) - cS * S(t),
    yv(t) = V(t),
    yf(t) = F(t),
    ys(t) = S(t)
)

const MODEL2_SI = StructuralIdentifiability.@ODEmodel(
    Ad'(t) = -ka * Ad(t) + u(t),
    Ac'(t) = ka * Ad(t) - Q / Vc * Ac(t) + Q / Vt * At(t) - Vmax * Ac(t) / (Km * Vc + Ac(t)),
    At'(t) = Q / Vc * Ac(t) - Q / Vt * At(t),
    U'(t) = -b * U(t) * V(t),
    E'(t) = b * U(t) * V(t) - cE * E(t),
    I'(t) = cE * E(t) - cI * I(t) - k * I(t) * A(t),
    V'(t) = (1 - Emax_V * (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))) / (C50_V + (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))))) * p * I(t) / (1 + kF * F(t)) - cV * V(t),
    LV'(t) = ((1 - Emax_V * (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))) / (C50_V + (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))))) * p * I(t) / (1 + kF * F(t)) - cV * V(t)) / (log10_base * (V(t) + 1)),
    F'(t) = (1 - Emax_F * (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))) / (C50_F + (fmax * At(t)^2 / (Vt * (At(t) + Vt * f50))))) * gF * LV(t) / (LV(t) + hV) * (Fmax - F(t)) - cF * F(t),
    A'(t) = LV(t) * F(t) / (LV(t) * F(t) + hF) + gA * A(t),
    S'(t) = gS * F(t) * LV(t) - cS * S(t),
    yv(t) = LV(t),
    yf(t) = F(t),
    ys(t) = S(t)
)

"""
Default settings for structural identifiability analysis.
"""
function default_structural_identifiability_settings()
    return (
        model_choices = ["model1", "model2"],
        prob_threshold = 0.99,
        loglevel = Logging.Warn,
        rationalize_tol = 1e-12,
        output_dir = repo_path("results", "output", "identifiability"),
        # Fixed value used in model2 auxiliary `LV` state:
        # LV = log10(V + 1)
        log10_base_value = log(10.0),
    )
end

"""
Map ODE parameter objects by name for robust lookup.
"""
function si_parameter_lookup(ode)
    return Dict(string(p) => p for p in ode.parameters)
end

"""
Build substitution dictionary for fixed parameters present in `ode`.
"""
function si_fixed_parameter_substitutions(
    ode,
    fixedpars::Dict{String,Float64};
    extra_fixed::Dict{String,Float64} = Dict{String,Float64}(),
    tol::Float64 = 1e-12,
)
    lookup = si_parameter_lookup(ode)
    subs = Dict{Any,Any}()
    used_names = String[]

    for (k, v) in fixedpars
        haskey(lookup, k) || continue
        subs[lookup[k]] = rationalize(v; tol = tol)
        push!(used_names, k)
    end
    for (k, v) in extra_fixed
        haskey(lookup, k) || continue
        subs[lookup[k]] = rationalize(v; tol = tol)
        push!(used_names, k)
    end

    return (subs = subs, used_names = sort(unique(used_names)))
end

"""
Resolve fitted parameters for identifiability checking.
"""
function si_funcs_to_check(ode, fitpar_names::Vector{String})
    lookup = si_parameter_lookup(ode)
    funcs = Any[]
    missing = String[]
    for pname in fitpar_names
        if haskey(lookup, pname)
            push!(funcs, lookup[pname])
        else
            push!(missing, pname)
        end
    end
    return (funcs = funcs, missing = missing)
end

"""
Normalize identifiability result dictionaries to `String => String`.
"""
function si_result_string_map(res)
    out = Dict{String,String}()
    for (k, v) in pairs(res)
        out[string(k)] = string(v)
    end
    return out
end

"""
Run structural identifiability for one model and return tabular summary.
"""
function run_structural_identifiability_for_model(model_choice::String, settings)
    model_choice in ("model1", "model2") || error("model_choice must be 'model1' or 'model2'.")

    model_config = build_model_config(model_choice)
    fixed_bundle = load_fixed_parameters(model_config.fixedpars_file)
    fixedpars = deepcopy(fixed_bundle.values)

    fitpar_names = copy(model_config.fitpar_order)

    ode_template = model_choice == "model1" ? MODEL1_SI : MODEL2_SI
    model_note = if model_choice == "model1"
        "Exact rational ODE form with virus output yv = V."
    else
        "Rational auxiliary-state form for log10(V+1): LV' = V'/(log10_base*(V+1)), with log10_base fixed."
    end
    extra_fixed = model_choice == "model2" ? Dict("log10_base" => Float64(settings.log10_base_value)) : Dict{String,Float64}()
    fixed_subs = si_fixed_parameter_substitutions(
        ode_template,
        fixedpars;
        extra_fixed = extra_fixed,
        tol = Float64(settings.rationalize_tol),
    )
    ode_reduced = StructuralIdentifiability.set_parameter_values(ode_template, fixed_subs.subs)

    funcs_info = si_funcs_to_check(ode_reduced, fitpar_names)
    isempty(funcs_info.funcs) && error("No fitted parameters were found in the symbolic model for $(model_choice).")

    @info "Running structural identifiability" model = model_choice n_fitted = length(funcs_info.funcs) missing_fitted = funcs_info.missing n_fixed_substituted = length(fixed_subs.used_names)

    local_res = StructuralIdentifiability.assess_local_identifiability(
        ode_reduced;
        funcs_to_check = funcs_info.funcs,
        prob_threshold = Float64(settings.prob_threshold),
        loglevel = settings.loglevel,
    )
    global_res = StructuralIdentifiability.assess_identifiability(
        ode_reduced;
        funcs_to_check = funcs_info.funcs,
        prob_threshold = Float64(settings.prob_threshold),
        loglevel = settings.loglevel,
    )

    local_map = si_result_string_map(local_res)
    global_map = si_result_string_map(global_res)

    rows = Vector{NamedTuple}()
    for pname in fitpar_names
        push!(
            rows,
            (
                model = model_choice,
                parameter = pname,
                local_identifiability = get(local_map, pname, "not_checked"),
                global_identifiability = get(global_map, pname, "not_checked"),
            ),
        )
    end
    summary_df = DataFrame(rows)

    return Dict(
        "model_choice" => model_choice,
        "model_note" => model_note,
        "fitpar_names" => fitpar_names,
        "missing_fitpar_in_symbolic_model" => funcs_info.missing,
        "fixed_parameters_substituted" => fixed_subs.used_names,
        "local_identifiability" => local_map,
        "global_identifiability" => global_map,
        "summary_df" => summary_df,
    )
end

"""
Top-level workflow: run structural identifiability for selected models.

Outputs per model
-----------------
- `results/output/identifiability/<model>-structural-identifiability.csv`

Combined output
---------------
- `results/output/identifiability/structural-identifiability-results.jld2`
"""
function run_structural_identifiability_workflow(settings = default_structural_identifiability_settings())
    start_wall = now()
    t_run = time()
    mkpath(settings.output_dir)

    out = Dict{String,Any}()
    for model_choice in settings.model_choices
        res = run_structural_identifiability_for_model(String(model_choice), settings)
        out[String(model_choice)] = res
        CSV.write(
            joinpath(settings.output_dir, "$(model_choice)-structural-identifiability.csv"),
            res["summary_df"],
        )
    end

    combined_path = joinpath(settings.output_dir, "structural-identifiability-results.jld2")
    save_julia_object(combined_path, out)

    elapsed_min = round((time() - t_run) / 60; digits = 2)
    elapsed = now() - start_wall
    @info "Structural identifiability workflow finished" output_dir = settings.output_dir combined_output = combined_path elapsed_min = elapsed_min elapsed = string(elapsed)
    return out
end
