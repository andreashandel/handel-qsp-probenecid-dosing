"""
Objective function helpers
==========================

This file centralizes objective-related logic so fitting, diagnostics, and app-
like inspections all rely on the same math.

Core pieces
-----------
1. `prepare_fast_objective_data`: precomputes per-scenario vectors to reduce
   repeated joins inside optimization loops.
2. `evaluate_objective_fast`: computes weighted Gaussian NLL objective directly
   from simulated trajectories and precomputed data.
3. `compute_objective_components`: slower but more explicit residual breakdown
   for diagnostics and tables.
"""

"""
Map scenario labels to numeric doses.
"""
function scenario_to_dose_map()
    return Dict(
        "NoTreatment" => 0.0,
        "PanCytoVir10mg" => 10.0,
        "PanCytoVir100mg" => 100.0,
    )
end

"""
Extract sigma parameters from fitted + fixed dictionaries.

When duplicate names exist, fitted values take precedence.
"""
function collect_sigma_pool(fitpars::Dict{String,Float64}, fixedpars::Dict{String,Float64})
    out = Dict{String,Float64}()
    for (k, v) in fixedpars
        if startswith(k, "sigma_add_") || startswith(k, "sigma_prop_")
            out[k] = v
        end
    end
    for (k, v) in fitpars
        if startswith(k, "sigma_add_") || startswith(k, "sigma_prop_")
            out[k] = v
        end
    end
    return out
end

"""
Precompute fast objective lookup structures.

`weight_mode` options:
- `"per_block"`: each scenario x quantity block weighted by `1/n_block`
- `"equal_quantity"`: each quantity weighted by `1/n_quantity`
"""
function prepare_fast_objective_data(fitdata::DataFrame; time_round::Int = 8, min_variance::Float64 = 1e-12, weight_mode::String = "per_block")
    weight_mode in ("per_block", "equal_quantity") || error("weight_mode must be 'per_block' or 'equal_quantity'")

    data = deepcopy(fitdata)
    data.Day = round.(Float64.(data.Day), digits = time_round)

    scenario_levels = ["NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg"]
    quantity_levels = [virus_quantity_name, "IL6", "WeightLossPerc"]
    quantity_id_map = Dict(q => i for (i, q) in enumerate(quantity_levels))

    # Precompute weights.
    weights = zeros(Float64, nrow(data))
    if weight_mode == "per_block"
        grouped = combine(groupby(data, [:Scenario, :Quantity]), nrow => :n)
        wmap = Dict((String(r.Scenario), String(r.Quantity)) => 1.0 / Float64(r.n) for r in eachrow(grouped))
        for i in 1:nrow(data)
            weights[i] = wmap[(String(data.Scenario[i]), String(data.Quantity[i]))]
        end
    else
        grouped = combine(groupby(data, :Quantity), nrow => :n)
        wmap = Dict(String(r.Quantity) => 1.0 / Float64(r.n) for r in eachrow(grouped))
        for i in 1:nrow(data)
            weights[i] = wmap[String(data.Quantity[i])]
        end
    end

    # Build scenario-specific vectors used in hot objective loop.
    data_by_scenario = Dict{String,NamedTuple}()
    for s in scenario_levels
        idx = findall(data.Scenario .== s)
        if isempty(idx)
            continue
        end
        qid = [quantity_id_map[String(data.Quantity[i])] for i in idx]
        data_by_scenario[s] = (
            time = Float64.(data.Day[idx]),
            quantity_id = qid,
            value = Float64.(data.Value[idx]),
            weight = Float64.(weights[idx]),
        )
    end

    return (
        data_by_scenario = data_by_scenario,
        scenario_levels = scenario_levels,
        quantity_levels = quantity_levels,
        quantity_id_map = quantity_id_map,
        min_variance = min_variance,
        weight_mode = weight_mode,
    )
end

"""
Evaluate weighted Gaussian NLL objective using precomputed objective data.

This is the high-performance objective path used during fitting.
"""
function evaluate_objective_fast(
    model_choice::String,
    fitpars::Dict{String,Float64},
    fixedpars::Dict{String,Float64},
    y0::Dict{String,Float64},
    doses::Vector{Float64},
    scenarios::Vector{String},
    solver_settings::NamedTuple,
    objective_data;
    return_block_breakdown::Bool = false,
)
    fitpars_ode = Dict(k => v for (k, v) in fitpars if !startswith(k, "sigma_"))
    fixedpars_ode = Dict(k => v for (k, v) in fixedpars if !startswith(k, "sigma_"))
    sigma_pool = collect_sigma_pool(fitpars, fixedpars)

    quantity_levels = objective_data.quantity_levels
    qmap = objective_data.quantity_id_map
    sigma_add_by_id = zeros(Float64, length(quantity_levels))
    sigma_prop_by_id = zeros(Float64, length(quantity_levels))

    for (i, qname) in enumerate(quantity_levels)
        sigma_add_by_id[i] = get(sigma_pool, "sigma_add_$(qname)", 0.0)
        sigma_prop_by_id[i] = get(sigma_pool, "sigma_prop_$(qname)", 0.0)
    end

    virus_id = qmap[virus_quantity_name]
    il6_id = qmap["IL6"]
    wl_id = qmap["WeightLossPerc"]

    objective_total = 0.0
    block_terms = Dict{Tuple{String,String},Float64}()

    for i in eachindex(doses)
        scenario = scenarios[i]
        haskey(objective_data.data_by_scenario, scenario) || continue
        scenario_data = objective_data.data_by_scenario[scenario]

        obs_times = sort(unique(vcat([0.0], scenario_data.time)))
        allpars = merged_parameters(fitpars_ode, fixedpars_ode)
        merge!(allpars, y0)
        allpars["Ad0"] = doses[i]
        allpars["txstart"] = 1.0
        allpars["txinterval"] = 0.5
        allpars["txend"] = 3.9
        allpars["tstart"] = 0.0
        allpars["tfinal"] = solver_settings.tfinal
        allpars["dt"] = solver_settings.dt
        allpars["tols"] = solver_settings.tols

        # Keep `allpars` strictly numeric; pass non-numeric solver metadata by keyword.
        sim = simulate_model(
            model_choice,
            allpars;
            times = obs_times,
            solvertype = String(solver_settings.solvertype),
        )
        if !sim.ok
            return return_block_breakdown ? (objective = Inf, breakdown_nll = DataFrame()) : Inf
        end

        sdf = sim.df
        target_time = maximum(obs_times)
        if nrow(sdf) == 0 || maximum(sdf.time) < (target_time - solver_settings.dt / 2)
            return return_block_breakdown ? (objective = Inf, breakdown_nll = DataFrame()) : Inf
        end

        time_to_idx = Dict(Float64(t) => idx for (idx, t) in enumerate(Float64.(sdf.time)))
        qid = scenario_data.quantity_id
        pred = zeros(Float64, length(qid))

        for j in eachindex(qid)
            t = Float64(scenario_data.time[j])
            idx = get(time_to_idx, t, 0)
            idx == 0 && return(return_block_breakdown ? (objective = Inf, breakdown_nll = DataFrame()) : Inf)
            if qid[j] == virus_id
                pred[j] = transform_virus(sdf.V[idx])
            elseif qid[j] == il6_id
                pred[j] = Float64(sdf.F[idx])
            elseif qid[j] == wl_id
                pred[j] = Float64(sdf.S[idx])
            end
        end

        add = sigma_add_by_id[qid]
        prop = sigma_prop_by_id[qid]
        variance = max.(add .^ 2 .+ (prop .* pred) .^ 2, objective_data.min_variance)
        residual = scenario_data.value .- pred
        nll_point = 0.5 .* (log.(variance) .+ (residual .^ 2) ./ variance)

        any(!isfinite, nll_point) && return(return_block_breakdown ? (objective = Inf, breakdown_nll = DataFrame()) : Inf)

        w = scenario_data.weight
        objective_total += sum(w .* nll_point)

        if return_block_breakdown
            for j in eachindex(qid)
                qname = quantity_levels[qid[j]]
                key = (qname, scenario)
                block_terms[key] = get(block_terms, key, 0.0) + w[j] * nll_point[j]
            end
        end
    end

    if !return_block_breakdown
        return objective_total
    end

    rows = Vector{NamedTuple{(:Quantity, :Scenario, :weighted_nll),Tuple{String,String,Float64}}}()
    for (key, val) in block_terms
        push!(rows, (Quantity = key[1], Scenario = key[2], weighted_nll = val))
    end
    return (objective = objective_total, breakdown_nll = DataFrame(rows))
end

"""
Build long-format predictions from simulation output.

This helper is mainly used for plotting diagnostics.
"""
function build_prediction_long(sim_df::DataFrame)
    required = ["time", "Scenario", "V", "F", "S"]
    all(in.(required, Ref(String.(names(sim_df))))) || error("sim_df missing required columns")

    out = DataFrame(Scenario = String[], Day = Float64[], Quantity = String[], Predicted = Float64[])
    for row in eachrow(sim_df)
        push!(out, (Scenario = String(row.Scenario), Day = Float64(row.time), Quantity = virus_quantity_name, Predicted = transform_virus(Float64(row.V))))
        push!(out, (Scenario = String(row.Scenario), Day = Float64(row.time), Quantity = "IL6", Predicted = Float64(row.F)))
        push!(out, (Scenario = String(row.Scenario), Day = Float64(row.time), Quantity = "WeightLossPerc", Predicted = Float64(row.S)))
    end
    return out
end

"""
Compute full residual/objective breakdown tables for diagnostics.

Compared to `evaluate_objective_fast`, this function prioritizes readability and
explicit tabular output over absolute speed.
"""
function compute_objective_components(
    fitdata::DataFrame,
    pred_long::DataFrame,
    sigma_pool::Dict{String,Float64};
    min_variance::Float64 = 1e-12,
    weight_mode::String = "per_block",
)
    weight_mode in ("per_block", "equal_quantity") || error("weight_mode must be 'per_block' or 'equal_quantity'")

    f = deepcopy(fitdata)
    p = deepcopy(pred_long)

    f.Day = Float64.(f.Day)
    p.Day = Float64.(p.Day)

    joined = innerjoin(
        select(f, [:Scenario, :Quantity, :Day, :Value]),
        select(p, [:Scenario, :Quantity, :Day, :Predicted]),
        on = [:Scenario, :Quantity, :Day],
    )
    if nrow(joined) == 0
        return (objective = Inf, breakdown_nll = DataFrame(), breakdown_ssr = DataFrame(), residuals = DataFrame())
    end

    if weight_mode == "per_block"
        counts = combine(groupby(f, [:Scenario, :Quantity]), nrow => :n)
        counts.weight = 1.0 ./ Float64.(counts.n)
        joined = leftjoin(joined, select(counts, [:Scenario, :Quantity, :weight]), on = [:Scenario, :Quantity])
    else
        counts = combine(groupby(f, :Quantity), nrow => :n)
        counts.weight = 1.0 ./ Float64.(counts.n)
        joined = leftjoin(joined, select(counts, [:Quantity, :weight]), on = :Quantity)
    end

    adds = [get(sigma_pool, "sigma_add_$(q)", 0.0) for q in joined.Quantity]
    props = [get(sigma_pool, "sigma_prop_$(q)", 0.0) for q in joined.Quantity]

    joined.add = adds
    joined.prop = props
    joined.Variance = max.(joined.add .^ 2 .+ (joined.prop .* joined.Predicted) .^ 2, min_variance)
    joined.Residual = joined.Value .- joined.Predicted
    joined.StdResid = joined.Residual ./ sqrt.(joined.Variance)
    joined.NLLPoint = 0.5 .* (log.(joined.Variance) .+ (joined.Residual .^ 2) ./ joined.Variance)
    joined.NLLWeightedResid = sqrt.(joined.weight) .* joined.StdResid

    any(!isfinite, joined.NLLPoint) && return (objective = Inf, breakdown_nll = DataFrame(), breakdown_ssr = DataFrame(), residuals = joined)

    g_nll = combine(groupby(joined, [:Quantity, :Scenario, :weight]), :NLLPoint => sum => :nll)
    g_nll.weighted_nll = g_nll.weight .* g_nll.nll

    g_ssr = combine(groupby(joined, [:Quantity, :Scenario, :weight]), :Residual => (x -> sum(x .^ 2)) => :ssr)
    g_ssr.weighted_ssr = g_ssr.weight .* g_ssr.ssr

    return (
        objective = sum(g_nll.weighted_nll),
        breakdown_nll = g_nll,
        breakdown_ssr = g_ssr,
        residuals = joined,
    )
end
