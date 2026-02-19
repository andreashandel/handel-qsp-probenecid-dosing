"""
Data loading and model configuration
====================================

This file mirrors the role of several R helper files:

- `fit-data-function.R`
- `fixed-parameters-function.R`
- `sigma-settings-function.R`
- `model-config-function.R`

The key design choice is to return plain Julia containers (`DataFrame`,
`Dict`, `NamedTuple`) instead of defining many custom structs. This keeps the
code approachable for users new to Julia and makes ad-hoc debugging easy in the
REPL.
"""

"""
Load and standardize processed fitting data.

Returned fields
---------------
- `fitdata`: DataFrame with transformed `Value`, ordered scenario/quantity labels,
  and `Dose`.
- `scenarios`: scenario labels in canonical order.
- `doses`: numeric doses in ascending order.
"""
function load_fit_data(; data_path::Union{Nothing,String} = nothing, time_round::Int = 8)
    path = isnothing(data_path) ? repo_path("data", "processed-data", "processeddata.csv") : data_path
    fitdata = CSV.read(path, DataFrame)

    if :Day in names(fitdata)
        fitdata.Day = round.(Float64.(fitdata.Day), digits = time_round)
    end

    scenario_levels = ["NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg"]
    quantity_levels = [virus_quantity_name, "IL6", "WeightLossPerc"]

    observed_scenarios = sort(unique(String.(fitdata.Scenario)))
    observed_quantities = sort(unique(String.(fitdata.Quantity)))
    unknown_scenarios = setdiff(observed_scenarios, scenario_levels)
    unknown_quantities = setdiff(observed_quantities, quantity_levels)

    !isempty(unknown_scenarios) && error("Unknown Scenario values: $(join(unknown_scenarios, ", "))")
    !isempty(unknown_quantities) && error("Unknown Quantity values: $(join(unknown_quantities, ", "))")

    fitdata.Scenario = String.(fitdata.Scenario)
    fitdata.Quantity = String.(fitdata.Quantity)

    fitdata.Value = [q == virus_quantity_name ? transform_virus(v) : Float64(v) for (q, v) in zip(fitdata.Quantity, fitdata.Value)]

    scenario_to_dose = Dict(
        "NoTreatment" => 0.0,
        "PanCytoVir10mg" => 10.0,
        "PanCytoVir100mg" => 100.0,
    )
    fitdata.Dose = [scenario_to_dose[s] for s in fitdata.Scenario]
    fitdata.xvals = copy(fitdata.Day)

    # Preserve canonical order in explicit vectors used downstream.
    scenarios = scenario_levels
    doses = sort(unique(fitdata.Dose))

    return (fitdata = fitdata, scenarios = scenarios, doses = doses)
end

"""
Load fixed parameters from a CSV file.

CSV format expected
-------------------
- `parname`: short machine-readable name.
- `parnamefull`: descriptive label.
- `value`: numeric value.
"""
function load_fixed_parameters(csv_path::String)
    df = CSV.read(csv_path, DataFrame)
    df.parname = strip.(String.(df.parname))
    df.parnamefull = strip.(String.(df.parnamefull))

    values = Dict{String,Float64}(String(row.parname) => Float64(row.value) for row in eachrow(df))
    labels = Dict{String,String}(String(row.parname) => String(row.parnamefull) for row in eachrow(df))

    return (values = values, labels = labels, raw = df)
end

"""
Compute sigma defaults and split into fitted vs fixed subsets.

The approach follows the R implementation exactly:
- additive sigma = sqrt(empirical variance per quantity)
- proportional sigma = 0 by default
"""
function compute_sigma_settings(fitdata::DataFrame; sigma_to_fit::Vector{String} = String[])
    qty_levels = [virus_quantity_name, "IL6", "WeightLossPerc"]

    var_map = Dict{String,Float64}()
    for q in qty_levels
        vals = Float64.(fitdata.Value[fitdata.Quantity .== q])
        v = isempty(vals) ? 0.0 : (length(vals) > 1 ? var(vals) : 0.0)
        var_map[q] = max(0.0, ifelse(isnan(v), 0.0, v))
    end

    sigma_all = Dict(
        "sigma_add_$(virus_quantity_name)" => sqrt(var_map[virus_quantity_name]),
        "sigma_prop_$(virus_quantity_name)" => 0.0,
        "sigma_add_IL6" => sqrt(var_map["IL6"]),
        "sigma_prop_IL6" => 0.0,
        "sigma_add_WeightLossPerc" => sqrt(var_map["WeightLossPerc"]),
        "sigma_prop_WeightLossPerc" => 0.0,
    )

    sigma_fit_ini = Dict{String,Float64}()
    sigma_fixed = Dict{String,Float64}()
    for (k, v) in sigma_all
        if k in sigma_to_fit
            sigma_fit_ini[k] = v
        else
            sigma_fixed[k] = v
        end
    end

    return (sigma_all = sigma_all, sigma_fit_ini = sigma_fit_ini, sigma_fixed = sigma_fixed)
end

"""
Build model configuration for `model1` or `model2`.

Fields in return value
----------------------
- `y0`: initial conditions as dictionary.
- `par_ini_full`, `lb`, `ub`: fitted-parameter defaults and bounds.
- `parlabels_full`: labels for tables and plots.
- `fixedpars_file`: path to fixed-parameter CSV.
- `state_symbols`: explicit state ordering used by simulators.
"""
function build_model_config(model_choice::String)
    model_choice in ("model1", "model2") || error("model_choice must be 'model1' or 'model2'.")

    fitpar_order = [
        "b", "k", "p", "kF", "cV", "gF", "hV",
        "Fmax", "hF", "gS", "cS", "Emax_F", "C50_F", "C50_V",
    ]

    par_ini_full = Dict(
        "b" => 1e-6,
        "k" => 1e-5,
        "p" => 1e4,
        "kF" => 0.1,
        "cV" => 100.0,
        "gF" => 1.0,
        "hV" => 1e3,
        "Fmax" => 2.0,
        "hF" => 1.0,
        "gS" => 10.0,
        "cS" => 1.0,
        "Emax_F" => 1.0,
        "C50_F" => 1e-5,
        "C50_V" => 1e-8,
    )

    lb = Dict(
        "b" => 1e-14,
        "k" => 1e-10,
        "p" => 1.0,
        "kF" => 1e-1,
        "cV" => 0.1,
        "gF" => 1e-3,
        "hV" => 1e-2,
        "Fmax" => 0.1,
        "hF" => 1e-5,
        "gS" => 1e-3,
        "cS" => 1e-5,
        "Emax_F" => 1e-5,
        "C50_F" => 1e-17,
        "C50_V" => 1e-17,
    )

    ub = Dict(
        "b" => 1e-1,
        "k" => 1.0,
        "p" => 1e10,
        "kF" => 1e5,
        "cV" => 1e5,
        "gF" => 1e3,
        "hV" => 1e5,
        "Fmax" => 1e4,
        "hF" => 1e3,
        "gS" => 1e3,
        "cS" => 1e3,
        "Emax_F" => 1.0,
        "C50_F" => 1e2,
        "C50_V" => 1e2,
    )

    parlabels_full = Dict(
        "b" => "Virus infection rate",
        "k" => "Adaptive response clearance rate",
        "p" => "Virus production rate",
        "kF" => "Innate response suppression strength",
        "cV" => "Virus removal rate",
        "gF" => "Maximum innate response induction",
        "hV" => "Virus half-maximum innate activation",
        "Fmax" => "Maximum innate response",
        "hF" => "Adaptive response half-maximum induction",
        "gS" => "Symptom induction rate",
        "cS" => "Symptom decay rate",
        "Emax_F" => "Maximum drug effect on innate response",
        "C50_F" => "Half maximum of innate response effect",
        "C50_V" => "Half maximum of virus suppression effect",
        "sigma_add_$(virus_quantity_name)" => "Sigma of VirusLoad",
        "sigma_prop_$(virus_quantity_name)" => "Proportional sigma of VirusLoad",
        "sigma_add_IL6" => "Sigma of IL6",
        "sigma_prop_IL6" => "Proportional sigma of IL6",
        "sigma_add_WeightLossPerc" => "Sigma of WeightLossPerc",
        "sigma_prop_WeightLossPerc" => "Proportional sigma of WeightLossPerc",
    )

    if model_choice == "model1"
        y0 = Dict("Ad" => 0.0, "Ac" => 0.0, "At" => 0.0, "U" => 1e7, "I" => 0.0, "V" => 1.0, "F" => 0.0, "A" => 0.0, "S" => 0.0)
        state_symbols = [:Ad, :Ac, :At, :U, :I, :V, :F, :A, :S]
        fixedpars_file = repo_path("data", "processed-data", "model1-fixed-parameters.csv")
    else
        y0 = Dict("Ad" => 0.0, "Ac" => 0.0, "At" => 0.0, "U" => 1e7, "E" => 0.0, "I" => 0.0, "V" => 1.0, "F" => 0.0, "A" => 0.0, "S" => 0.0)
        state_symbols = [:Ad, :Ac, :At, :U, :E, :I, :V, :F, :A, :S]
        fixedpars_file = repo_path("data", "processed-data", "model2-fixed-parameters.csv")
    end

    return (
        model_choice = model_choice,
        y0 = y0,
        state_symbols = state_symbols,
        fitpar_order = fitpar_order,
        par_ini_full = par_ini_full,
        lb = lb,
        ub = ub,
        parlabels_full = parlabels_full,
        fixedpars_file = fixedpars_file,
    )
end
