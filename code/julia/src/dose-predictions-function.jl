"""
Dose-prediction workflow (Julia equivalent of `run-dose-predictions.R`)
========================================================================

This module runs forward simulations over dose grids and schedule variants,
producing two classes of outputs:

1. Summary AUC dose-response tables.
2. Full time-series trajectories for a selected subset of doses.

The implementation mirrors the R workflow, but stores outputs in `.jld2`.
"""

"""
Default settings for dose-prediction runs.
"""
function default_dose_prediction_settings(model_choice::String = "model1")
    dose_grid = 10 .^ range(-2, 4; length = 50)
    ts_doses = [0.0, 1.0, 10.0, 100.0, 1_000.0, 10_000.0]
    all_doses = sort(unique(vcat(ts_doses, dose_grid)))

    schedule_defs = [
        Dict("txstart" => 1.0, "txend" => 3.9, "txinterval" => 0.5, "name" => "s1", "label" => "baseline"),
        Dict("txstart" => 2.0, "txend" => 4.9, "txinterval" => 0.5, "name" => "s2", "label" => "d2 start"),
        Dict("txstart" => 3.0, "txend" => 5.9, "txinterval" => 0.5, "name" => "s3", "label" => "d3 start"),
        Dict("txstart" => 1.0, "txend" => 3.9, "txinterval" => 1.0, "name" => "s4", "label" => "daily tx"),
        Dict("txstart" => 1.0, "txend" => 1.0, "txinterval" => 1.0, "name" => "s5", "label" => "single tx"),
    ]

    return (
        model_choice = model_choice,
        bestfit_file = repo_path("results", "output", "$(model_choice)-bestfit-sample.jld2"),
        fallback_bestfit_file = repo_path("results", "output", "$(model_choice)-bestfit-multistart.jld2"),
        timeseries_doses = ts_doses,
        all_doses = all_doses,
        schedule_defs = schedule_defs,
        solver_settings = (solvertype = "rodas5p", tols = 1e-10, dt = 0.02, tfinal = 7.0),
        output_dir = repo_path("results", "output"),
    )
end

"""
Compute percent reduction columns using Dose=0 baseline within each schedule.
"""
function add_percent_reduction(df::DataFrame)
    out = deepcopy(df)
    out.perc_AUCV = fill(NaN, nrow(out))
    out.perc_AUCF = fill(NaN, nrow(out))
    out.perc_AUCS = fill(NaN, nrow(out))

    for g in groupby(out, :Schedule)
        base_rows = findall(g.Dose .== 0.0)
        isempty(base_rows) && continue
        base_idx = base_rows[1]
        base_v = g.AUCV[base_idx]
        base_f = g.AUCF[base_idx]
        base_s = g.AUCS[base_idx]
        g.perc_AUCV .= (base_v .- g.AUCV) ./ base_v .* 100.0
        g.perc_AUCF .= (base_f .- g.AUCF) ./ base_f .* 100.0
        g.perc_AUCS .= (base_s .- g.AUCS) ./ base_s .* 100.0
    end
    return out
end

"""
Simulate one fitted sample over all doses and schedules.
"""
function simulate_dose_predictions_for_sample(bestfit::Dict{String,Any}, settings)
    fitpars = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["fitpars"])
    fixedpars = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["fixedpars"])
    y0 = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["Y0"])

    fitpars_ode = Dict(k => v for (k, v) in fitpars if !startswith(k, "sigma_"))
    fixedpars_ode = Dict(k => v for (k, v) in fixedpars if !startswith(k, "sigma_"))

    output_times = collect(0.0:settings.solver_settings.dt:settings.solver_settings.tfinal)

    all_summary_rows = Vector{NamedTuple}()
    all_ts = DataFrame()
    failures = DataFrame(Schedule = String[], Dose = Float64[], Message = String[])
    warnings = DataFrame(Schedule = String[], Dose = Float64[], Message = String[])

    for sched in settings.schedule_defs
        sname = String(sched["name"])
        for dose in settings.all_doses
            allpars = merged_parameters(fitpars_ode, fixedpars_ode)
            merge!(allpars, y0)
            allpars["Ad0"] = Float64(dose)
            allpars["txstart"] = Float64(sched["txstart"])
            allpars["txend"] = Float64(sched["txend"])
            allpars["txinterval"] = Float64(sched["txinterval"])
            allpars["tstart"] = 0.0
            allpars["tfinal"] = settings.solver_settings.tfinal
            allpars["dt"] = settings.solver_settings.dt
            allpars["tols"] = settings.solver_settings.tols

            # Keep `allpars` numeric for type stability; pass solver kind separately.
            sim = simulate_model(
                settings.model_choice,
                allpars;
                times = output_times,
                solvertype = String(settings.solver_settings.solvertype),
            )
            if !sim.ok
                push!(failures, (sname, Float64(dose), "Solver failed with retcode $(sim.retcode)"))
                continue
            end
            sdf = sim.df
            if nrow(sdf) == 0 || maximum(sdf.time) < (settings.solver_settings.tfinal - settings.solver_settings.dt / 2)
                push!(failures, (sname, Float64(dose), "Simulation returned before tfinal"))
                continue
            end

            Ct = sdf.At ./ get(allpars, "Vt", 1.0)
            fu = get(allpars, "fmax", 0.0) .* Ct ./ (get(allpars, "f50", 1.0) .+ Ct)
            Cu = fu .* Ct
            fV = get(allpars, "Emax_V", 0.0) .* Cu ./ (get(allpars, "C50_V", 1.0) .+ Cu)
            fF = get(allpars, "Emax_F", 0.0) .* Cu ./ (get(allpars, "C50_F", 1.0) .+ Cu)

            aucv = trapz_area(sdf.time, transform_virus(sdf.V))
            aucf = trapz_area(sdf.time, sdf.F)
            aucs = trapz_area(sdf.time, sdf.S)
            push!(all_summary_rows, (Dose = Float64(dose), AUCV = aucv, AUCF = aucf, AUCS = aucs, Schedule = sname))

            if Float64(dose) in settings.timeseries_doses
                ts = deepcopy(sdf)
                ts.Dose = fill(Float64(dose), nrow(ts))
                ts.Schedule = fill(sname, nrow(ts))
                ts.Cu = Cu
                ts.fV = fV
                ts.fF = fF
                all_ts = isempty(all_ts) ? ts : vcat(all_ts, ts)
            end
        end
    end

    summary_df = DataFrame(all_summary_rows)
    summary_df = add_percent_reduction(summary_df)

    return Dict(
        "all_results_df" => summary_df,
        "timeseries_df" => all_ts,
        "failures" => failures,
        "warnings" => warnings,
    )
end

"""
Top-level dose-prediction workflow.
"""
function run_dose_prediction_workflow(settings = default_dose_prediction_settings())
    start_wall = now()
    t_run = time()
    @info "Starting Julia run-dose-predictions" model = settings.model_choice

    bestfit_path = isfile(settings.bestfit_file) ? settings.bestfit_file : settings.fallback_bestfit_file
    isfile(bestfit_path) || error("Bestfit file not found: $(settings.bestfit_file) (and fallback missing)")

    bestfit_list = load_julia_object(bestfit_path)
    bestfit_list isa AbstractVector || error("Bestfit file must contain a vector of bestfit objects")

    simres_list = Vector{Dict{String,Any}}()
    for (i, bf) in enumerate(bestfit_list)
        @info "Simulating dose predictions for sample" sample = i total = length(bestfit_list)
        push!(simres_list, simulate_dose_predictions_for_sample(bf, settings))
    end

    out_obj = Dict(
        "simres_list" => simres_list,
        "ts_doses" => settings.timeseries_doses,
        "bestfit_file" => bestfit_path,
    )

    mkpath(settings.output_dir)
    output_path = joinpath(settings.output_dir, "$(settings.model_choice)-dose-response-results.jld2")
    save_julia_object(output_path, out_obj)

    # Persist warning/error table as CSV for quick inspection.
    msg_rows = DataFrame(Sample = Int[], Type = String[], Schedule = String[], Dose = Float64[], Message = String[])
    for (i, item) in enumerate(simres_list)
        fdf = item["failures"]
        for row in eachrow(fdf)
            push!(msg_rows, (i, "error", String(row.Schedule), Float64(row.Dose), String(row.Message)))
        end
        wdf = item["warnings"]
        for row in eachrow(wdf)
            push!(msg_rows, (i, "warning", String(row.Schedule), Float64(row.Dose), String(row.Message)))
        end
    end
    if nrow(msg_rows) > 0
        CSV.write(joinpath(settings.output_dir, "$(settings.model_choice)-dose-response-messages.csv"), msg_rows)
    end

    elapsed_min = round((time() - t_run) / 60; digits = 2)
    elapsed = now() - start_wall
    @info "Julia run-dose-predictions finished" path = output_path elapsed_min = elapsed_min elapsed = string(elapsed)
    return nothing
end
