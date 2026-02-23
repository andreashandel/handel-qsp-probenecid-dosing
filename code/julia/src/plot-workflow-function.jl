"""
Plotting workflow (Julia equivalent of `run-plotting.R`)
=========================================================

This file turns fitted and dose-response outputs into publication-style figures
and parameter tables.

Scope
-----
To keep the first Julia implementation maintainable, plotting is intentionally
straightforward and explicit:

- Time-series figures for schedules.
- Best-fit overlays for baseline doses.
- Diagnostic residual plots.
- Dose-response curves.
- Parameter tables as CSV.
"""

"""
Default settings for plotting workflow.
"""
function default_plot_settings(model_choice::String = "model1")
    return (
        model_choice = model_choice,
        nsamp = 1,
        make_timeseries_figures = true,
        make_diagnostic_figures = true,
        make_dose_response_figures = true,
        make_parameter_tables = true,
        output_fig_dir = repo_path("results", "figures", "julia"),
        output_table_dir = repo_path("results", "tables"),
        dose_x_limits_baseline = (1e-2, 1e5),
        dose_x_limits_txstart = (1e-2, 1e5),
        dose_x_limits_txinterval = (1e-2, 1e5),
        dose_x_breaks = [1e-3, 1e-1, 1e1, 1e3],
        dose_y_limits = (0.0, 100.0),
        results_file = repo_path("results", "output", "$(model_choice)-dose-response-results.jld2"),
        bestfit_file = repo_path("results", "output", "$(model_choice)-bestfit-sample.jld2"),
        fallback_bestfit_file = repo_path("results", "output", "$(model_choice)-bestfit-multistart.jld2"),
    )
end

"""
Compute axis limits following the R plotting workflow formulas.
"""
function compute_r_style_limits(vals::Vector{Float64}; logy::Bool = false, clamp_zero::Bool = false)
    finite_vals = [v for v in vals if isfinite(v)]
    isempty(finite_vals) && return nothing

    if logy
        pos_vals = [v for v in finite_vals if v > 0.0]
        isempty(pos_vals) && return (1e-8, 1.0)

        min_val = minimum(pos_vals)
        max_val = maximum(pos_vals)
        if min_val == max_val
            lower = min_val / 2
            upper = max_val * 2
        else
            lower = min_val / 1.5
            upper = max_val * 1.25
        end
        lower = max(lower, min_val / 10)
        lower = max(lower, 1e-8)
        upper = max(upper, lower * 1.1)
        return (lower, upper)
    end

    min_val = minimum(finite_vals)
    max_val = maximum(finite_vals)
    if min_val == max_val
        pad = min_val == 0 ? 1.0 : abs(min_val) * 0.1
    else
        pad = (max_val - min_val) * 0.1
    end
    if !isfinite(pad) || pad == 0
        pad = 1.0
    end
    lower = min_val - pad
    upper = max_val + pad
    if clamp_zero && min_val >= 0
        lower = 0.0
    end
    if lower == upper
        lower -= 1.0
        upper += 1.0
    end
    return (lower, upper)
end

"""
Create a 3x3 time-series panel matching the R plotting workflow.

Panels (left to right, top to bottom):
- U, I, V
- F, A, S
- Ad, Ac, At

If `obs_data` is provided, observed points are overlaid for V/F/S.
"""
function plot_timeseries_panel(df::DataFrame; title::String, outfile::String, obs_data::Union{Nothing,DataFrame} = nothing, tmax::Float64 = 7.0)
    specs = [
        (col = "U", label = "Uninfected Cells", logy = true, quantity = nothing),
        (col = "I", label = "Infected Cells", logy = true, quantity = nothing),
        (col = "V", label = "Virus Load", logy = true, quantity = virus_quantity_name),
        (col = "F", label = "Innate Response (IL6)", logy = false, quantity = "IL6"),
        (col = "A", label = "Adaptive Response", logy = true, quantity = nothing),
        (col = "S", label = "Morbidity (weight loss)", logy = false, quantity = "WeightLossPerc"),
        (col = "Ad", label = "Drug depot", logy = true, quantity = nothing),
        (col = "Ac", label = "Drug central", logy = true, quantity = nothing),
        (col = "At", label = "Drug target", logy = true, quantity = nothing),
    ]

    format_dose_label(d::Real) = isapprox(Float64(d), 0.0; atol = 1e-12) ? "no drug" : "$(Float64(d)) mg/kg"
    safe_log(v::AbstractVector{<:Real}) = [x > 0 ? Float64(x) : NaN for x in v]

    dose_vals = sort(unique(Float64.(df.Dose)))
    palette = [:black, :dodgerblue, :orange, :forestgreen, :firebrick, :purple, :brown, :goldenrod]
    dose_colors = Dict(d => palette[mod1(i, length(palette))] for (i, d) in enumerate(dose_vals))

    # Prepare observed-data table only when provided and non-empty.
    obs_prepped = nothing
    if !isnothing(obs_data) && nrow(obs_data) > 0
        odf = deepcopy(obs_data)
        if "Dose" ∉ names(odf) && "Scenario" ∈ names(odf)
            scenario_to_dose = Dict("NoTreatment" => 0.0, "PanCytoVir10mg" => 10.0, "PanCytoVir100mg" => 100.0)
            odf.Dose = [get(scenario_to_dose, String(s), NaN) for s in odf.Scenario]
        end
        odf = odf[in.(Float64.(odf.Dose), Ref(dose_vals)), :]
        odf = odf[in.(String.(odf.Quantity), Ref([virus_quantity_name, "IL6", "WeightLossPerc"])), :]
        if nrow(odf) > 0
            vals = Float64.(odf.Value)
            q = String.(odf.Quantity)
            odf.ValuePlot = [q[i] == virus_quantity_name ? inverse_transform_virus(vals[i]) : vals[i] for i in eachindex(vals)]
            odf.Day = Float64.(odf.Day)
            odf.Dose = Float64.(odf.Dose)
            obs_prepped = odf
        end
    end

    fig = Figure(size = (1200, 850))
    Label(fig[0, 1:3], title, fontsize = 20)

    legend_axis = nothing
    for (idx, spec) in enumerate(specs)
        row = fld(idx - 1, 3) + 1
        col = mod1(idx, 3)
        yscale_mode = spec.logy ? log10 : identity
        ax = Axis(
            fig[row, col],
            title = spec.label,
            xlabel = row == 3 ? "Time (days)" : "",
            yscale = yscale_mode,
        )

        if legend_axis === nothing
            legend_axis = ax
        end

        if spec.col ∉ names(df)
            text!(ax, 0.5, 0.5, text = "Missing column: $(spec.col)", space = :relative)
            xlims!(ax, 0.0, tmax)
            if spec.logy
                ylims!(ax, 1e-8, 1.0)
            end
            continue
        end

        panel_vals = Float64[]
        for d in dose_vals
            sub = df[df.Dose .== d, :]
            isempty(sub) && continue
            yvals = Float64.(sub[!, spec.col])
            yplot = spec.logy ? safe_log(yvals) : yvals
            lines!(ax, Float64.(sub.time), yplot, color = dose_colors[d], linewidth = 2, label = format_dose_label(d))
            append!(panel_vals, [y for y in yvals if isfinite(y)])
        end

        # Overlay observed data for fitted quantities only (V/F/S).
        if !isnothing(obs_prepped) && !isnothing(spec.quantity)
            od = obs_prepped[String.(obs_prepped.Quantity) .== spec.quantity, :]
            for d in dose_vals
                dsub = od[od.Dose .== d, :]
                isempty(dsub) && continue
                yvals = Float64.(dsub.ValuePlot)
                yplot = spec.logy ? safe_log(yvals) : yvals
                scatter!(ax, dsub.Day, yplot, color = dose_colors[d], markersize = 8)
                append!(panel_vals, [y for y in yvals if isfinite(y)])
            end
        end

        xlims!(ax, 0.0, tmax)
        clamp_zero = spec.col in ("F", "S")
        lims = compute_r_style_limits(panel_vals; logy = spec.logy, clamp_zero = clamp_zero)
        if !isnothing(lims)
            ylims!(ax, lims[1], lims[2])
        end
    end

    !isnothing(legend_axis) && axislegend(legend_axis, position = :rb)
    save(outfile, fig)
end

"""
Build diagnostic residual plots for one sample.
"""
function plot_diagnostics(bestfit::Dict{String,Any}, sim_df::DataFrame, outfile_prefix::String)
    dose_levels = [0.0, 10.0, 100.0]
    scenario_map = Dict(0.0 => "NoTreatment", 10.0 => "PanCytoVir10mg", 100.0 => "PanCytoVir100mg")

    sdf = sim_df[(sim_df.Schedule .== "s1") .& in.(Float64.(sim_df.Dose), Ref(dose_levels)), :]
    sdf = deepcopy(sdf)
    sdf.Scenario = [scenario_map[Float64(d)] for d in sdf.Dose]

    pred_long = build_prediction_long(sdf)
    fitpars = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["fitpars"])
    fixedpars = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["fixedpars"])
    sigma_pool = collect_sigma_pool(fitpars, fixedpars)
    fitdata = bestfit["fitdata"]
    components = compute_objective_components(fitdata, pred_long, sigma_pool; weight_mode = "equal_quantity")
    resid = components.residuals
    nrow(resid) > 0 || return

    # Residual vs day.
    fig1 = Figure(size = (900, 500))
    ax = Axis(fig1[1, 1], xlabel = "Day", ylabel = "NLL-weighted standardized residual", title = "Residual diagnostics")
    for q in unique(String.(resid.Quantity))
        sub = resid[resid.Quantity .== q, :]
        scatter!(ax, sub.Day, sub.NLLWeightedResid, label = q, markersize = 8)
    end
    hlines!(ax, [0.0], color = :black, linestyle = :dash)
    resid_abs_max = maximum(abs.(Float64.(resid.NLLWeightedResid)))
    if !isfinite(resid_abs_max) || resid_abs_max == 0.0
        resid_abs_max = 1.0
    end
    resid_pad = 0.1 * resid_abs_max
    y_lim_comb = resid_abs_max + resid_pad
    ylims!(ax, -y_lim_comb, y_lim_comb)
    axislegend(ax)
    save("$(outfile_prefix)-residuals.png", fig1)

    # Predicted vs observed.
    fig2 = Figure(size = (700, 600))
    ax2 = Axis(fig2[1, 1], xlabel = "Observed", ylabel = "Predicted", title = "Predicted vs observed")
    for q in unique(String.(resid.Quantity))
        sub = resid[resid.Quantity .== q, :]
        scatter!(ax2, sub.Value, sub.Predicted, label = q, markersize = 8)
    end
    mins = minimum(vcat(Float64.(resid.Value), Float64.(resid.Predicted)))
    maxs = maximum(vcat(Float64.(resid.Value), Float64.(resid.Predicted)))
    pad = 0.05 * (maxs - mins + ((maxs - mins) == 0.0 ? 1.0 : 0.0))
    low = mins - pad
    high = maxs + pad
    lines!(ax2, [mins, maxs], [mins, maxs], color = :black, linestyle = :dash)
    xlims!(ax2, low, high)
    ylims!(ax2, low, high)
    axislegend(ax2)
    save("$(outfile_prefix)-gof.png", fig2)
end

"""
Plot dose-response summaries for selected schedules.
"""
function plot_dose_response(
    results_list,
    schedules::Vector{String},
    outfile::String;
    x_limits::Tuple{Float64,Float64} = (1e-2, 1e5),
    x_breaks::Vector{Float64} = [1e-3, 1e-1, 1e1, 1e3],
    y_limits::Tuple{Float64,Float64} = (0.0, 100.0),
)
    fig = Figure(size = (1200, 400))
    ax1 = Axis(fig[1, 1], xlabel = "Dose", ylabel = "Log Viral Load Reduction (%)", xscale = log10)
    ax2 = Axis(fig[1, 2], xlabel = "Dose", ylabel = "Innate Response Reduction (%)", xscale = log10)
    ax3 = Axis(fig[1, 3], xlabel = "Dose", ylabel = "Morbidity Reduction (%)", xscale = log10)

    schedule_labels = Dict(
        "s1" => "baseline",
        "s2" => "d2 start",
        "s3" => "d3 start",
        "s4" => "daily tx",
        "s5" => "single tx",
    )
    base_palette = ["#0072B2", "#009E73", "#D55E00"]
    line_styles = [:solid, :dash, :dot]
    colors = [base_palette[mod1(i, length(base_palette))] for i in eachindex(schedules)]
    styles = [line_styles[mod1(i, length(line_styles))] for i in eachindex(schedules)]
    color_map = Dict(String(schedules[i]) => colors[i] for i in eachindex(schedules))
    style_map = Dict(String(schedules[i]) => styles[i] for i in eachindex(schedules))

    perc_cols = [:perc_AUCV, :perc_AUCF, :perc_AUCS]
    use_perc = all(all(c in names(sample["all_results_df"]) for c in perc_cols) for sample in results_list)
    metric_specs = [
        (ax = ax1, col = use_perc ? :perc_AUCV : :AUCV),
        (ax = ax2, col = use_perc ? :perc_AUCF : :AUCF),
        (ax = ax3, col = use_perc ? :perc_AUCS : :AUCS),
    ]

    for sched in String.(schedules)
        c = color_map[sched]
        ls = style_map[sched]
        label = get(schedule_labels, sched, sched)

        # Thin overlays for each sample (R sample_display = "lines"-like).
        for sample in results_list
            df = sample["all_results_df"]
            sub = df[df.Schedule .== sched, :]
            sub = sub[(sub.Dose .> 0.0) .& isfinite.(Float64.(sub.Dose)), :]
            isempty(sub) && continue
            sub = sort(sub, :Dose)

            for spec in metric_specs
                y = Float64.(sub[!, spec.col])
                x = Float64.(sub.Dose)
                mask = isfinite.(x) .& isfinite.(y)
                any(mask) || continue
                lines!(spec.ax, x[mask], y[mask], color = c, linestyle = ls, linewidth = 0.8, alpha = 0.45)
            end
        end

        # Thick baseline curve from first sample.
        if !isempty(results_list)
            base_df = results_list[1]["all_results_df"]
            base_sub = base_df[base_df.Schedule .== sched, :]
            base_sub = base_sub[(base_sub.Dose .> 0.0) .& isfinite.(Float64.(base_sub.Dose)), :]
            if !isempty(base_sub)
                base_sub = sort(base_sub, :Dose)
                for (j, spec) in enumerate(metric_specs)
                    y = Float64.(base_sub[!, spec.col])
                    x = Float64.(base_sub.Dose)
                    mask = isfinite.(x) .& isfinite.(y)
                    any(mask) || continue
                    if j == 1
                        lines!(spec.ax, x[mask], y[mask], color = c, linestyle = ls, linewidth = 2.0, label = label)
                    else
                        lines!(spec.ax, x[mask], y[mask], color = c, linestyle = ls, linewidth = 2.0)
                    end
                end
            end
        end
    end

    valid_breaks = [b for b in x_breaks if x_limits[1] <= b <= x_limits[2]]
    if !isempty(valid_breaks)
        tick_labels = [string(b) for b in valid_breaks]
        ax1.xticks = (valid_breaks, tick_labels)
        ax2.xticks = (valid_breaks, tick_labels)
        ax3.xticks = (valid_breaks, tick_labels)
    end
    xlims!(ax1, x_limits[1], x_limits[2])
    xlims!(ax2, x_limits[1], x_limits[2])
    xlims!(ax3, x_limits[1], x_limits[2])
    ylims!(ax1, y_limits[1], y_limits[2])
    ylims!(ax2, y_limits[1], y_limits[2])
    ylims!(ax3, y_limits[1], y_limits[2])
    vlines!(ax1, [10.0, 100.0], color = :black, linestyle = :dash)
    vlines!(ax2, [10.0, 100.0], color = :black, linestyle = :dash)
    vlines!(ax3, [10.0, 100.0], color = :black, linestyle = :dash)
    Legend(fig[0, 1:3], ax1, orientation = :horizontal)

    save(outfile, fig)
end

"""
Write parameter table for one sample as CSV.
"""
function write_parameter_table(bestfit::Dict{String,Any}, outfile::String)
    namesv = String.(bestfit["fitparnames"])
    fitpars = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["fitpars"])
    parlabels = Dict{String,String}(String(k) => String(v) for (k, v) in bestfit["parlabels"])
    vals = [Float64(fitpars[n]) for n in namesv]
    labels = [get(parlabels, n, n) for n in namesv]
    df = DataFrame(Parameter = namesv, Value = vals, Label = labels)
    CSV.write(outfile, df)
end

"""
Write parameter table figure (PNG) for one sample.
"""
function write_parameter_table_figure(bestfit::Dict{String,Any}, outfile::String)
    namesv = String.(bestfit["fitparnames"])
    fitpars = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in bestfit["fitpars"])
    parlabels = Dict{String,String}(String(k) => String(v) for (k, v) in bestfit["parlabels"])
    vals = [Float64(fitpars[n]) for n in namesv]
    labels = [get(parlabels, n, n) for n in namesv]

    nrows_tbl = length(namesv)
    fig_h = max(260, 80 + 26 * (nrows_tbl + 1))
    fig = Figure(size = (1200, fig_h))
    gl = fig[1, 1] = GridLayout()

    Label(gl[1, 1], "Parameter", fontsize = 14, halign = :left)
    Label(gl[1, 2], "Value", fontsize = 14, halign = :right)
    Label(gl[1, 3], "Label", fontsize = 14, halign = :left)

    for i in 1:nrows_tbl
        Label(gl[i + 1, 1], namesv[i], fontsize = 12, halign = :left)
        Label(gl[i + 1, 2], @sprintf("%.6g", vals[i]), fontsize = 12, halign = :right)
        Label(gl[i + 1, 3], labels[i], fontsize = 12, halign = :left)
    end

    colgap!(gl, 18)
    rowgap!(gl, 6)
    save(outfile, fig)
end

"""
Top-level plotting workflow.
"""
function run_plot_workflow(settings = default_plot_settings())
    start_wall = now()
    t_run = time()
    @info "Starting Julia run-plotting" model = settings.model_choice

    isfile(settings.results_file) || error("Dose-response results not found: $(settings.results_file)")

    bestfit_path = isfile(settings.bestfit_file) ? settings.bestfit_file : settings.fallback_bestfit_file
    isfile(bestfit_path) || error("Bestfit file not found: $(settings.bestfit_file) (and fallback missing)")

    results_obj = load_julia_object(settings.results_file)
    simres_list = results_obj["simres_list"]
    bestfit_list = load_julia_object(bestfit_path)

    nsamp = min(settings.nsamp, length(bestfit_list), length(simres_list))
    mkpath(settings.output_fig_dir)
    mkpath(settings.output_table_dir)

    if settings.make_timeseries_figures
        @info "Generating Julia time-series figures"
        first_ts = simres_list[1]["timeseries_df"]
        for sched in ["s1", "s2", "s3", "s4", "s5"]
            sub = first_ts[first_ts.Schedule .== sched, :]
            isempty(sub) && continue
            plot_timeseries_panel(
                sub;
                title = "$(settings.model_choice) $(sched)",
                outfile = joinpath(settings.output_fig_dir, "$(settings.model_choice)-timeseries-$(sched)-julia.png"),
            )
        end

        # Baseline best-fit overlays with observed data (same 9-panel layout).
        for i in 1:nsamp
            simdf = simres_list[i]["timeseries_df"]
            sub = simdf[(simdf.Schedule .== "s1") .& in.(Float64.(simdf.Dose), Ref([0.0, 10.0, 100.0])), :]
            objective_label = "NA"
            if haskey(bestfit_list[i], "objective")
                objective_val = try
                    Float64(bestfit_list[i]["objective"])
                catch
                    NaN
                end
                if isfinite(objective_val)
                    objective_label = @sprintf("%.6g", objective_val)
                end
            end
            plot_timeseries_panel(
                sub;
                title = "$(settings.model_choice) baseline bestfit sample $(i) | objective: $(objective_label)",
                outfile = joinpath(settings.output_fig_dir, "$(settings.model_choice)-bestfit$(i)-julia.png"),
                obs_data = bestfit_list[i]["fitdata"],
            )
        end
    end

    if settings.make_diagnostic_figures
        @info "Generating Julia diagnostic figures"
        for i in 1:nsamp
            plot_diagnostics(
                bestfit_list[i],
                simres_list[i]["timeseries_df"],
                joinpath(settings.output_fig_dir, "$(settings.model_choice)-sample$(i)-julia"),
            )
        end
    end

    if settings.make_dose_response_figures
        @info "Generating Julia dose-response figures"
        plot_dose_response(
            simres_list[1:nsamp],
            ["s1"],
            joinpath(settings.output_fig_dir, "$(settings.model_choice)-dose-response-baseline-julia.png");
            x_limits = (Float64(settings.dose_x_limits_baseline[1]), Float64(settings.dose_x_limits_baseline[2])),
            x_breaks = Float64.(settings.dose_x_breaks),
            y_limits = (Float64(settings.dose_y_limits[1]), Float64(settings.dose_y_limits[2])),
        )
        plot_dose_response(
            simres_list[1:nsamp],
            ["s1", "s2", "s3"],
            joinpath(settings.output_fig_dir, "$(settings.model_choice)-dose-response-txstart-julia.png");
            x_limits = (Float64(settings.dose_x_limits_txstart[1]), Float64(settings.dose_x_limits_txstart[2])),
            x_breaks = Float64.(settings.dose_x_breaks),
            y_limits = (Float64(settings.dose_y_limits[1]), Float64(settings.dose_y_limits[2])),
        )
        plot_dose_response(
            simres_list[1:nsamp],
            ["s1", "s4", "s5"],
            joinpath(settings.output_fig_dir, "$(settings.model_choice)-dose-response-txinterval-julia.png");
            x_limits = (Float64(settings.dose_x_limits_txinterval[1]), Float64(settings.dose_x_limits_txinterval[2])),
            x_breaks = Float64.(settings.dose_x_breaks),
            y_limits = (Float64(settings.dose_y_limits[1]), Float64(settings.dose_y_limits[2])),
        )
    end

    if settings.make_parameter_tables
        @info "Writing Julia parameter tables"
        for i in 1:nsamp
            write_parameter_table(
                bestfit_list[i],
                joinpath(settings.output_table_dir, "$(settings.model_choice)-parametertab$(i)-julia.csv"),
            )
            write_parameter_table_figure(
                bestfit_list[i],
                joinpath(settings.output_fig_dir, "$(settings.model_choice)-parametertab$(i)-julia.png"),
            )
        end
    end

    elapsed_min = round((time() - t_run) / 60; digits = 2)
    elapsed = now() - start_wall
    @info "Julia run-plotting finished" elapsed_min = elapsed_min elapsed = string(elapsed)
    return nothing
end
