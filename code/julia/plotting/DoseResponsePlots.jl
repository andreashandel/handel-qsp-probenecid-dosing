module DoseResponsePlots

using CairoMakie
using Colors
using DataFrames
using Statistics

export plot_outcomes, save_dose_response_figures

const BASE_PALETTE = [colorant"#0072B2", colorant"#009E73", colorant"#D55E00"]
const LINESTYLES = Dict(1 => :solid, 2 => :dash, 3 => :dot)

function _stack_reduction(dose_response_list::Vector{Dict})
    frames = DataFrame[]
    for (idx, res) in enumerate(dose_response_list)
        df = copy(res[:reduction_df])
        df.rep = fill(idx, nrow(df))
        push!(frames, df)
    end
    return reduce(vcat, frames)
end

"""
    plot_outcomes(dose_response_list; scenarios) -> Figure

Create the multi-panel dose-response figure that mirrors
`code/plotting-code/dose-response-plot-function.R`.  The input is the list of
simulation dictionaries produced by `run_dose_predictions` and the `scenarios`
argument selects which scheduling labels to visualise.
"""
function plot_outcomes(dose_response_list::Vector{Dict};
    scenarios::Vector{String})
    metrics = [:perc_AUCV, :perc_AUCF, :perc_AUCS]
    ylabels = Dict(
        :perc_AUCV => "Log Viral Load Reduction (%)",
        :perc_AUCF => "Innate Response Reduction (%)",
        :perc_AUCS => "Morbidity Reduction (%)",
    )

    df_all = _stack_reduction(dose_response_list)
    df_all = filter(:Scenario => in(scenarios), df_all)
    df_all.Scenario = categorical(df_all.Scenario; ordered = true, levels = scenarios)

    baseline_df = copy(dose_response_list[1][:reduction_df])
    baseline_df = filter(:Scenario => in(scenarios), baseline_df)

    long = stack(df_all, metrics; variable_name = :metric, value_name = :value)
    baseline_long = stack(baseline_df, metrics; variable_name = :metric, value_name = :baseline)

    summ = combine(groupby(long, [:Scenario, :Dose, :metric]),
        :value => mean => :mean,
        :value => x -> quantile(skipmissing(x), 0.025) => :lower,
        :value => x -> quantile(skipmissing(x), 0.975) => :upper)

    baseline_long = rename!(baseline_long, :value => :baseline)
    summ = leftjoin(summ, baseline_long, on = [:Scenario, :Dose, :metric])

    fig = Figure(resolution = (1200, 400))
    legend_entries = []
    axes = Axis[]

    for (idx, metric) in enumerate(metrics)
        ax = Axis(fig[1, idx];
            xlabel = idx == 2 ? "Dose" : "",
            ylabel = ylabels[metric],
            xscale = log10,
            xticks = ([1e-3, 1e-1, 1e1, 1e3], string.([1e-3, 1e-1, 1e1, 1e3])),
            yticksvisible = true)
        push!(axes, ax)

        for (j, scen) in enumerate(scenarios)
            colour = BASE_PALETTE[j]
            ls = LINESTYLES[1 + (j - 1) % length(LINESTYLES)]
            dfm = filter(row -> row.metric == metric && row.Scenario == scen, summ)
            if nrow(dfm) == 0
                continue
            end
            if idx == 1
                push!(legend_entries, (LineElement(color = colour, linestyle = ls, linewidth = 3), scen))
            end
            if :lower in names(dfm)
                band!(ax, dfm.Dose, dfm.lower, dfm.upper; color = (colour, 0.2), transparency = true)
            end
            if :baseline in names(dfm)
                lines!(ax, dfm.Dose, dfm.baseline; color = colour, linestyle = ls, linewidth = 3)
            end
        end
        vlines!(ax, [10, 100]; linestyle = :dash, color = :black)
    end

    legend = Legend(fig[0, :], legend_entries; orientation = :horizontal, tellwidth = false)
    legend.padding = 5
    for ax in axes
        ax.topspinevisible = false
        ax.rightspinevisible = false
    end

    return fig
end

"""
    save_dose_response_figures(dose_response_list, output_path; scenarios)

Convenience wrapper that calls [`plot_outcomes`] and writes the figure to
`output_path`.
"""
function save_dose_response_figures(dose_response_list::Vector{Dict},
    output_path::AbstractString; scenarios::Vector{String})
    fig = plot_outcomes(dose_response_list; scenarios = scenarios)
    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

end # module
