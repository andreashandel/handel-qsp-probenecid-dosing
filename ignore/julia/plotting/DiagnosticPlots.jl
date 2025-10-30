module DiagnosticPlots

using CairoMakie
using Colors
using DataFrames

include("../analysis/FitModel.jl")
using .ModelFitting: FitResult

export residual_grid_plot, residual_combined_plot

const SCENARIO_LABELS = Dict(
    "NoTreatment" => "No Treatment",
    "PanCytoVir10mg" => "10 mg/kg",
    "PanCytoVir100mg" => "100 mg/kg",
)

const QUANTITY_LABELS = Dict(
    "LogVirusLoad" => "Log Virus Load",
    "IL6" => "IL-6",
    "Weight" => "Weight",
)

const PALETTE = [colorant"#0072B2", colorant"#009E73", colorant"#D55E00"]
const MARKERS = [:circle, :utriangle, :rect]

function _compute_residuals(bestfit::FitResult)
    dat = bestfit.fitdata
    sim = bestfit.simresult
    times_to_keep = unique(dat.xvals)
    sim = filter(:time => t -> t in times_to_keep, sim)

    sim_long = select(sim, :time, :Scenario, :V, :F, :S)
    sim_long.V = log10.(sim_long.V .+ 1e-12)
    rename!(sim_long, :time => :Day, :V => :LogVirusLoad, :F => :IL6, :S => :Weight)
    sim_long = stack(sim_long, [:LogVirusLoad, :IL6, :Weight];
        variable_name = :Quantity, value_name = :Predicted)

    resid_df = innerjoin(dat, sim_long,
        on = [:Scenario, :Day => :Day, :Quantity])
    resid_df.Residual = resid_df.Value .- resid_df.Predicted
    return resid_df
end

"""
    residual_grid_plot(bestfit) -> Figure

Replicates the 3×3 residual diagnostic figure produced by
`code/plotting-code/make-diagnostic-figures.R`.
"""
function residual_grid_plot(bestfit::FitResult)
    resid_df = _compute_residuals(bestfit)
    scenarios = ["NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg"]
    quantities = ["LogVirusLoad", "IL6", "Weight"]

    fig = Figure(resolution = (900, 900))

    for (i, qty) in enumerate(quantities)
        for (j, scen) in enumerate(scenarios)
            ax = Axis(fig[i, j];
                xlabel = i == length(quantities) ? "Time (days)" : "",
                ylabel = j == 1 ? QUANTITY_LABELS[qty] : "",
                title = i == 1 ? SCENARIO_LABELS[scen] : "",
                xticks = collect(0:floor(Int, maximum(resid_df.Day))))
            hlines!(ax, [0.0]; linestyle = :dash, color = :gray)
            sub = filter(row -> row.Quantity == qty && row.Scenario == scen, resid_df)
            scatter!(ax, sub.Day, sub.Residual; color = PALETTE[j], markersize = 8)
        end
    end

    return fig
end

"""
    residual_combined_plot(bestfit) -> Figure

Create the overlaid residual plot analogous to
`code/plotting-code/make-combined-diagnostic-figures.R`.
"""
function residual_combined_plot(bestfit::FitResult)
    resid_df = _compute_residuals(bestfit)
    scenarios = ["NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg"]
    quantities = ["LogVirusLoad", "IL6", "Weight"]
    color_map = Dict(scenarios .=> PALETTE)
    marker_map = Dict(scenarios .=> MARKERS)

    fig = Figure(resolution = (900, 350))

    for (idx, qty) in enumerate(quantities)
        ax = Axis(fig[1, idx];
            xlabel = idx == 2 ? "Time (days)" : "",
            ylabel = QUANTITY_LABELS[qty],
            xticks = collect(0:floor(Int, maximum(resid_df.Day))))
        hlines!(ax, [0.0]; linestyle = :dash, color = :gray)

        for scen in scenarios
            sub = filter(row -> row.Quantity == qty && row.Scenario == scen, resid_df)
            scatter!(ax, sub.Day, sub.Residual;
                color = color_map[scen], marker = marker_map[scen], markersize = 8)
        end
    end

    legend = Legend(fig[0, :],
        [MarkerElement(color = color_map[scen], marker = marker_map[scen], markersize = 10)
         for scen in scenarios],
        [SCENARIO_LABELS[scen] for scen in scenarios]; orientation = :horizontal)
    legend.padding = 5

    return fig
end

end # module
