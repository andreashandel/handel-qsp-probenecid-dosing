module TimeSeriesPlots

using CairoMakie
using Colors
using DataFrames

export plot_timeseries, save_timeseries_plot

const COLORS = [colorant"#0072B2", colorant"#009E73", colorant"#D55E00"]
const LINESTYLES = [:solid, :dash, :dot]
const MARKERS = [:circle, :utriangle, :rect]

function _dose_mapping(doses::Vector{Float64}, labels::Vector{String})
    if length(labels) != length(doses)
        error("Dose labels must match number of dose levels")
    end
    return Dict(doses .=> labels)
end

function _prepare_measurements(data::Union{Nothing, DataFrame}, dose_map::Dict{Float64, String})
    if data === nothing
        return Dict{String, DataFrame}()
    end
    df = copy(data)
    df.DoseLabel = [dose_map[Float64(d)] for d in df.Dose]
    virus = filter(:Quantity => ==("LogVirusLoad"), df)
    virus.Value .= 10 .^ virus.Value
    inn = filter(:Quantity => ==("IL6"), df)
    sym = filter(:Quantity => ==("Weight"), df)
    return Dict(
        "Virus" => virus,
        "Innate" => inn,
        "Symptom" => sym,
    )
end

"""
    plot_timeseries(; data, modelfit, tmax, dose_levels) -> Figure

Recreate the multi-panel time-series visualisation from
`code/plotting-code/timeseries-plot-function.R`.  `modelfit` must contain the
model trajectories (one row per time point and dose level) while `data`
optionally overlays experimental observations.
"""
function plot_timeseries(; data::Union{Nothing, DataFrame} = nothing,
    modelfit::DataFrame, tmax::Float64 = 7.0,
    dose_levels::Vector{String})

    unique_doses = sort(unique(Float64.(modelfit.Dose)))
    dose_map = _dose_mapping(unique_doses, dose_levels)
    color_map = Dict(dose_levels[i] => COLORS[i] for i in eachindex(dose_levels))
    marker_map = Dict(dose_levels[i] => MARKERS[i] for i in eachindex(dose_levels))

    df = copy(modelfit)
    df.DoseLabel = [dose_map[Float64(d)] for d in df.Dose]

    measurements = _prepare_measurements(data, dose_map)

    panels = [
        (:U, "Uninfected Cells", true, (1e2, 1e7)),
        (:I, "Infected Cells", true, (1e-2, 1e10)),
        (:V, "Virus Load", true, (1e-2, 1e10)),
        (:F, "Innate Response (IL-6)", false, (0, 2)),
        (:A, "Adaptive Response", true, (1e-2, 1e10)),
        (:S, "Morbidity (weight loss)", false, (0, 30)),
        (:Ad, "Drug depot compartment", true, (1e-5, 1e1)),
        (:Ac, "Drug central compartment", true, (1e-5, 1e0)),
        (:At, "Drug target compartment", true, (1e-6, 1e-1)),
    ]

    fig = Figure(resolution = (1200, 900))
    legend_handles = Tuple{Any, String}[]

    for (idx, (var, ylabel, logy, ylims)) in enumerate(panels)
        row = Int(ceil(idx / 3))
        col = idx - (row - 1) * 3
        axis_kwargs = (; ylabel = ylabel, xlabel = row == 3 ? "Days" : "")
        if logy
            axis = Axis(fig[row, col]; axis_kwargs..., yscale = log10)
        else
            axis = Axis(fig[row, col]; axis_kwargs...)
        end
        axis.ylimits = ylims
        axis.xlimits = (0, tmax)
        axis.xticks = collect(0:floor(Int, tmax))

        for (i, dose) in enumerate(unique_doses)
            label = dose_map[dose]
            colour = COLORS[i]
            ls = LINESTYLES[i]
            sub = filter(:Dose => ==(dose), df)
            lines!(axis, sub.time, sub[!, var]; color = colour, linestyle = ls, linewidth = 2)
            if idx == 1
                push!(legend_handles, (LineElement(color = colour, linestyle = ls, linewidth = 2), label))
            end
        end

        if var == :V && haskey(measurements, "Virus")
            meas = measurements["Virus"]
            scatter!(axis, meas.xvals, meas.Value;
                color = map(label -> color_map[label], meas.DoseLabel),
                marker = map(label -> marker_map[label], meas.DoseLabel),
                markersize = 8, strokewidth = 0.5, transparency = true)
        elseif var == :F && haskey(measurements, "Innate")
            meas = measurements["Innate"]
            scatter!(axis, meas.xvals, meas.Value;
                color = map(label -> color_map[label], meas.DoseLabel),
                marker = map(label -> marker_map[label], meas.DoseLabel),
                markersize = 8, strokewidth = 0.5, transparency = true)
        elseif var == :S && haskey(measurements, "Symptom")
            meas = measurements["Symptom"]
            scatter!(axis, meas.xvals, meas.Value;
                color = map(label -> color_map[label], meas.DoseLabel),
                marker = map(label -> marker_map[label], meas.DoseLabel),
                markersize = 8, strokewidth = 0.5, transparency = true)
        end
    end

    legend = Legend(fig[0, :], legend_handles; orientation = :horizontal)
    legend.padding = 5

    return fig
end

"""
    save_timeseries_plot(output_path; kwargs...)

Wrapper that calls [`plot_timeseries`] with the supplied keyword arguments and
persists the figure to `output_path`.
"""
function save_timeseries_plot(output_path::AbstractString; kwargs...)
    fig = plot_timeseries(; kwargs...)
    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

end # module
