# Plotting utilities (Julia)

These Julia scripts replicate the figure generation routines from `code/plotting-code`.  The plotting code is implemented with [CairoMakie.jl](https://makie.juliaplots.org/stable/) to provide high-quality vector and raster outputs similar to the original ggplot2 + patchwork workflow.

## Files

- `DoseResponsePlots.jl` – reproduces `dose-response-plot-function.R` and `make-dosing-figures.R`.
- `TimeSeriesPlots.jl` – covers `timeseries-plot-function.R`, `make-bestfit-figure.R`, and `make-timeseries-figures.R`.
- `DiagnosticPlots.jl` – implements the single- and multi-panel residual visualisations from `make-diagnostic-figures.R` and `make-combined-diagnostic-figures.R`.
- `TableExporter.jl` – the Julia analogue of `make-tables.R`, producing formatted tables via PrettyTables.

Refer to `code/julia/scripts` for ready-to-run entry points that call into these modules and save the corresponding figures/tables under `results/`.
