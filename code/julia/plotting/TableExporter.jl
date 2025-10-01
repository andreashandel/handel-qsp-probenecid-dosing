module TableExporter

using DataFrames
using PrettyTables
using Printf

include("../analysis/FitModel.jl")
using .ModelFitting: FitResult

export export_parameter_tables

"""
    export_parameter_tables(bestfit_list; output_dir)

Generate markdown tables equivalent to the `gt` outputs from the R workflow.
Each element of `bestfit_list` produces a file `parametertabX.md` in
`output_dir`.
"""
function export_parameter_tables(bestfit_list::Vector{FitResult};
    output_dir::AbstractString = "results/tables")
    mkpath(output_dir)
    for (idx, bestfit) in enumerate(bestfit_list)
        df = DataFrame(
            Parameter = string.(bestfit.fitparnames),
            Value = bestfit.solution,
            Label = bestfit.parlabels,
        )
        outfile = joinpath(output_dir, "parametertab$(idx).md")
        open(outfile, "w") do io
            pretty_table(io, df; backend = Val(:markdown),
                formatters = ((v, i, j) -> j == 2 ? @sprintf("%.6g", v) : string(v)))
        end
    end
    return nothing
end

end # module
