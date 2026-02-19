"""
Path and serialization utilities
================================

This file centralizes two practical concerns:

1. Resolving repository-relative paths from any entry script.
2. Saving/loading workflow outputs in a Julia-native format.

The R workflow uses `.Rds`; in Julia we use `.jld2` to preserve structured
objects (vectors, dictionaries, DataFrames, nested containers) without manual
JSON flattening.
"""

"""
Return the repository root based on the location of this source file.

The path is computed as:
`code/julia/src/paths-function.jl` -> go up three levels.
"""
function project_root()::String
    return normpath(joinpath(@__DIR__, "..", "..", ".."))
end

"""
Resolve a repository-relative path.

Example
-------
`repo_path("data", "processed-data", "processeddata.csv")`
"""
function repo_path(parts::AbstractString...)
    return joinpath(project_root(), parts...)
end

"""
Save any Julia object to disk in `.jld2` format.

Notes for readers new to Julia
------------------------------
`JLD2.jldsave` writes variables by name. We store the payload under the key
`"payload"` so the read side stays simple and stable.
"""
function save_julia_object(path::AbstractString, obj)
    mkpath(dirname(path))
    JLD2.jldsave(path; payload = obj)
    return path
end

"""
Load an object previously saved with `save_julia_object`.
"""
function load_julia_object(path::AbstractString)
    return JLD2.load(path, "payload")
end
