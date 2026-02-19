"""
Virus transform helpers
=======================

The R workflow transforms viral load before fitting. Keeping the transform in
one location avoids subtle mismatches across fitting, diagnostics, and plotting.

Current default
---------------
`log10_p1`: `log10(max(x, 0) + 1)`
"""

const virus_quantity_name = "VirusLoad"
const virus_transform_mode = "log10_p1"

"""
Apply the configured virus transform to a scalar.
"""
function transform_virus(x::Real; mode::String = virus_transform_mode)
    xval = Float64(x)
    if mode == "identity"
        return xval
    elseif mode == "log10_p1"
        return log10(max(xval, 0.0) + 1.0)
    elseif mode == "log10_clamp1"
        return log10(max(xval, 1.0))
    else
        error("Unknown virus transform mode: $mode")
    end
end

"""
Vector overload for convenience and readability.
"""
function transform_virus(x::AbstractVector{<:Real}; mode::String = virus_transform_mode)
    return [transform_virus(xi; mode = mode) for xi in x]
end

"""
Inverse virus transform.

Useful for plotting trajectories on the natural scale if needed.
"""
function inverse_transform_virus(x::Real; mode::String = virus_transform_mode)
    xval = Float64(x)
    if mode == "identity"
        return xval
    elseif mode == "log10_p1"
        return max(0.0, 10.0^xval - 1.0)
    elseif mode == "log10_clamp1"
        return 10.0^xval
    else
        error("Unknown virus transform mode: $mode")
    end
end
