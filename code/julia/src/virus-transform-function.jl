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
    if mode == "identity"
        return x
    elseif mode == "log10_p1"
        # Keep numeric type generic (e.g., ForwardDiff.Dual in stiff ODE solvers).
        return log10(max(x, zero(x)) + one(x))
    elseif mode == "log10_clamp1"
        return log10(max(x, one(x)))
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
    if mode == "identity"
        return x
    elseif mode == "log10_p1"
        return max(zero(x), (10.0^x) - one(x))
    elseif mode == "log10_clamp1"
        return 10.0^x
    else
        error("Unknown virus transform mode: $mode")
    end
end
