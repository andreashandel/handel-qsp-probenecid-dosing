# -----------------------------------------------------------------------------
# virus-transform-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Centralize virus load transformations so data, model predictions, and plots
#   stay consistent. Adjust the transform in ONE place below.
# -----------------------------------------------------------------------------

# Name used for virus quantity throughout the codebase.
virus_quantity_name <- "VirusLoad"

# Transformation mode for virus values. Options:
#   - "identity": no transform (keep in natural units)
#   - "log10_p1": log10(V + 1)
#   - "log10_clamp1": log10(max(V, 1))
virus_transform_mode <- "log10_p1"

# Apply the selected virus transform.
transform_virus <- function(x, mode = virus_transform_mode) {
  x <- as.numeric(x)

  switch(
    mode,
    identity = x,
    log10_p1 = log10(pmax(x, 0) + 1),
    log10_clamp1 = log10(pmax(x, 1)),
    stop("Unknown virus_transform_mode: ", mode)
  )
}

# Invert the selected virus transform (for plotting raw values).
inverse_transform_virus <- function(x, mode = virus_transform_mode) {
  x <- as.numeric(x)

  switch(
    mode,
    identity = x,
    log10_p1 = pmax(0, 10^x - 1),
    log10_clamp1 = 10^x,
    stop("Unknown virus_transform_mode: ", mode)
  )
}
