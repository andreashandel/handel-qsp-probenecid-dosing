# -----------------------------------------------------------------------------
# fixed-parameter-sampling-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Create a list of fixed-parameter samples using Latin Hypercube Sampling
#   (LHS). These samples are used to explore uncertainty by refitting the model
#   across multiple fixed-parameter scenarios.
# -----------------------------------------------------------------------------

library(lhs) # Latin Hypercube Sampling.

#' Create a list of fixed-parameter samples.
#'
#' @param fixedpars Named numeric vector of baseline fixed parameters.
#' @param nsamp Number of random samples to generate in addition to the baseline.
#' @param seed RNG seed for reproducibility.
#' @param lower_factor Multiply baseline by this value for lower bounds.
#' @param upper_factor Multiply baseline by this value for upper bounds.
#' @return A list of named numeric vectors. The first element is always the
#'   baseline fixed parameter set, followed by nsamp samples (if nsamp > 0).
#'
make_fixed_parameter_samples <- function(
  fixedpars,
  nsamp,
  seed = 1234,
  lower_factor = 0.5,
  upper_factor = 2
) {
  # Always include the baseline fixed parameter set first.
  samples_list <- list(fixedpars)

  # If no sampling requested, return baseline only.
  if (nsamp <= 0) {
    return(samples_list)
  }

  # Seed for reproducibility of the LHS.
  set.seed(seed)

  # Define the sampling bounds around each fixed parameter.
  lower <- lower_factor * fixedpars
  upper <- upper_factor * fixedpars

  # Latin Hypercube in [0,1], then scale to [lower, upper].
  U <- lhs::randomLHS(nsamp, length(fixedpars))
  M <- sweep(U, 2, (upper - lower), "*")
  M <- sweep(M, 2, lower, "+")
  colnames(M) <- names(fixedpars)

  # Convert each row into a named vector matching the baseline parameters.
  samples_list_2 <- lapply(seq_len(nsamp), function(i) {
    setNames(as.numeric(M[i, ]), names(fixedpars))
  })

  # Return baseline + sampled parameter sets.
  c(samples_list, samples_list_2)
}
