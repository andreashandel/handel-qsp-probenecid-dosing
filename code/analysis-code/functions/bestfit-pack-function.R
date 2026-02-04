# -----------------------------------------------------------------------------
# bestfit-pack-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Attach consistent metadata to nloptr fit results so all downstream scripts
#   (plots, tables, diagnostics) can reuse a standard bestfit structure.
#
# WHY THIS FILE EXISTS
#   The original fitting scripts packaged results in a particular shape. This
#   helper replicates that shape so newer scripts remain backward compatible.
# -----------------------------------------------------------------------------

#' Attach standardized metadata to an nloptr fit object.
#'
#' @param fit_result Object returned by nloptr::nloptr.
#' @param fitparnames Character vector of fitted parameter names.
#' @param fixedpars Named vector of fixed parameters (including fixed sigmas).
#' @param Y0 Named vector of initial conditions.
#' @param fitdata The processed fitting data (data.frame).
#' @param parlabels Named vector of parameter labels for tables.
#' @param algorithm Name of the optimization algorithm.
#' @param logfit Logical flag: TRUE if parameters were fit in log-space.
#' @return The fit_result augmented with standard fields.
#'
pack_bestfit <- function(
  fit_result,
  fitparnames,
  fixedpars,
  Y0,
  fitdata,
  parlabels,
  algorithm,
  logfit = TRUE
) {
  # Start from the solution vector produced by nloptr.
  params <- fit_result$solution
  names(params) <- fitparnames

  # Transform back to natural scale if the optimizer ran in log-space.
  if (isTRUE(logfit)) {
    params <- exp(params)
  }

  # Attach the standard fields expected by the plotting and table scripts.
  fit_result$parstring <- paste0("c(", paste(as.numeric(params), collapse = ", "), ")")
  fit_result$fitpars <- params
  fit_result$fitparnames <- fitparnames
  fit_result$fixedpars <- fixedpars
  fit_result$Y0 <- Y0
  fit_result$fitdata <- fitdata
  fit_result$parlabels <- parlabels
  fit_result$algorithm <- algorithm

  # Return the augmented object.
  fit_result
}
