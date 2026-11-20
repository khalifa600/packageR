#' Summary of an analysis
#'
#' Prints a readable summary: overall equation, enthalpy, classification,
#' coefficients, steps and atom check (when available).
#'
#' @param result Result of [hess_run()], [hess_combine()] or any list
#'   containing `delta_h` (and ideally `reactants` and `products`).
#' @return `result`, invisibly.
#' @examples
#' hess_summary(list(delta_h = -100, reactants = c(A = 1), products = c(B = 1)))
#' @export
hess_summary <- function(result) {

  if (!is.list(result) || is.null(result$delta_h)) {
    stop("Resultat HessR invalide.")
  }

  line <- strrep("=", 48)

  cat("\n", line, "\n", sep = "")
  cat("                 RESUME HESSR\n")
  cat(line, "\n", sep = "")

  if (!is.null(result$reactants) && !is.null(result$products)) {
    cat("Equation bilan  : ",
        .hess_equation(result$reactants, result$products), "\n", sep = "")
  }

  cat(sprintf("Enthalpie (dH)  : %.2f kJ/mol\n", result$delta_h))
  cat("Classification  : ", hess_classify(result$delta_h), "\n", sep = "")

  if (!is.null(result$reference_delta_h)) {
    cat(sprintf("Valeur cible    : %.2f kJ/mol (ecart %.2e)\n",
                result$reference_delta_h,
                result$delta_h - result$reference_delta_h))
  }

  if (!is.null(result$coefficients)) {
    cat("\nCoefficients de combinaison :\n")
    cf <- result$coefficients
    nm <- names(cf)
    if (is.null(nm)) nm <- paste("R", seq_along(cf))
    for (i in seq_along(cf)) {
      cat(sprintf("  %-28s x %s\n", nm[i], format(cf[[i]], digits = 6)))
    }
  }

  if (!is.null(result$steps) && length(result$steps) > 0) {
    cat("\nContribution de chaque etape (kJ/mol) :\n")
    st <- result$steps
    nm <- names(st)
    if (is.null(nm)) nm <- paste("Etape", seq_along(st))
    for (i in seq_along(st)) {
      cat(sprintf("  %-28s %+10.2f\n", nm[i], st[[i]]))
    }
  }

  if (!is.null(result$atom_check)) {
    cat("\nBilan atomique  : ", result$atom_check$message, "\n", sep = "")
  }

  cat(line, "\n", sep = "")

  invisible(result)
}
