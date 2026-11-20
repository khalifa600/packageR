#' Total enthalpy from coefficients and enthalpies
#'
#' Computes the weighted sum `sum(coefficients * enthalpies)`.
#'
#' @param coefficients Numeric vector of combination coefficients.
#' @param enthalpies Numeric vector of enthalpies (kJ/mol), same length.
#' @return A number: the total enthalpy.
#' @examples
#' hess_calculate(c(1, 2), c(-100, -50))
#' @export
hess_calculate <- function(coefficients, enthalpies) {

  if (!is.numeric(coefficients) || !is.numeric(enthalpies)) {
    stop("Les coefficients et les enthalpies doivent etre numeriques.")
  }

  if (length(coefficients) != length(enthalpies)) {
    stop("Les deux vecteurs doivent avoir la meme longueur.")
  }

  if (length(coefficients) == 0 ||
      any(!is.finite(coefficients)) ||
      any(!is.finite(enthalpies))) {
    stop("Les vecteurs doivent contenir des valeurs finies.")
  }

  sum(coefficients * enthalpies)
}
