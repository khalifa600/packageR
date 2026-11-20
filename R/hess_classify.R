#' Classify a reaction according to the sign of its enthalpy
#'
#' @param delta_h Reaction enthalpy (finite number).
#' @param tolerance "Thermoneutrality" zone around zero (>= 0).
#' @return `"Exothermique"` (exothermic), `"Endothermique"` (endothermic) or
#'   `"Thermoneutre"` (thermoneutral).
#' @examples
#' hess_classify(-890.3)
#' @export
hess_classify <- function(delta_h, tolerance = 1e-8) {

  if (!is.numeric(delta_h) ||
      length(delta_h) != 1 ||
      !is.finite(delta_h)) {
    stop("delta_h doit etre un nombre fini.")
  }

  if (!is.numeric(tolerance) ||
      length(tolerance) != 1 ||
      !is.finite(tolerance) ||
      tolerance < 0) {
    stop("La tolerance doit etre positive ou nulle.")
  }

  if (delta_h < -tolerance) {
    return("Exothermique")
  }

  if (delta_h > tolerance) {
    return("Endothermique")
  }

  "Thermoneutre"
}
