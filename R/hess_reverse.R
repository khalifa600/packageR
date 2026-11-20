#' Reverse a reaction
#'
#' Swaps reactants and products and changes the sign of the enthalpy.
#'
#' @param reaction List with `reactants`, `products` (named vectors) and
#'   `delta_h` (number).
#' @return The reversed reaction (list).
#' @examples
#' r <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -50)
#' hess_reverse(r)
#' @export
hess_reverse <- function(reaction) {

  if (!is.list(reaction) ||
      is.null(reaction$reactants) ||
      is.null(reaction$products) ||
      is.null(reaction$delta_h)) {
    stop("La reaction doit contenir reactants, products et delta_h.")
  }

  if (!is.numeric(reaction$delta_h) ||
      length(reaction$delta_h) != 1 ||
      !is.finite(reaction$delta_h)) {
    stop("delta_h doit etre un nombre fini.")
  }

  list(
    reactants = reaction$products,
    products = reaction$reactants,
    delta_h = -reaction$delta_h
  )
}
