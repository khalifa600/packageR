#' Combine several reactions (Hess's law)
#'
#' Adds up the reactions weighted by `coefficients`, cancels the species
#' appearing on both sides and sums the enthalpies.
#'
#' @param reactions List of reactions (`reactants`, `products`, `delta_h`).
#' @param coefficients Combination coefficients (negative values allowed).
#'   Defaults to 1 for each reaction.
#' @return A `hess_reaction` object (list) with `reactants`, `products`,
#'   `delta_h` and `balance` (net balance per species).
#' @examples
#' r1 <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -50)
#' r2 <- list(reactants = c(B = 1), products = c(C = 1), delta_h = -20)
#' hess_combine(list(r1, r2))$delta_h
#' @export
hess_combine <- function(reactions, coefficients = NULL) {

  if (!is.list(reactions) || length(reactions) == 0) {
    stop("Fournissez une liste de reactions.")
  }

  if (is.null(coefficients)) {
    coefficients <- rep(1, length(reactions))
  }

  if (!is.numeric(coefficients) ||
      length(coefficients) != length(reactions) ||
      any(!is.finite(coefficients))) {
    stop("Les coefficients doivent etre numeriques et de bonne longueur.")
  }

  species_balance <- numeric(0)
  total_h <- 0

  for (i in seq_along(reactions)) {

    r <- reactions[[i]]
    .hess_validate_reaction(r, paste("Reaction", i))

    for (sp in names(r$reactants)) {
      if (!sp %in% names(species_balance)) {
        species_balance[sp] <- 0
      }
      species_balance[sp] <-
        species_balance[sp] - coefficients[i] * r$reactants[[sp]]
    }

    for (sp in names(r$products)) {
      if (!sp %in% names(species_balance)) {
        species_balance[sp] <- 0
      }
      species_balance[sp] <-
        species_balance[sp] + coefficients[i] * r$products[[sp]]
    }

    total_h <- total_h + coefficients[i] * r$delta_h
  }

  tol <- 1e-10
  species_balance[abs(species_balance) < tol] <- 0

  reactants <- -species_balance[species_balance < 0]
  products <- species_balance[species_balance > 0]

  structure(
    list(
      reactants = reactants,
      products = products,
      delta_h = as.numeric(total_h),
      balance = species_balance
    ),
    class = "hess_reaction"
  )
}
