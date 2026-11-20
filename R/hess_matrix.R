#' Build the stoichiometric matrix
#'
#' Automatically builds, from the reactions (and an optional target reaction),
#' the matrix whose columns are the reactions and whose rows are the species
#' (products positive, reactants negative).
#'
#' @param reactions List of reactions.
#' @param target Target reaction (list with `reactants` and `products`),
#'   optional.
#' @return A list: `matrix`, `target` (vector or `NULL`) and `species`.
#' @examples
#' ex <- hess_example("formation_methane")
#' hess_matrix(ex$reactions, ex$target)$matrix
#' @export
hess_matrix <- function(reactions, target = NULL) {

  if (.hess_is_reaction(reactions)) {
    reactions <- list(reactions)
  }

  if (!is.list(reactions) || length(reactions) == 0) {
    stop("Fournissez une liste de reactions.")
  }

  for (i in seq_along(reactions)) {
    .hess_validate_reaction(reactions[[i]], paste("Reaction", i), need_h = FALSE)
  }

  all_reactions <- reactions
  if (!is.null(target)) {
    .hess_validate_reaction(target, "Cible", need_h = FALSE)
    all_reactions[[length(all_reactions) + 1L]] <- target
  }

  species <- unique(unlist(lapply(all_reactions, function(r) {
    c(names(r$reactants), names(r$products))
  })))

  net_vector <- function(r) {
    v <- numeric(length(species))
    names(v) <- species
    v[names(r$reactants)] <- v[names(r$reactants)] - r$reactants
    v[names(r$products)] <- v[names(r$products)] + r$products
    v
  }

  mat <- matrix(
    unlist(lapply(reactions, net_vector)),
    nrow = length(species),
    dimnames = list(species, NULL)
  )

  cn <- names(reactions)
  if (is.null(cn)) cn <- rep("", length(reactions))
  empty <- is.na(cn) | cn == ""
  cn[empty] <- paste0("R", which(empty))
  colnames(mat) <- cn

  list(
    matrix = mat,
    target = if (!is.null(target)) net_vector(target) else NULL,
    species = species
  )
}
