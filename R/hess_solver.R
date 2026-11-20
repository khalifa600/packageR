#' Find the combination coefficients
#'
#' Solves `reactions %*% x = target` using the pseudo-inverse (SVD). Works for
#' square, overdetermined and underdetermined systems (minimum-norm solution).
#' The `exact` field indicates whether the target is reproduced.
#'
#' @param reactions Numeric matrix (one column per reaction).
#' @param target Target vector (one value per row of the matrix).
#' @param tolerance Tolerance on the residual (positive number).
#' @return List: `coefficients`, `residual`, `exact`, `rank`, `unique`.
#' @examples
#' hess_solver(matrix(c(1, 0, 0, 1), nrow = 2), c(2, 3))$coefficients
#' @export
hess_solver <- function(reactions, target, tolerance = 1e-8) {

  if (!is.matrix(reactions) || !is.numeric(reactions)) {
    stop("reactions doit etre une matrice numerique.")
  }

  if (!is.numeric(target) ||
      length(target) != nrow(reactions)) {
    stop("target doit avoir une valeur par ligne de la matrice.")
  }

  if (any(!is.finite(reactions)) ||
      any(!is.finite(target))) {
    stop("La matrice et la cible doivent contenir des valeurs finies.")
  }

  if (!is.numeric(tolerance) ||
      length(tolerance) != 1 ||
      !is.finite(tolerance) ||
      tolerance <= 0) {
    stop("tolerance doit etre un nombre positif.")
  }

  sv <- svd(reactions)
  rank_tol <- max(dim(reactions)) * .Machine$double.eps * max(sv$d, 0)
  keep <- sv$d > rank_tol
  rank <- sum(keep)

  if (rank == 0) {
    solution <- numeric(ncol(reactions))
  } else {
    u <- sv$u[, keep, drop = FALSE]
    v <- sv$v[, keep, drop = FALSE]
    solution <- as.vector(v %*% ((t(u) %*% as.numeric(target)) / sv$d[keep]))
  }

  solution[abs(solution) < 1e-10] <- 0
  solution <- round(solution, 10)

  residual <- as.vector(reactions %*% solution - as.numeric(target))
  exact <- sqrt(sum(residual^2)) <= tolerance

  if (!exact) {
    warning("La cible n'est pas reproduite exactement.")
  }

  names(solution) <- colnames(reactions)

  list(
    coefficients = solution,
    residual = residual,
    exact = exact,
    rank = rank,
    unique = (rank == ncol(reactions))
  )
}
