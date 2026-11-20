# Fonctions internes (non exportees) ---------------------------------------

.hess_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(name, " doit etre TRUE ou FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

.hess_is_reaction <- function(r) {
  is.list(r) && !is.null(r$reactants) && !is.null(r$products) &&
    !is.null(r$delta_h)
}

.hess_validate_reaction <- function(r, label = "Reaction", need_h = TRUE) {

  if (!is.list(r) || is.null(r$reactants) || is.null(r$products) ||
      (need_h && is.null(r$delta_h))) {
    msg <- if (need_h) "reactants, products et delta_h" else "reactants et products"
    stop(label, " : format invalide (", msg, " requis).", call. = FALSE)
  }

  for (side in c("reactants", "products")) {
    v <- r[[side]]
    if (!is.numeric(v) || length(v) == 0 || is.null(names(v)) ||
        any(is.na(names(v))) || any(names(v) == "") ||
        anyDuplicated(names(v)) > 0 ||
        any(!is.finite(v)) || any(v <= 0)) {
      stop(label, " : ", side,
           " doit etre un vecteur numerique nomme ",
           "(noms uniques, coefficients positifs).", call. = FALSE)
    }
  }

  if (need_h || !is.null(r$delta_h)) {
    if (!is.numeric(r$delta_h) || length(r$delta_h) != 1 ||
        !is.finite(r$delta_h)) {
      stop(label, " : delta_h doit etre un nombre fini.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

.hess_side <- function(x) {
  if (length(x) == 0) return("(aucun)")
  coef <- vapply(x, function(v) {
    if (abs(v - 1) < 1e-12) "" else paste0(format(v, digits = 6, trim = TRUE), " ")
  }, character(1))
  paste0(coef, names(x), collapse = " + ")
}

.hess_equation <- function(reactants, products) {
  paste(.hess_side(reactants), "->", .hess_side(products))
}

# Analyse d'une formule chimique ("Ca(OH)2(s)") -> vecteur nomme d'atomes
.hess_parse_formula <- function(formula) {

  f <- sub("\\((s|l|g|aq)\\)$", "", formula)
  tokens <- regmatches(
    f, gregexpr("[A-Z][a-z]?|[0-9]+(\\.[0-9]+)?|\\(|\\)", f)
  )[[1]]

  if (!identical(paste(tokens, collapse = ""), f)) {
    stop("Formule chimique non reconnue : ", formula, call. = FALSE)
  }

  add <- function(vec, el, k) {
    if (el %in% names(vec)) vec[el] <- vec[el] + k else vec[el] <- k
    vec
  }

  stack <- list(numeric(0))
  n <- length(tokens)
  i <- 1L

  while (i <= n) {
    tk <- tokens[i]

    if (grepl("^[A-Z]", tk)) {
      k <- 1
      if (i < n && grepl("^[0-9]", tokens[i + 1L])) {
        k <- as.numeric(tokens[i + 1L])
        i <- i + 1L
      }
      top <- length(stack)
      stack[[top]] <- add(stack[[top]], tk, k)

    } else if (tk == "(") {
      stack[[length(stack) + 1L]] <- numeric(0)

    } else if (tk == ")") {
      if (length(stack) < 2) {
        stop("Parentheses non equilibrees dans : ", formula, call. = FALSE)
      }
      inner <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      k <- 1
      if (i < n && grepl("^[0-9]", tokens[i + 1L])) {
        k <- as.numeric(tokens[i + 1L])
        i <- i + 1L
      }
      top <- length(stack)
      for (el in names(inner)) {
        stack[[top]] <- add(stack[[top]], el, inner[[el]] * k)
      }

    } else {
      stop("Formule chimique non reconnue : ", formula, call. = FALSE)
    }

    i <- i + 1L
  }

  if (length(stack) != 1) {
    stop("Parentheses non equilibrees dans : ", formula, call. = FALSE)
  }

  stack[[1]]
}
