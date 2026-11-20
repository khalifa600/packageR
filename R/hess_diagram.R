#' Colored energy diagram
#'
#' Draws an energy level diagram (ggplot2): reactants in blue, intermediates
#' in purple, products in green (exothermic), orange (endothermic) or gray
#' (thermoneutral). If `steps` contains at least two steps whose sum equals
#' `delta_h`, the diagram shows each step of the Hess combination.
#'
#' @param delta_h Reaction enthalpy (kJ/mol).
#' @param reactants_energy Energy level of the reactants (default 0).
#' @param steps Numeric (named) vector of step enthalpies, optional.
#' @param reactants_label,products_label Labels of the extreme levels.
#' @param title Plot title.
#' @return A `ggplot` object (display it with `print()`).
#' @examples
#' p <- hess_diagram(-890.3)
#' \donttest{
#' print(p)
#' }
#' @export
hess_diagram <- function(delta_h, reactants_energy = 0, steps = NULL,
                         reactants_label = "R\u00e9actifs",
                         products_label = "Produits",
                         title = "Diagramme \u00e9nerg\u00e9tique \u2013 Loi de Hess") {

  if (!is.numeric(delta_h) ||
      length(delta_h) != 1 ||
      !is.finite(delta_h)) {
    stop("delta_h doit etre un nombre fini.")
  }

  if (!is.numeric(reactants_energy) ||
      length(reactants_energy) != 1 ||
      !is.finite(reactants_energy)) {
    stop("reactants_energy doit etre un nombre fini.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Le package ggplot2 est requis : install.packages(\"ggplot2\")")
  }

  if (!is.null(steps)) {
    if (!is.numeric(steps) || any(!is.finite(steps))) {
      stop("steps doit etre un vecteur numerique fini.")
    }
    if (length(steps) < 2) {
      steps <- NULL
    } else if (abs(sum(steps) - delta_h) > 1e-6 * max(1, abs(delta_h))) {
      warning("La somme des etapes ne correspond pas a delta_h : diagramme simplifie.")
      steps <- NULL
    }
  }

  if (is.null(steps)) {
    energies <- c(reactants_energy, reactants_energy + delta_h)
    labels <- c(reactants_label, products_label)
  } else {
    energies <- reactants_energy + c(0, cumsum(as.numeric(steps)))
    sn <- names(steps)
    if (is.null(sn)) sn <- paste("\u00c9tape", seq_along(steps))
    k <- length(steps)
    labels <- c(reactants_label, paste("Apr\u00e8s", sn[-k]), products_label)
  }

  n <- length(energies)

  # Palette ---------------------------------------------------------------
  col_react <- "#3B82F6"
  col_mid <- "#8B5CF6"
  col_exo <- "#10B981"
  col_endo <- "#F97316"
  col_neutral <- "#64748B"

  class_txt <- hess_classify(delta_h)
  col_final <- switch(class_txt,
                      Exothermique = col_exo,
                      Endothermique = col_endo,
                      col_neutral)
  nature <- switch(class_txt,
                   Exothermique = "r\u00e9action exothermique \u2013 lib\u00e8re de l'\u00e9nergie",
                   Endothermique = "r\u00e9action endothermique \u2013 absorbe de l'\u00e9nergie",
                   "r\u00e9action thermoneutre")

  # Donnees ---------------------------------------------------------------
  half <- 0.22

  levels_df <- data.frame(
    x = seq_len(n),
    y = energies,
    label = labels,
    txt = sprintf("%.1f", energies),
    col = c(col_react, rep(col_mid, n - 2), col_final),
    stringsAsFactors = FALSE
  )

  dh_steps <- diff(energies)

  seg_df <- data.frame(
    x = levels_df$x[-n] + half,
    xend = levels_df$x[-1] - half,
    y = energies[-n],
    yend = energies[-1],
    dh = dh_steps,
    col = ifelse(dh_steps < 0, col_exo, ifelse(dh_steps > 0, col_endo, col_neutral)),
    vj = ifelse(dh_steps < 0, -0.4, 1.4),
    stringsAsFactors = FALSE
  )
  seg_df$mx <- (seg_df$x + seg_df$xend) / 2
  seg_df$my <- (seg_df$y + seg_df$yend) / 2
  seg_df$txt <- sprintf("%+.1f kJ/mol", seg_df$dh)

  rng <- range(energies)
  pad <- 0.28 * max(diff(rng), 1)

  # Graphique -------------------------------------------------------------
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = reactants_energy,
                        linetype = "dotted", colour = "#94A3B8") +
    ggplot2::geom_segment(
      data = seg_df,
      mapping = ggplot2::aes(x = .data$x, xend = .data$xend,
                             y = .data$y, yend = .data$yend,
                             colour = .data$col),
      linewidth = 1.3,
      arrow = grid::arrow(length = grid::unit(0.16, "inches"), type = "closed")
    ) +
    ggplot2::geom_segment(
      data = levels_df,
      mapping = ggplot2::aes(x = .data$x - half, xend = .data$x + half,
                             y = .data$y, yend = .data$y,
                             colour = .data$col),
      linewidth = 3.2, lineend = "round"
    ) +
    ggplot2::geom_text(
      data = levels_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                             label = .data$txt, colour = .data$col),
      vjust = -1.1, fontface = "bold", size = 4.2
    ) +
    ggplot2::geom_text(
      data = seg_df,
      mapping = ggplot2::aes(x = .data$mx, y = .data$my,
                             label = .data$txt, colour = .data$col,
                             vjust = .data$vj),
      hjust = -0.1, fontface = "bold", size = 3.9
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      breaks = levels_df$x, labels = levels_df$label,
      guide = ggplot2::guide_axis(n.dodge = if (n > 3) 2 else 1)
    ) +
    ggplot2::coord_cartesian(xlim = c(0.4, n + 0.6),
                             ylim = c(rng[1] - pad, rng[2] + pad)) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("\u0394H = %+.1f kJ/mol  \u2022  %s", delta_h, nature),
      x = NULL,
      y = "\u00c9nergie (kJ/mol)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17,
                                         colour = "#0F172A"),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 12,
                                            colour = col_final,
                                            margin = ggplot2::margin(b = 12)),
      axis.title.y = ggplot2::element_text(colour = "#475569", face = "bold"),
      axis.text.x = ggplot2::element_text(face = "bold", colour = "#1E293B",
                                          size = 12),
      axis.text.y = ggplot2::element_text(colour = "#475569"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#E2E8F0",
                                                 linewidth = 0.4),
      plot.background = ggplot2::element_rect(fill = "#F8FAFC", colour = NA),
      plot.margin = ggplot2::margin(16, 22, 12, 16)
    )
}
