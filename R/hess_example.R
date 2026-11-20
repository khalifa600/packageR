#' Example datasets
#'
#' Provides ready-to-use reactions (usual values in kJ/mol).
#'
#' * `"methane"`: combustion of methane (-890.3 kJ/mol).
#' * `"formation_methane"`: formation of methane from combustions
#'   (expected result: -74.8 kJ/mol).
#' * `"monoxyde_carbone"`: formation of CO from two combustions
#'   (expected result: -110.5 kJ/mol).
#'
#' @param case Name of the example case.
#' @return List: `reactions`, `target` (or `NULL`) and `description`.
#' @examples
#' ex <- hess_example("formation_methane")
#' res <- hess_run(ex$reactions, target = ex$target,
#'                 show_diagram = FALSE, verbose = FALSE)
#' res$delta_h
#' @export
hess_example <- function(case = c("methane", "formation_methane",
                                  "monoxyde_carbone")) {

  case <- match.arg(case)

  combustion_methane <- list(
    reactants = c("CH4(g)" = 1, "O2(g)" = 2),
    products = c("CO2(g)" = 1, "H2O(l)" = 2),
    delta_h = -890.3
  )

  combustion_carbone <- list(
    reactants = c("C(s)" = 1, "O2(g)" = 1),
    products = c("CO2(g)" = 1),
    delta_h = -393.5
  )

  formation_eau <- list(
    reactants = c("H2(g)" = 1, "O2(g)" = 0.5),
    products = c("H2O(l)" = 1),
    delta_h = -285.8
  )

  combustion_co <- list(
    reactants = c("CO(g)" = 1, "O2(g)" = 0.5),
    products = c("CO2(g)" = 1),
    delta_h = -283.0
  )

  switch(case,
    methane = list(
      reactions = list("Combustion du methane" = combustion_methane),
      target = NULL,
      description = "Combustion complete du methane (eau liquide)."
    ),
    formation_methane = list(
      reactions = list(
        "Combustion du carbone" = combustion_carbone,
        "Formation de l'eau" = formation_eau,
        "Combustion du methane" = combustion_methane
      ),
      target = list(
        reactants = c("C(s)" = 1, "H2(g)" = 2),
        products = c("CH4(g)" = 1)
      ),
      description = "Formation du methane : C + 2 H2 -> CH4."
    ),
    monoxyde_carbone = list(
      reactions = list(
        "Combustion du carbone" = combustion_carbone,
        "Combustion du CO" = combustion_co
      ),
      target = list(
        reactants = c("C(s)" = 1, "O2(g)" = 0.5),
        products = c("CO(g)" = 1)
      ),
      description = "Formation du monoxyde de carbone : C + 1/2 O2 -> CO."
    )
  )
}
