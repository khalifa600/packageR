# HessR

R package to apply **Hess's law** to chemical reactions (combustion, formation, etc.).

**A single function is enough: `hess_run()`.** As soon as the engineer provides their reactions, all the functions are activated together: coefficient solving, reversal, multiplication, combination, ΔH calculation, atom check, classification, summary and colored diagram.

## Installation (RStudio)

1. Unzip the `HessR` folder, then double-click `HessR.Rproj`.
2. In the console:

```r
install.packages(c("devtools", "testthat", "ggplot2", "rlang"))
devtools::document()   # generates the help files (man/)
devtools::test()       # runs the tests
devtools::install()    # installs the package
library(HessR)
```

## Quick start

### 1. A single reaction (combustion of methane)

```r
library(HessR)

methane <- list(
  reactants = c("CH4(g)" = 1, "O2(g)" = 2),
  products  = c("CO2(g)" = 1, "H2O(l)" = 2),
  delta_h   = -890.3            # kJ/mol
)

hess_run(methane)
```

### 2. Several reactions + target reaction (coefficients found automatically)

```r
ex <- hess_example("formation_methane")   # C + 2 H2 -> CH4

res <- hess_run(ex$reactions, target = ex$target)
res$delta_h        # -74.8 kJ/mol
res$coefficients   # 1, 2, -1
```

### 3. Coefficients chosen by hand

```r
hess_run(ex$reactions, coefficients = c(1, 2, -1))
```

## What `hess_run()` does

| Step | Function |
|---|---|
| Stoichiometric matrix | `hess_matrix()` |
| Coefficients (if a target is given) | `hess_solver()` |
| Reversal (coefficient < 0) | `hess_reverse()` |
| Multiplication | `hess_multiply()` |
| Overall equation | `hess_combine()` |
| Total enthalpy | `hess_calculate()` |
| Atom check | `hess_check_balance()` |
| Classification | `hess_classify()` |
| Console summary | `hess_summary()` |
| Colored diagram | `hess_diagram()` |

Useful options: `show_diagram = FALSE`, `verbose = FALSE`, `check_balance = FALSE`, `reactants_energy = 0`.

The result contains `delta_h`, `classification`, `equation`, `coefficients`, `steps`, `atom_check`, `target_reached` and `diagram` (a `ggplot` object, which can be modified or saved with `ggplot2::ggsave()`).

## Reaction format

```r
list(
  reactants = c("A" = 1, "B" = 2),   # named vector, coefficients > 0
  products  = c("C" = 1),
  delta_h   = -100                    # finite number, kJ/mol
)
```

Writing the formulas with their physical state (`"H2O(l)"`, `"CO2(g)"`, `"C(s)"`) enables the automatic atom check. The names of the reaction list are used as labels in the summary and in the diagram.

## Units

Enthalpies must be expressed in a consistent unit (kJ/mol of reaction as written).

## Limitations

HessR computes enthalpy changes. It does not compute the real efficiency of an engine or of an installation (heat losses, operating conditions). The physical states must be consistent with the data used.
