# Complete HessR example --------------------------------------------------
library(HessR)

# 1. Combustion of methane (one reaction)
ex1 <- hess_example("methane")
hess_run(ex1$reactions)

# 2. Formation of methane: the coefficients are found automatically
ex2 <- hess_example("formation_methane")
res <- hess_run(ex2$reactions, target = ex2$target)
res$delta_h
res$coefficients

# 3. Formation of carbon monoxide
ex3 <- hess_example("monoxyde_carbone")
hess_run(ex3$reactions, target = ex3$target)

# 4. Save the plot
# ggplot2::ggsave("diagramme_hess.png", res$diagram, width = 9, height = 6, dpi = 300)
