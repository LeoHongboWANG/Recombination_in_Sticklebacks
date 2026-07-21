```r
library(glmmTMB)
library(MuMIn)
library(dplyr)
library(ggplot2)
library(tibble)
library(car)

options(na.action = "na.fail")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Load input data
rrate.2.pop_avg <- read.csv("data/rrate_2_pop_avg.csv")
rrate.2.eco_avg <- read.csv("data/rrate_2_eco_avg.csv")

rrate.2.pop_avg.fw <- rrate.2.pop_avg %>%
  subset(Ecotype == "Freshwater")

rrate.2.pop_avg.ma <- rrate.2.pop_avg %>%
  subset(Ecotype == "Marine")


# -----------------------------
# Helper functions
# -----------------------------

get_r2 <- function(model) {
  r2 <- suppressWarnings(r.squaredGLMM(model))
  data.frame(
    R2m = as.numeric(r2[1, "R2m"]),
    R2c = as.numeric(r2[1, "R2c"])
  )
}

plot_residuals <- function(model, outfile, title) {
  png(outfile, width = 1200, height = 900, res = 150)
  plot(
    fitted(model),
    residuals(model),
    xlab = "Fitted values",
    ylab = "Residuals",
    main = title
  )
  abline(h = 0, col = "red", lwd = 2)
  dev.off()
}

make_formula <- function(response, terms, random = NULL) {
  if (length(terms) == 0) {
    rhs <- "1"
  } else {
    rhs <- paste(terms, collapse = " + ")
  }

  if (!is.null(random)) {
    rhs <- paste(rhs, random, sep = " + ")
  }

  as.formula(paste(response, "~", rhs))
}

partition_r2 <- function(data, response, main_effects, interactions = character(0), random = NULL) {
  main_formula <- make_formula(response, main_effects, random)
  main_model <- glmmTMB(
    main_formula,
    family = gaussian,
    data = data,
    na.action = na.fail
  )

  r2_main <- get_r2(main_model)$R2m

  full_terms <- c(main_effects, interactions)
  full_formula <- make_formula(response, full_terms, random)
  full_model <- glmmTMB(
    full_formula,
    family = gaussian,
    data = data,
    na.action = na.fail
  )

  r2_full <- get_r2(full_model)$R2m

  results <- data.frame(
    Term = character(),
    R2_increase = numeric(),
    stringsAsFactors = FALSE
  )

  for (term in main_effects) {
    reduced_terms <- setdiff(main_effects, term)

    reduced_formula <- make_formula(response, reduced_terms, random)

    reduced_model <- glmmTMB(
      reduced_formula,
      family = gaussian,
      data = data,
      na.action = na.fail
    )

    r2_reduced <- get_r2(reduced_model)$R2m
    r2_increase <- max(0, r2_main - r2_reduced)

    results <- rbind(
      results,
      data.frame(
        Term = term,
        R2_increase = r2_increase
      )
    )
  }

  for (term in interactions) {
    interaction_terms <- c(main_effects, term)

    interaction_formula <- make_formula(response, interaction_terms, random)

    interaction_model <- glmmTMB(
      interaction_formula,
      family = gaussian,
      data = data,
      na.action = na.fail
    )

    r2_interaction <- get_r2(interaction_model)$R2m
    r2_increase <- max(0, r2_interaction - r2_main)

    results <- rbind(
      results,
      data.frame(
        Term = term,
        R2_increase = r2_increase
      )
    )
  }

  results <- results %>%
    arrange(desc(R2_increase)) %>%
    mutate(
      R2_percent = R2_increase * 100,
      Relative_contribution = R2_increase / sum(R2_increase) * 100
    )

  attr(results, "R2m_main") <- r2_main
  attr(results, "R2m_full") <- r2_full

  return(results)
}


# -----------------------------
# Final models
# -----------------------------

best_model_fw <- glmmTMB(
  pi_scale ~ cpg_scale + ravg_scale + (1 | Population),
  family = gaussian,
  data = rrate.2.pop_avg.fw,
  na.action = na.fail
)

best_model_ma <- glmmTMB(
  pi_scale ~ cpg_scale + gene_density_scale + ravg_scale +
    cpg_scale:gene_density_scale +
    cpg_scale:ravg_scale +
    gene_density_scale:ravg_scale +
    (1 | Population),
  family = gaussian,
  data = rrate.2.pop_avg.ma,
  na.action = na.fail
)

best_model_eco <- glmmTMB(
  scale(pi) ~ cpg_scale + Ecotype + gene_density_scale + ravg_scale +
    cpg_scale:Ecotype +
    cpg_scale:ravg_scale +
    Ecotype:ravg_scale,
  family = gaussian,
  data = rrate.2.eco_avg,
  na.action = na.fail
)


# -----------------------------
# Model summaries
# -----------------------------

sink("results/model_summaries.txt")

cat("Freshwater model\n")
cat("================\n")
print(summary(best_model_fw))
cat("\nVariance components\n")
print(VarCorr(best_model_fw))
cat("\nR2\n")
print(get_r2(best_model_fw))

cat("\n\nMarine model\n")
cat("============\n")
print(summary(best_model_ma))
cat("\nVariance components\n")
print(VarCorr(best_model_ma))
cat("\nR2\n")
print(get_r2(best_model_ma))

cat("\n\nEcotype-level model\n")
cat("===================\n")
print(summary(best_model_eco))
cat("\nR2\n")
print(get_r2(best_model_eco))

sink()


# -----------------------------
# Residual plots
# -----------------------------

plot_residuals(
  best_model_fw,
  "figures/residuals_fw.png",
  "Freshwater residuals"
)

plot_residuals(
  best_model_ma,
  "figures/residuals_ma.png",
  "Marine residuals"
)

plot_residuals(
  best_model_eco,
  "figures/residuals_eco.png",
  "Ecotype-level residuals"
)


# -----------------------------
# R2 partitioning
# -----------------------------

fw_main_effects <- c(
  "ravg_scale",
  "cpg_scale"
)

fw_interactions <- character(0)

ma_main_effects <- c(
  "ravg_scale",
  "cpg_scale",
  "gene_density_scale"
)

ma_interactions <- c(
  "cpg_scale:gene_density_scale",
  "cpg_scale:ravg_scale",
  "gene_density_scale:ravg_scale"
)

eco_main_effects <- c(
  "ravg_scale",
  "cpg_scale",
  "gene_density_scale",
  "Ecotype"
)

eco_interactions <- c(
  "cpg_scale:Ecotype",
  "cpg_scale:ravg_scale",
  "Ecotype:ravg_scale"
)

fw_r2 <- partition_r2(
  data = rrate.2.pop_avg.fw,
  response = "pi_scale",
  main_effects = fw_main_effects,
  interactions = fw_interactions,
  random = "(1 | Population)"
)

ma_r2 <- partition_r2(
  data = rrate.2.pop_avg.ma,
  response = "pi_scale",
  main_effects = ma_main_effects,
  interactions = ma_interactions,
  random = "(1 | Population)"
)

eco_r2 <- partition_r2(
  data = rrate.2.eco_avg,
  response = "scale(pi)",
  main_effects = eco_main_effects,
  interactions = eco_interactions,
  random = NULL
)

write.csv(fw_r2, "results/r2_partition_fw.csv", row.names = FALSE)
write.csv(ma_r2, "results/r2_partition_ma.csv", row.names = FALSE)
write.csv(eco_r2, "results/r2_partition_eco.csv", row.names = FALSE)


# -----------------------------
# VIF checks
# -----------------------------

fw_lm <- lm(
  pi_scale ~ cpg_scale + ravg_scale,
  data = rrate.2.pop_avg.fw,
  na.action = na.fail
)

ma_lm <- lm(
  pi_scale ~ cpg_scale + gene_density_scale + ravg_scale,
  data = rrate.2.pop_avg.ma,
  na.action = na.fail
)

eco_lm <- lm(
  scale(pi) ~ cpg_scale + Ecotype + gene_density_scale + ravg_scale,
  data = rrate.2.eco_avg,
  na.action = na.fail
)

sink("results/vif_values.txt")

cat("Freshwater VIF\n")
cat("================\n")
print(vif(fw_lm))

cat("\nMarine VIF\n")
cat("==========\n")
print(vif(ma_lm))

cat("\nEcotype-level VIF\n")
cat("=================\n")
print(vif(eco_lm))

sink()


# -----------------------------
# Variance explained ring plot
# -----------------------------

make_variance_data <- function(r2_table, model, group_name, level) {
  r2 <- get_r2(model)

  random_r2 <- max(0, r2$R2c - r2$R2m)

  out <- bind_rows(
    tibble(
      Term = "Population",
      R2_increase = random_r2
    ),
    as_tibble(r2_table[, c("Term", "R2_increase")])
  ) %>%
    mutate(
      Group = group_name,
      Level = level,
      Variance = R2_increase * 100
    )

  return(out)
}

variance_data <- bind_rows(
  make_variance_data(
    r2_table = ma_r2,
    model = best_model_ma,
    group_name = "Marine",
    level = 3
  ),
  make_variance_data(
    r2_table = fw_r2,
    model = best_model_fw,
    group_name = "Freshwater",
    level = 2
  )
) %>%
  group_by(Group) %>%
  mutate(
    Variance_norm = Variance / sum(Variance) * 100,
    Label = paste0(Term, "\n", round(Variance_norm, 1), "%")
  ) %>%
  ungroup()

p_ring <- ggplot(
  variance_data,
  aes(x = Level, y = Variance_norm, fill = Group)
) +
  geom_bar(
    stat = "identity",
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 2.6,
    fontface = "bold"
  ) +
  coord_polar("y", start = 0) +
  scale_x_continuous(
    breaks = c(2, 3),
    limits = c(0, 4)
  ) +
  scale_fill_manual(
    values = c(
      "Marine" = "#3b83a1",
      "Freshwater" = "#f2b764"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Variance explained",
    fill = "Group"
  )

ggsave(
  "figures/variance_explained_ring.png",
  p_ring,
  width = 7,
  height = 7,
  dpi = 300
)


# -----------------------------
# Save session information
# -----------------------------

sink("results/sessionInfo.txt")
sessionInfo()
sink()
