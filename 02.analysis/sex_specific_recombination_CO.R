library(glmmTMB)
library(MuMIn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)

options(na.action = "na.fail")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# -----------------------------
# Input files
# -----------------------------

recom_sexdifferent <- read.table(
  "data/recom_sexdifferent.txt",
  header = FALSE,
  sep = "\t"
)

colnames(recom_sexdifferent)[1:7] <- c(
  "Population",
  "Chr",
  "Chr_length",
  "male_cM",
  "male_r",
  "female_cM",
  "female_r"
)

Family_r <- read.csv("data/Family_r.csv", header = TRUE)
eco <- read.csv("data/eco.csv", header = TRUE)

recom_time <- read.csv("data/CO_all_ecotype.csv", header = TRUE)

CO_sexdifferent <- read.table(
  "data/CO_sexdifferent.txt",
  header = FALSE,
  sep = "\t"
)

colnames(CO_sexdifferent) <- c(
  "Population",
  "Family",
  "Chr",
  "cM",
  "Sex"
)

# -----------------------------
# Chromosome length table
# -----------------------------

chr_length <- data.frame(
  Chr = c(
    "chr1", "chr10", "chr11", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr2", "chr20",
    "chr21", "chr3", "chr4", "chr5", "chr6", "chr7",
    "chr8", "chr9", "chrX", "chrY", "chrMT"
  ),
  Chr_length = c(
    30989408, 18147295, 19641151, 23954481, 18126596, 20609778,
    20106581, 21952404, 17740740, 21319505, 25303411, 21753483,
    16290672, 19092019, 34316924, 16393524, 21123412, 18824730,
    21585973, 22055034, 35320918, 46320123, 16580
  )
)

# -----------------------------
# Sex-specific recombination rate per chromosome
# -----------------------------

pop_eco <- Family_r %>%
  select(Population, Ecotype) %>%
  distinct()

r_chr <- recom_sexdifferent %>%
  filter(Population != "Marine_map", Population != "Pond_map") %>%
  left_join(pop_eco, by = "Population") %>%
  pivot_longer(
    cols = c(male_cM, male_r, female_cM, female_r),
    names_to = c("Sex", ".value"),
    names_pattern = "(female|male)_(cM|r)"
  ) %>%
  filter(!is.na(r), !is.na(cM), !is.na(Chr_length))

r_chr$Chr_length_scale <- as.numeric(scale(r_chr$Chr_length))

model_chr_r <- glmmTMB(
  r ~ (Sex + Ecotype + Chr_length_scale)^3 + (1 | Population),
  family = gaussian,
  data = r_chr,
  na.action = na.fail
)

model_set_chr_r <- dredge(model_chr_r)
model_chr_r_best <- get.models(model_set_chr_r, subset = 1)[[1]]

sink("results/sex_specific_chr_recombination_model.txt")

cat("Global model: chromosome recombination rate\n")
cat("===========================================\n")
print(summary(model_chr_r))

cat("\nBest model: chromosome recombination rate\n")
cat("=========================================\n")
print(summary(model_chr_r_best))

cat("\nR2\n")
print(r.squaredGLMM(model_chr_r_best))

sink()

write.csv(
  as.data.frame(model_set_chr_r),
  "results/model_selection_chr_recombination_rate.csv"
)

r_chr_eco <- r_chr %>%
  mutate(Group = interaction(Ecotype, Sex, sep = "_")) %>%
  group_by(Group, Chr) %>%
  summarise(
    mean.r = mean(r, na.rm = TRUE),
    mean.cM = mean(cM, na.rm = TRUE),
    mean.chrlength = mean(Chr_length, na.rm = TRUE),
    .groups = "drop"
  )

p_chr_r <- ggplot(
  r_chr_eco,
  aes(
    x = mean.chrlength / 1e6,
    y = mean.r,
    group = Group,
    color = Group,
    fill = Group
  )
) +
  geom_point(size = 0.6) +
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    linewidth = 0.5,
    alpha = 0.2
  ) +
  theme_classic(base_size = 10) +
  xlab("Chromosome length (Mb)") +
  ylab("Chromosome recombination rate (cM/Mb)") +
  scale_color_manual(
    values = c(
      "Marine_male" = "#528fad",
      "Marine_female" = "#528fad",
      "Freshwater_male" = "#ffd06f",
      "Freshwater_female" = "#ffd06f"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Marine_male" = "#528fad",
      "Marine_female" = "#528fad",
      "Freshwater_male" = "#ffd06f",
      "Freshwater_female" = "#ffd06f"
    )
  ) +
  guides(color = guide_legend(title = "Ecotype_Sex"))

ggsave(
  "figures/chromosome_recombination_rate_by_sex.png",
  p_chr_r,
  width = 5,
  height = 4,
  dpi = 300
)

p_chr_cM <- ggplot(
  r_chr_eco,
  aes(
    x = mean.chrlength / 1e6,
    y = mean.cM,
    group = Group,
    color = Group,
    fill = Group
  )
) +
  geom_point(size = 0.6) +
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    linewidth = 0.5,
    alpha = 0.2
  ) +
  theme_classic(base_size = 10) +
  xlab("Chromosome length (Mb)") +
  ylab("Genetic length (cM)") +
  scale_color_manual(
    values = c(
      "Marine_male" = "#528fad",
      "Marine_female" = "#528fad",
      "Freshwater_male" = "#ffd06f",
      "Freshwater_female" = "#ffd06f"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Marine_male" = "#528fad",
      "Marine_female" = "#528fad",
      "Freshwater_male" = "#ffd06f",
      "Freshwater_female" = "#ffd06f"
    )
  ) +
  guides(color = guide_legend(title = "Ecotype_Sex"))

ggsave(
  "figures/chromosome_genetic_length_by_sex.png",
  p_chr_cM,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# Crossover number per individual
# -----------------------------

recom_time$Population <- factor(
  recom_time$Population,
  levels = c("TVA", "POR", "UME", "RAA", "KRK", "PYO", "RYT", "BYN")
)

CO.indv <- recom_time %>%
  group_by(Sample_ID) %>%
  summarise(
    sum.CO = sum(CO, na.rm = TRUE),
    Family = first(Family),
    Population = first(Population),
    Ecotype = first(Ecotype),
    .groups = "drop"
  )

p_CO_indv <- CO.indv %>%
  filter(Population != "PYO") %>%
  ggplot(aes(x = Ecotype, y = sum.CO, fill = Ecotype)) +
  geom_boxplot(linewidth = 0.8, outlier.size = 0.5) +
  geom_jitter(size = 0.5, alpha = 0.5, width = 0.15) +
  labs(x = "Ecotype", y = "Number of crossovers") +
  scale_fill_manual(values = c("Marine" = "#528fad", "Freshwater" = "#ffd06f")) +
  theme_classic(base_size = 12)

ggsave(
  "figures/crossover_number_per_individual.png",
  p_CO_indv,
  width = 4,
  height = 4,
  dpi = 300
)

CO_indv_test <- CO.indv %>%
  t_test(sum.CO ~ Ecotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

write.csv(
  CO_indv_test,
  "results/crossover_number_per_individual_ttest.csv",
  row.names = FALSE
)

CO_family_test <- CO.indv %>%
  group_by(Family) %>%
  summarise(
    mean.CO = mean(sum.CO, na.rm = TRUE),
    Population = first(Population),
    Ecotype = first(Ecotype),
    .groups = "drop"
  ) %>%
  t_test(mean.CO ~ Ecotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

write.csv(
  CO_family_test,
  "results/crossover_number_family_mean_ttest.csv",
  row.names = FALSE
)

# -----------------------------
# Population-level crossover number
# -----------------------------

CO.indv_pop <- CO.indv %>%
  group_by(Population) %>%
  summarise(
    sd = sd(sum.CO, na.rm = TRUE),
    n = n(),
    avg = mean(sum.CO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se = sd / sqrt(n),
    lower.ci = avg - qt(1 - 0.05 / 2, n - 1) * se,
    upper.ci = avg + qt(1 - 0.05 / 2, n - 1) * se
  ) %>%
  left_join(eco, by = "Population")

write.csv(
  CO.indv_pop,
  "results/crossover_number_population_summary.csv",
  row.names = FALSE
)

p_CO_pop <- ggplot(
  CO.indv_pop,
  aes(
    x = avg,
    y = factor(
      Population,
      levels = c("RAA", "TVA", "POR", "UME", "KRK", "BYN", "RYT", "PYO")
    )
  )
) +
  geom_point(
    aes(color = Ecotype),
    size = 0.8,
    position = position_dodge(0.1)
  ) +
  geom_errorbar(
    aes(xmin = lower.ci, xmax = upper.ci, color = Ecotype),
    position = position_dodge(0.1),
    width = 0.02,
    alpha = 0.5
  ) +
  labs(x = "Number of crossovers", y = "Population") +
  theme_classic(base_size = 12) +
  guides(color = "none", fill = "none") +
  xlim(28, 42) +
  scale_color_manual(values = c("Freshwater" = "#cf8e2a", "Marine" = "#5c9dd5")) +
  scale_y_discrete(limits = rev)

ggsave(
  "figures/crossover_number_population.png",
  p_CO_pop,
  width = 4.5,
  height = 4,
  dpi = 300
)

# -----------------------------
# Mean crossover number per chromosome
# -----------------------------

chr_CO <- recom_time %>%
  group_by(Ecotype, Chr) %>%
  summarise(
    mean.CO = mean(CO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(chr_length, by = "Chr")

write.csv(
  chr_CO,
  "results/crossover_number_per_chromosome.csv",
  row.names = FALSE
)

p_chr_CO <- ggplot(
  chr_CO,
  aes(x = Chr_length / 1e6, y = mean.CO)
) +
  geom_point(aes(color = Ecotype), size = 0.6) +
  geom_smooth(
    aes(color = Ecotype, fill = Ecotype),
    method = "glm",
    linewidth = 0.5,
    alpha = 0.2
  ) +
  theme_classic(base_size = 10) +
  xlab("Chromosome length (Mb)") +
  ylab("Mean number of crossovers") +
  scale_color_manual(values = c("Marine" = "#528fad", "Freshwater" = "#ffd06f")) +
  scale_fill_manual(values = c("Marine" = "#528fad", "Freshwater" = "#ffd06f"))

ggsave(
  "figures/mean_crossover_number_per_chromosome.png",
  p_chr_CO,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# Sex-specific crossover number per chromosome
# -----------------------------

Offspring_count <- CO.indv %>%
  group_by(Ecotype, Population, Family) %>%
  summarise(
    offspring_num = n(),
    .groups = "drop"
  )

family_pop_eco <- Family_r %>%
  select(Family, Population, Ecotype) %>%
  distinct()

CO_sexdifferent <- CO_sexdifferent %>%
  group_by(Family, Sex) %>%
  count(Chr, name = "sum.CO") %>%
  left_join(family_pop_eco, by = "Family") %>%
  mutate(sum.CO = sum.CO - 1) %>%
  left_join(
    Offspring_count %>% select(Family, offspring_num),
    by = "Family"
  ) %>%
  mutate(avg.CO = sum.CO / offspring_num) %>%
  left_join(chr_length, by = "Chr") %>%
  filter(!is.na(avg.CO), !is.na(Chr_length))

sex_CO_summary <- CO_sexdifferent %>%
  group_by(Family, Sex) %>%
  summarise(
    sum = sum(avg.CO, na.rm = TRUE),
    Ecotype = first(Ecotype),
    .groups = "drop"
  ) %>%
  group_by(Sex, Ecotype) %>%
  summarise(
    mean = mean(sum, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  sex_CO_summary,
  "results/sex_specific_crossover_summary.csv",
  row.names = FALSE
)

perindv_chr_CO <- CO_sexdifferent %>%
  mutate(Group = interaction(Ecotype, Sex, sep = "_")) %>%
  group_by(Group, Chr) %>%
  summarise(
    mean.CO = mean(avg.CO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(chr_length, by = "Chr")

p_sex_chr_CO <- ggplot(
  perindv_chr_CO,
  aes(
    x = Chr_length / 1e6,
    y = mean.CO,
    group = Group,
    color = Group,
    fill = Group
  )
) +
  geom_point(size = 0.6) +
  geom_smooth(
    method = "lm",
    formula = y ~ log(x),
    linewidth = 0.5,
    alpha = 0.2
  ) +
  theme_classic(base_size = 10) +
  xlab("Chromosome length (Mb)") +
  ylab("Per-chromosome number of crossovers") +
  scale_color_manual(
    values = c(
      "Marine_male" = "#528fad",
      "Marine_female" = "#528fad",
      "Freshwater_male" = "#ffd06f",
      "Freshwater_female" = "#ffd06f"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Marine_male" = "#528fad",
      "Marine_female" = "#528fad",
      "Freshwater_male" = "#ffd06f",
      "Freshwater_female" = "#ffd06f"
    )
  ) +
  guides(color = guide_legend(title = "Ecotype_Sex"))

ggsave(
  "figures/sex_specific_crossover_number_per_chromosome.png",
  p_sex_chr_CO,
  width = 5,
  height = 4,
  dpi = 300
)

CO_sexdifferent$Chr_length_scale <- as.numeric(scale(CO_sexdifferent$Chr_length))

model_chr_co <- glmmTMB(
  avg.CO ~ (Sex + Ecotype + Chr_length_scale)^3 + (1 | Population),
  family = gaussian,
  data = CO_sexdifferent,
  na.action = na.fail
)

model_set_chr_co <- dredge(model_chr_co)
model_chr_co_best <- get.models(model_set_chr_co, subset = 1)[[1]]

sink("results/sex_specific_chr_crossover_model.txt")

cat("Global model: chromosome crossover number\n")
cat("=========================================\n")
print(summary(model_chr_co))

cat("\nBest model: chromosome crossover number\n")
cat("=======================================\n")
print(summary(model_chr_co_best))

cat("\nR2\n")
print(r.squaredGLMM(model_chr_co_best))

sink()

write.csv(
  as.data.frame(model_set_chr_co),
  "results/model_selection_chr_crossover_number.csv"
)

# -----------------------------
# Save session information
# -----------------------------

sink("results/sessionInfo_sex_specific_recombination_CO.txt")
sessionInfo()
sink()
