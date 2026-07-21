library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(readxl)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MuMIn)
library(car)
library(GGally)
library(lavaan)
library(AICcmodavg)

options(na.action = "na.fail")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# -----------------------------
# Basic settings
# -----------------------------

populations <- c("BYN", "KRK", "PYO", "POR", "RAA", "RYT", "TVA", "UME")

sex_ratios <- c(2, 1.18, 0.96, 0.739, 0.774, 0.862, 2.06, 1.23)

ecotypes <- c(
  "Freshwater", "Freshwater", "Freshwater",
  "Marine", "Marine", "Freshwater", "Marine", "Marine"
)

sex_ratio_data <- data.frame(
  Population = populations,
  Sex_Ratio = sex_ratios,
  Ecotype = ecotypes
)

eco <- data.frame(
  Population = c("TVA", "UME", "POR", "RAA", "PYO", "RYT", "BYN", "KRK"),
  Ecotype = rep(c("Marine", "Freshwater"), each = 4)
)

chr_order <- c(
  "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
  "chr8", "chr9", "chr10", "chr11", "chr13", "chr14",
  "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
  "chr21", "chrX"
)

pop_order <- c("TVA", "POR", "UME", "RAA", "KRK", "PYO", "RYT", "BYN")

pop_colors <- c(
  "#081d58", "#253494", "#225ea8", "#1d91c0",
  "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1"
)

names(pop_colors) <- pop_order

eco_colors <- c(
  "Marine" = "#528fad",
  "Freshwater" = "#ffd06f",
  "Pond" = "#ffd06f"
)

sex_colors <- c(
  "female" = "#e09351",
  "male" = "#94b594",
  "Female" = "#A94A4A",
  "Male" = "#547792"
)

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
  rhs <- ifelse(length(terms) == 0, "1", paste(terms, collapse = " + "))

  if (!is.null(random)) {
    rhs <- paste(rhs, random, sep = " + ")
  }

  as.formula(paste(response, "~", rhs))
}

partition_r2 <- function(data, response, main_effects, interactions = character(0), random = NULL) {
  main_model <- glmmTMB(
    make_formula(response, main_effects, random),
    family = gaussian,
    data = data,
    na.action = na.fail
  )

  r2_main <- get_r2(main_model)$R2m

  full_terms <- c(main_effects, interactions)

  full_model <- glmmTMB(
    make_formula(response, full_terms, random),
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

    reduced_model <- glmmTMB(
      make_formula(response, reduced_terms, random),
      family = gaussian,
      data = data,
      na.action = na.fail
    )

    r2_reduced <- get_r2(reduced_model)$R2m

    results <- rbind(
      results,
      data.frame(
        Term = term,
        R2_increase = max(0, r2_main - r2_reduced)
      )
    )
  }

  for (term in interactions) {
    interaction_model <- glmmTMB(
      make_formula(response, c(main_effects, term), random),
      family = gaussian,
      data = data,
      na.action = na.fail
    )

    r2_interaction <- get_r2(interaction_model)$R2m

    results <- rbind(
      results,
      data.frame(
        Term = term,
        R2_increase = max(0, r2_interaction - r2_main)
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

  results
}

multiplesheets <- function(fname) {
  sheets <- readxl::excel_sheets(fname)
  out <- lapply(sheets, function(x) as.data.frame(readxl::read_excel(fname, sheet = x)))
  names(out) <- sheets
  out
}

# -----------------------------
# 5 Mb window data
# -----------------------------

gene_data <- read.table(
  "data/gene_density_5Mb.txt",
  header = FALSE,
  sep = "\t"
)

cpg_data <- read.table(
  "data/cpg_5Mb.txt",
  header = FALSE,
  sep = "\t"
)

colnames(gene_data) <- c(
  "chr", "start_bp", "gene_end_bp", "number_gene",
  "gene_length", "window_length", "gene_density"
)

colnames(cpg_data) <- c(
  "chr", "start_bp", "cpg_end_bp", "cpg_content"
)

gene_data <- gene_data %>%
  filter(chr != "LG12")

cpg_data <- cpg_data %>%
  filter(chr != "LG12")

all_data <- list()

for (pop in populations) {
  recomb_data <- read.table(
    paste0("data/", pop, "_5Mb_5kb"),
    header = TRUE,
    sep = "\t"
  )

  pi_data <- read.table(
    paste0("data/", pop, "_pi_5Mb_5kb.windowed.pi"),
    header = TRUE,
    sep = "\t"
  )

  recomb_data <- recomb_data %>%
    filter(chr != "LG12")

  pi_data <- pi_data %>%
    filter(CHROM != "LG12")

  colnames(pi_data) <- c(
    "chr", "start_bp", "end_bp", "N_VARIANTS", "PI"
  )

  pi_data$start_bp <- pi_data$start_bp - 1

  recomb_data$chr <- as.character(recomb_data$chr)
  pi_data$chr <- as.character(pi_data$chr)

  merged_data <- recomb_data %>%
    left_join(pi_data, by = c("chr", "start_bp", "end_bp")) %>%
    left_join(gene_data, by = c("chr", "start_bp")) %>%
    left_join(cpg_data, by = c("chr", "start_bp"))

  merged_data$Population <- pop
  merged_data$sex_ratio <- sex_ratio_data$Sex_Ratio[sex_ratio_data$Population == pop]
  merged_data$Ecotype <- sex_ratio_data$Ecotype[sex_ratio_data$Population == pop]

  all_data[[pop]] <- merged_data
}

combined_data <- bind_rows(all_data)

combined_data$end_bp <- combined_data$gene_end_bp
combined_data$gene_end_bp <- NULL

combined_data <- combined_data %>%
  mutate(
    rho_scaled = as.numeric(scale(rho)),
    gene_density_scaled = as.numeric(scale(gene_density)),
    cpg_content_scaled = as.numeric(scale(cpg_content)),
    sex_ratio_scaled = as.numeric(scale(sex_ratio))
  )

write.csv(
  combined_data,
  "results/combined_5Mb_window_data.csv",
  row.names = FALSE
)

top_5_percent <- combined_data %>%
  group_by(Ecotype) %>%
  filter(rho >= quantile(rho, 0.95, na.rm = TRUE)) %>%
  ungroup()

low_5_percent <- combined_data %>%
  filter(rho != 0) %>%
  group_by(Ecotype) %>%
  filter(rho <= quantile(rho, 0.05, na.rm = TRUE)) %>%
  ungroup()

write.csv(top_5_percent, "results/top_5_percent_rho_5Mb.csv", row.names = FALSE)
write.csv(low_5_percent, "results/low_5_percent_rho_5Mb.csv", row.names = FALSE)

p_rho_hist <- ggplot(combined_data, aes(x = rho, fill = Ecotype)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_classic(base_size = 10) +
  scale_fill_manual(values = eco_colors) +
  xlab("Recombination rate") +
  ylab("Count")

ggsave(
  "figures/5Mb_rho_histogram.png",
  p_rho_hist,
  width = 5,
  height = 4,
  dpi = 300
)

p_top5 <- ggplot(top_5_percent, aes(x = Ecotype, y = rho, fill = Ecotype)) +
  geom_violin(outlier.shape = NA) +
  theme_classic(base_size = 10) +
  xlab("Ecotype") +
  ylab("Top 5% recombination rate") +
  stat_compare_means(method = "t.test", size = 3.5) +
  scale_fill_manual(values = eco_colors) +
  guides(fill = "none")

ggsave(
  "figures/5Mb_top5_rho_ecotype.png",
  p_top5,
  width = 4,
  height = 4,
  dpi = 300
)

p_low5 <- ggplot(low_5_percent, aes(x = Ecotype, y = rho, fill = Ecotype)) +
  geom_violin(outlier.shape = NA) +
  theme_classic(base_size = 10) +
  xlab("Ecotype") +
  ylab("Low 5% recombination rate") +
  stat_compare_means(method = "t.test", size = 3.5) +
  scale_fill_manual(values = eco_colors) +
  guides(fill = "none")

ggsave(
  "figures/5Mb_low5_rho_ecotype.png",
  p_low5,
  width = 4,
  height = 4,
  dpi = 300
)

model_5Mb_pi <- lmer(
  PI ~ rho * Ecotype + gene_density + cpg_content + sex_ratio + (1 | Population),
  data = combined_data
)

sink("results/model_5Mb_pi_summary.txt")
print(summary(model_5Mb_pi))
cat("\nCorrelation between CpG content and recombination rate\n")
print(cor.test(combined_data$cpg_content, combined_data$rho, method = "pearson"))
sink()

# -----------------------------
# Supplementary family-level recombination rate
# -----------------------------

supp <- multiplesheets("data/Supplementary_data.xlsx")
Family_r <- as.data.frame(supp[["Sheet1"]])

p_family_sex <- Family_r %>%
  filter(Ref == "PYOref") %>%
  ggplot(aes(x = Sex, y = ats_r, fill = Ecotype)) +
  geom_boxplot(linewidth = 0.8) +
  labs(x = "Sex", y = "Family-level recombination rate (cM/Mb)") +
  stat_compare_means(aes(group = Ecotype), label = "p.format", label.y = 6.6, size = 3) +
  stat_compare_means(aes(group = Sex), label = "p.format", label.y = 7, size = 3) +
  scale_fill_manual(values = eco_colors) +
  theme_classic(base_size = 12)

ggsave(
  "figures/family_recombination_rate_by_sex.png",
  p_family_sex,
  width = 4,
  height = 4,
  dpi = 300
)

p_family_age <- Family_r %>%
  filter(Ref == "TVAref") %>%
  ggplot(aes(x = age, y = ats_r, color = Sex, fill = Sex)) +
  geom_point(size = 0.8, position = position_dodge(0.1)) +
  geom_smooth(method = "glm", alpha = 0.2) +
  scale_color_manual(values = sex_colors, name = NULL) +
  scale_fill_manual(values = sex_colors, name = NULL) +
  labs(x = "Parental age", y = "Recombination rate (cM/Mb)") +
  theme_classic(base_size = 12) +
  guides(color = "none", fill = "none")

ggsave(
  "figures/family_recombination_rate_by_age.png",
  p_family_age,
  width = 4,
  height = 4,
  dpi = 300
)

family_test <- Family_r %>%
  group_by(Ref, Family) %>%
  summarise(
    mean_atsr = mean(ats_r, na.rm = TRUE),
    Ecotype = first(Ecotype),
    .groups = "drop"
  ) %>%
  t_test(mean_atsr ~ Ecotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

write.csv(
  family_test,
  "results/family_recombination_rate_ttest.csv",
  row.names = FALSE
)

# -----------------------------
# 1 Mb / 10 kb data
# -----------------------------

pi <- read.csv("data/pi_1mb10kb.txt", sep = "\t", header = FALSE)

colnames(pi) <- c(
  "Population", "Chr", "Start", "End",
  "Marker_number", "pi", "sex"
)

pi$Start <- pi$Start - 1
pi$Chr <- factor(pi$Chr, levels = chr_order)

rrate.eco <- read.csv(
  "data/ecotype_level_1mb10kb.txt",
  sep = "\t",
  header = FALSE
)

rrate.eco <- rrate.eco[, c(1, 3, 4, 5, 6, 7, 8)]

colnames(rrate.eco) <- c(
  "Ecotype", "Sex", "Chr", "Start",
  "End", "r", "Marker_number"
)

rrate.eco$Ecotype <- recode(rrate.eco$Ecotype, "Pond" = "Freshwater")
rrate.eco$Chr <- factor(rrate.eco$Chr, levels = chr_order)

rrate.eco <- rrate.eco %>%
  left_join(
    pi,
    by = c(
      "Ecotype" = "Population",
      "Chr",
      "Start",
      "End"
    )
  )

rrate.pop <- read.csv(
  "data/population_level_1mb10kb.txt",
  sep = "\t",
  header = FALSE
)

rrate.pop <- rrate.pop[, c(1, 3, 4, 5, 6, 7, 8)]

colnames(rrate.pop) <- c(
  "Population", "Sex", "Chr", "Start",
  "End", "r", "Marker_number"
)

rrate.pop$Population <- factor(rrate.pop$Population, levels = pop_order)
rrate.pop$Chr <- factor(rrate.pop$Chr, levels = chr_order)

cpg <- read.csv("data/m_1mb10kb.txt", sep = "\t", header = FALSE)
cpg <- cpg[, c(1, 2, 3, 7)]
colnames(cpg) <- c("Chr", "Start", "End", "cpg")
cpg$Chr <- factor(cpg$Chr, levels = chr_order)

gene <- read.csv("data/gene_density_1mb10kb.txt", sep = "\t", header = FALSE)
gene <- gene[, c(1, 2, 3, 7)]
colnames(gene) <- c("Chr", "Start", "End", "gene_density")
gene$Chr <- factor(gene$Chr, levels = chr_order)

rrate.eco <- rrate.eco %>%
  left_join(gene, by = c("Chr", "Start", "End")) %>%
  left_join(cpg, by = c("Chr", "Start", "End"))

rrate.pop <- rrate.pop %>%
  left_join(pi, by = c("Population", "Chr", "Start", "End")) %>%
  left_join(gene, by = c("Chr", "Start", "End")) %>%
  left_join(cpg, by = c("Chr", "Start", "End"))

write.csv(rrate.eco, "results/rrate_ecotype_1mb10kb_merged.csv", row.names = FALSE)
write.csv(rrate.pop, "results/rrate_population_1mb10kb_merged.csv", row.names = FALSE)

p_eco_r_pi <- rrate.eco %>%
  filter(Sex == "male", r < 20) %>%
  ggplot(aes(x = r, y = pi, color = Ecotype)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  labs(x = "Recombination rate (cM/Mb)", y = "Nucleotide diversity") +
  scale_color_manual(values = eco_colors)

ggsave(
  "figures/1mb10kb_ecotype_male_r_pi.png",
  p_eco_r_pi,
  width = 5,
  height = 4,
  dpi = 300
)

p_pop_r_pi <- rrate.pop %>%
  filter(Sex == "female", r < 30) %>%
  ggplot(aes(x = r, y = pi, color = Population)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  labs(x = "Recombination rate (cM/Mb)", y = "Nucleotide diversity") +
  scale_color_manual(values = pop_colors, name = NULL)

ggsave(
  "figures/1mb10kb_population_female_r_pi.png",
  p_pop_r_pi,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# 1 Mb / 200 kb data
# -----------------------------

pi.2 <- read.csv("data/pi_1mb200kb.txt", sep = "\t", header = FALSE)

colnames(pi.2) <- c(
  "Population", "Chr", "Start", "End",
  "Marker_number", "pi"
)

pi.2$Start <- pi.2$Start - 1
pi.2$Chr <- factor(pi.2$Chr, levels = chr_order)

rrate.2 <- read.csv(
  "data/r_perchr_1mb200kb.txt",
  sep = "\t",
  header = FALSE
)

colnames(rrate.2) <- c(
  "Population", "Sex", "Chr", "Start",
  "End", "r", "Marker_number"
)

rrate.2$Chr <- factor(rrate.2$Chr, levels = chr_order)

rrate.2 <- rrate.2 %>%
  left_join(pi.2, by = c("Population", "Chr", "Start", "End"))

cpg.2 <- read.csv("data/m_1mb200kb.txt", sep = "\t", header = FALSE)
cpg.2 <- cpg.2[, c(1, 2, 3, 4, 7)]
colnames(cpg.2) <- c("Chr", "Start", "End", "count_cpg", "cpg_score")
cpg.2$cpg <- cpg.2$count_cpg / 1000000
cpg.2$Chr <- factor(cpg.2$Chr, levels = chr_order)

gene.2 <- read.csv("data/gene_density_1mb200kb.txt", sep = "\t", header = FALSE)
gene.2 <- gene.2[, c(1, 2, 3, 4, 7)]
colnames(gene.2) <- c("Chr", "Start", "End", "count_gene", "gene_density")
gene.2$Chr <- factor(gene.2$Chr, levels = chr_order)

repeat.2 <- read.csv("data/repeat_1mb200kb.txt", sep = "\t", header = FALSE)
repeat.2 <- repeat.2[, c(1, 2, 3, 4, 7)]
colnames(repeat.2) <- c("Chr", "Start", "End", "count_rp", "repeat_density")
repeat.2$Chr <- factor(repeat.2$Chr, levels = chr_order)

dnm.2 <- read.csv("data/DNM_tvaref_1mb200kb.txt", sep = "\t", header = FALSE)
dnm.2 <- dnm.2[, c(1, 2, 3, 4, 7)]
colnames(dnm.2) <- c("Chr", "Start", "End", "count_dnm", "dnm_density")

rrate.2 <- rrate.2 %>%
  left_join(cpg.2, by = c("Chr", "Start", "End")) %>%
  left_join(gene.2, by = c("Chr", "Start", "End")) %>%
  left_join(repeat.2, by = c("Chr", "Start", "End")) %>%
  left_join(dnm.2, by = c("Chr", "Start", "End"))

write.csv(
  rrate.2,
  "results/rrate_1mb200kb_merged.csv",
  row.names = FALSE
)

# -----------------------------
# Correlation between ecotype maps and sex maps
# -----------------------------

marine_pond_corr <- rrate.2 %>%
  filter(
    repeat_density < 0.2,
    Population == "Marine"
  ) %>%
  select(Sex, Chr, Start, End, r) %>%
  left_join(
    rrate.2 %>%
      filter(Population == "Pond") %>%
      select(Sex, Chr, Start, End, r),
    by = c("Sex", "Chr", "Start", "End"),
    suffix = c(".marine", ".pond")
  ) %>%
  drop_na(r.marine, r.pond) %>%
  filter(r.marine <= 30, r.pond <= 30)

sex_corr <- rrate.2 %>%
  filter(Population %in% c("Marine", "Pond")) %>%
  select(Population, Sex, Chr, Start, End, r) %>%
  filter(Sex == "female") %>%
  left_join(
    rrate.2 %>%
      filter(Population %in% c("Marine", "Pond"), Sex == "male") %>%
      select(Population, Chr, Start, End, r),
    by = c("Population", "Chr", "Start", "End"),
    suffix = c(".female", ".male")
  ) %>%
  drop_na(r.female, r.male) %>%
  filter(r.female <= 30, r.male <= 30)

sink("results/recombination_rate_correlations.txt")

cat("Marine vs Freshwater recombination rate correlation\n")
print(cor.test(marine_pond_corr$r.marine, marine_pond_corr$r.pond, method = "pearson"))

cat("\nFemale vs male recombination rate correlation\n")
print(cor.test(sex_corr$r.female, sex_corr$r.male, method = "pearson"))

sink()

p_corr_eco <- ggplot(marine_pond_corr, aes(x = r.marine, y = r.pond)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_classic(base_size = 10) +
  labs(x = "Marine recombination rate", y = "Freshwater recombination rate")

ggsave(
  "figures/marine_freshwater_recombination_correlation.png",
  p_corr_eco,
  width = 4,
  height = 4,
  dpi = 300
)

p_corr_sex <- ggplot(sex_corr, aes(x = r.female, y = r.male)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_classic(base_size = 10) +
  labs(x = "Female recombination rate", y = "Male recombination rate")

ggsave(
  "figures/female_male_recombination_correlation.png",
  p_corr_sex,
  width = 4,
  height = 4,
  dpi = 300
)

# -----------------------------
# Window-level plots
# -----------------------------

p_r_cpg <- rrate.2 %>%
  filter(
    Population %in% c("Marine", "Pond"),
    repeat_density < 0.2,
    r < 25,
    r > 0
  ) %>%
  ggplot(aes(x = cpg, y = r, color = Sex)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  ylim(0, 20) +
  labs(x = "CpG content", y = "Recombination rate (cM/Mb)") +
  scale_color_manual(values = sex_colors)

ggsave(
  "figures/1mb200kb_r_vs_cpg_by_sex.png",
  p_r_cpg,
  width = 5,
  height = 4,
  dpi = 300
)

p_r_gene <- rrate.2 %>%
  filter(
    Population %in% c("Marine", "Pond"),
    repeat_density < 0.2,
    r < 25,
    r > 0
  ) %>%
  ggplot(aes(x = gene_density, y = r, color = Sex)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  ylim(0, 20) +
  labs(x = "Gene density", y = "Recombination rate (cM/Mb)") +
  scale_color_manual(values = sex_colors)

ggsave(
  "figures/1mb200kb_r_vs_gene_density_by_sex.png",
  p_r_gene,
  width = 5,
  height = 4,
  dpi = 300
)

p_pi_r_eco <- rrate.2 %>%
  filter(
    Population %in% c("Marine", "Pond"),
    repeat_density < 0.2,
    r < 25,
    r > 0
  ) %>%
  mutate(Ecotype = recode(Population, "Pond" = "Freshwater")) %>%
  ggplot(aes(x = r, y = pi, color = Ecotype)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  labs(x = "Recombination rate (cM/Mb)", y = "Nucleotide diversity") +
  scale_color_manual(values = eco_colors)

ggsave(
  "figures/1mb200kb_pi_vs_r_ecotype.png",
  p_pi_r_eco,
  width = 5,
  height = 4,
  dpi = 300
)

p_pi_cpg <- rrate.2 %>%
  filter(
    Population %in% c("Marine", "Pond"),
    repeat_density < 0.2
  ) %>%
  mutate(Ecotype = recode(Population, "Pond" = "Freshwater")) %>%
  ggplot(aes(x = cpg, y = pi, color = Ecotype)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  labs(x = "CpG content", y = "Nucleotide diversity") +
  scale_color_manual(values = eco_colors)

ggsave(
  "figures/1mb200kb_pi_vs_cpg_ecotype.png",
  p_pi_cpg,
  width = 5,
  height = 4,
  dpi = 300
)

p_pi_gene <- rrate.2 %>%
  filter(
    Population %in% c("Marine", "Pond"),
    repeat_density < 0.2
  ) %>%
  mutate(Ecotype = recode(Population, "Pond" = "Freshwater")) %>%
  ggplot(aes(x = gene_density, y = pi, color = Ecotype)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(method = "glm", alpha = 0.5) +
  theme_minimal() +
  labs(x = "Gene density", y = "Nucleotide diversity") +
  scale_color_manual(values = eco_colors)

ggsave(
  "figures/1mb200kb_pi_vs_gene_density_ecotype.png",
  p_pi_gene,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# Hot and cold windows
# -----------------------------

rrate.2.norepeat <- rrate.2 %>%
  filter(
    Population %in% c("Marine", "Pond"),
    repeat_density < 0.2
  ) %>%
  mutate(Ecotype = recode(Population, "Pond" = "Freshwater"))

top5.marine <- rrate.2.norepeat %>%
  filter(Population == "Marine") %>%
  filter(r >= quantile(r, 0.95, na.rm = TRUE))

top5.pond <- rrate.2.norepeat %>%
  filter(Population == "Pond") %>%
  filter(r >= quantile(r, 0.95, na.rm = TRUE))

lower5.marine <- rrate.2.norepeat %>%
  filter(Population == "Marine", r == 0)

lower5.pond <- rrate.2.norepeat %>%
  filter(Population == "Pond", r == 0)

top5.marine$Type <- "Hot"
top5.pond$Type <- "Hot"
lower5.marine$Type <- "Cold"
lower5.pond$Type <- "Cold"

hot_cold <- bind_rows(
  top5.marine,
  top5.pond,
  lower5.marine,
  lower5.pond
)

write.csv(
  hot_cold,
  "results/hot_cold_recombination_windows.csv",
  row.names = FALSE
)

sink("results/hot_cold_window_tests.txt")
cat("CpG content\n")
print(t.test(cpg ~ Type, hot_cold))
cat("\nGene density\n")
print(t.test(gene_density ~ Type, hot_cold))
sink()

p_hot_hist <- bind_rows(top5.marine, top5.pond) %>%
  mutate(Ecotype = recode(Population, "Pond" = "Freshwater")) %>%
  ggplot(aes(x = r, fill = Ecotype)) +
  geom_histogram(bins = 15, alpha = 0.6, position = "identity") +
  theme_classic(base_size = 10) +
  scale_fill_manual(values = eco_colors) +
  labs(x = "Recombination rate", y = "Count")

ggsave(
  "figures/hot_window_r_histogram.png",
  p_hot_hist,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# Population-level GLMMs for pi
# -----------------------------

rrate.2.pop <- rrate.2 %>%
  filter(Population %in% eco$Population) %>%
  left_join(eco, by = "Population") %>%
  filter(
    !is.na(pi),
    !is.na(Sex),
    repeat_density < 0.2
  ) %>%
  mutate(
    cpg_scale = as.numeric(scale(cpg)),
    gene_density_scale = as.numeric(scale(gene_density)),
    r_scale = as.numeric(scale(r)),
    pi_scale = as.numeric(scale(pi)),
    pi_log = log(pi)
  )

write.csv(
  rrate.2.pop,
  "results/rrate_1mb200kb_population_model_data.csv",
  row.names = FALSE
)

sink("results/population_window_correlations.txt")

cat("VIF for log(pi) model\n")
print(vif(lm(log(pi) ~ r + cpg + gene_density, data = rrate.2.pop)))

cat("\nSpearman correlations\n")
print(
  rrate.2.pop %>%
    select(pi, r, cpg, gene_density) %>%
    cor_test(method = "spearman") %>%
    adjust_pvalue(method = "bonferroni")
)

sink()

rrate.2.pop_avg <- rrate.2.pop %>%
  filter(Sex == "female") %>%
  select(
    Population, Ecotype, Chr, Start, End,
    pi, cpg, gene_density, repeat_density,
    cpg_scale, gene_density_scale, r
  ) %>%
  left_join(
    rrate.2.pop %>%
      filter(Sex == "male") %>%
      select(Population, Chr, Start, End, r),
    by = c("Population", "Chr", "Start", "End"),
    suffix = c(".female", ".male")
  ) %>%
  drop_na(r.female, r.male, pi, cpg, gene_density) %>%
  mutate(
    ravg = (r.female + r.male) / 2,
    ravg_scale = as.numeric(scale(ravg)),
    pi_scale = as.numeric(scale(pi)),
    pi_log = log(pi)
  )

write.csv(
  rrate.2.pop_avg,
  "results/rrate_1mb200kb_population_sexavg.csv",
  row.names = FALSE
)

model_B <- glmmTMB(
  pi_scale ~ (ravg_scale + cpg_scale + gene_density_scale + Ecotype)^2 +
    (1 | Population),
  family = gaussian,
  data = rrate.2.pop_avg,
  na.action = na.fail
)

model_set_B <- dredge(model_B)
best_model_pop <- get.models(model_set_B, subset = 1)[[1]]

sink("results/model_population_pi_sexavg.txt")
cat("Global model\n")
print(summary(model_B))
cat("\nBest model\n")
print(summary(best_model_pop))
cat("\nR2\n")
print(r.squaredGLMM(best_model_pop))
sink()

write.csv(
  as.data.frame(model_set_B),
  "results/model_selection_population_pi_sexavg.csv"
)

plot_residuals(
  best_model_pop,
  "figures/residuals_population_pi_sexavg.png",
  "Population sex-averaged pi model"
)

main_effects_pop <- c(
  "ravg_scale",
  "cpg_scale",
  "gene_density_scale",
  "Ecotype"
)

interactions_pop <- c(
  "cpg_scale:Ecotype",
  "cpg_scale:ravg_scale",
  "Ecotype:ravg_scale",
  "gene_density_scale:ravg_scale",
  "Ecotype:gene_density_scale"
)

r2_pop_partition <- partition_r2(
  data = rrate.2.pop_avg,
  response = "pi_scale",
  main_effects = main_effects_pop,
  interactions = interactions_pop,
  random = "(1 | Population)"
)

write.csv(
  r2_pop_partition,
  "results/r2_partition_population_pi_sexavg.csv",
  row.names = FALSE
)

p_pop_model <- rrate.2.pop_avg %>%
  ggplot(aes(x = ravg_scale, y = pi_scale, color = Ecotype, fill = Ecotype)) +
  geom_point(size = 0.8, alpha = 0.4) +
  geom_smooth(method = "glm", alpha = 0.2) +
  scale_fill_manual(values = eco_colors) +
  scale_color_manual(values = eco_colors) +
  labs(
    x = "Sex-averaged recombination rate",
    y = "Nucleotide diversity"
  ) +
  theme_classic(base_size = 12)

ggsave(
  "figures/population_pi_vs_sexavg_r.png",
  p_pop_model,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# Ecotype-level model
# -----------------------------

rrate.2.eco_avg <- rrate.2 %>%
  filter(Population %in% c("Marine", "Pond"), repeat_density < 0.2) %>%
  mutate(Ecotype = recode(Population, "Pond" = "Freshwater")) %>%
  filter(Sex == "female") %>%
  select(
    Ecotype, Chr, Start, End,
    pi, cpg, gene_density, repeat_density, r
  ) %>%
  left_join(
    rrate.2 %>%
      filter(Population %in% c("Marine", "Pond"), Sex == "male") %>%
      mutate(Ecotype = recode(Population, "Pond" = "Freshwater")) %>%
      select(Ecotype, Chr, Start, End, r),
    by = c("Ecotype", "Chr", "Start", "End"),
    suffix = c(".female", ".male")
  ) %>%
  drop_na(pi, cpg, gene_density, r.female, r.male) %>%
  mutate(
    ravg = (r.female + r.male) / 2,
    ravg_scale = as.numeric(scale(ravg)),
    cpg_scale = as.numeric(scale(cpg)),
    gene_density_scale = as.numeric(scale(gene_density)),
    pi_scale = as.numeric(scale(pi))
  )

write.csv(
  rrate.2.eco_avg,
  "results/rrate_1mb200kb_ecotype_sexavg.csv",
  row.names = FALSE
)

model_eco_pi <- glmmTMB(
  pi_scale ~ (ravg_scale + cpg_scale + gene_density_scale + Ecotype)^2,
  family = gaussian,
  data = rrate.2.eco_avg,
  na.action = na.fail
)

model_set_eco_pi <- dredge(model_eco_pi)
best_model_eco_pi <- get.models(model_set_eco_pi, subset = 1)[[1]]

sink("results/model_ecotype_pi_sexavg.txt")
cat("Global model\n")
print(summary(model_eco_pi))
cat("\nBest model\n")
print(summary(best_model_eco_pi))
cat("\nR2\n")
print(r.squaredGLMM(best_model_eco_pi))
sink()

write.csv(
  as.data.frame(model_set_eco_pi),
  "results/model_selection_ecotype_pi_sexavg.csv"
)

plot_residuals(
  best_model_eco_pi,
  "figures/residuals_ecotype_pi_sexavg.png",
  "Ecotype sex-averaged pi model"
)

# -----------------------------
# Zero and positive recombination models
# -----------------------------

rrate.2.pop <- rrate.2.pop %>%
  mutate(
    zero_indicator = ifelse(r == 0, 1, 0),
    r_positive = ifelse(r > 0, 1, 0)
  )

zero_model <- glmmTMB(
  zero_indicator ~ (Sex + cpg + gene_density)^3,
  family = binomial,
  data = rrate.2.pop,
  na.action = na.fail
)

positive_model <- glmmTMB(
  r_positive ~ (Sex + cpg + gene_density)^3 + (1 | Population),
  family = binomial,
  data = rrate.2.pop,
  na.action = na.fail
)

model_set_zero <- dredge(zero_model)
model_set_positive <- dredge(positive_model)

best_model_zero <- get.models(model_set_zero, subset = 1)[[1]]
best_model_positive <- get.models(model_set_positive, subset = 1)[[1]]

sink("results/model_recombination_zero_positive.txt")

cat("Zero model\n")
print(summary(best_model_zero))
cat("\nR2\n")
print(r.squaredGLMM(best_model_zero))

cat("\nPositive recombination model\n")
print(summary(best_model_positive))
cat("\nR2\n")
print(r.squaredGLMM(best_model_positive))

sink()

write.csv(
  as.data.frame(model_set_zero),
  "results/model_selection_recombination_zero.csv"
)

write.csv(
  as.data.frame(model_set_positive),
  "results/model_selection_recombination_positive.csv"
)

# -----------------------------
# Correlation matrix plot
# -----------------------------

rrate.2$r_log <- log(rrate.2$r * rrate.2$r + 1e-06)

p_pairs <- rrate.2 %>%
  filter(Population %in% c("Marine", "Pond")) %>%
  select(pi, r_log, cpg, gene_density) %>%
  ggpairs(
    lower = list(continuous = wrap("smooth", alpha = 0.05, size = 0.1)),
    diag = list(continuous = wrap("densityDiag")),
    upper = list(continuous = wrap("cor", size = 4))
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.background = element_blank()
  )

ggsave(
  "figures/pairwise_pi_r_cpg_gene_density.png",
  p_pairs,
  width = 7,
  height = 7,
  dpi = 300
)

# -----------------------------
# GWAS summary
# -----------------------------

if (file.exists("data/GWAS/CO.FarmCPU.csv")) {
  gwas <- read.csv("data/GWAS/CO.FarmCPU.csv")
  gwas2 <- read.csv("data/GWAS/CO.MLM.csv")
  gwas3 <- read.csv("data/GWAS/CO.GLM.csv")

  gwas$CHROM <- factor(gwas$CHROM, levels = chr_order)

  p_gwas <- gwas %>%
    filter(CO.FarmCPU < 0.01) %>%
    ggplot(aes(x = POS / 1e06, y = -log10(CO.FarmCPU), color = CHROM)) +
    geom_point(size = 0.5) +
    facet_grid(. ~ CHROM, space = "free_x", scales = "free_x") +
    theme_classic(base_size = 10) +
    xlab("Map position") +
    ylab("-log10(p-value)")

  ggsave(
    "figures/gwas_CO_FarmCPU.png",
    p_gwas,
    width = 8,
    height = 3,
    dpi = 300
  )

  common_gwas <- Reduce(
    intersect,
    list(
      subset(gwas, CO.FarmCPU < 1e-04)[, 1],
      subset(gwas2, CO.MLM < 1e-04)[, 1],
      subset(gwas3, CO.GLM < 1e-04)[, 1]
    )
  )

  write.csv(
    data.frame(SNP = common_gwas),
    "results/common_gwas_hits_1e-4.csv",
    row.names = FALSE
  )
}

# -----------------------------
# Linkage map plots
# -----------------------------

if (file.exists("data/8pops_ats_map_sexavg")) {
  pops_LM <- read.table(
    "data/8pops_ats_map_sexavg",
    header = FALSE,
    sep = "\t"
  )

  colnames(pops_LM) <- c("Population", "Chr", "Physical", "Genetic")

  pops_LM$Physical <- as.numeric(pops_LM$Physical)
  pops_LM$Chr <- factor(pops_LM$Chr, levels = chr_order)
  pops_LM$Population <- factor(
    pops_LM$Population,
    levels = c("RAA", "TVA", "POR", "UME", "KRK", "BYN", "RYT", "PYO")
  )

  p_linkage_pop <- ggplot(
    pops_LM,
    aes(
      x = Physical / 1e06,
      y = Genetic,
      color = Population,
      group = interaction(Population, Chr)
    )
  ) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ Chr, nrow = 5, ncol = 4, scales = "free_x") +
    theme_bw() +
    labs(
      x = "Physical position (Mb)",
      y = "Genetic position (cM)",
      color = "Population"
    ) +
    scale_color_manual(
      values = c(
        "#1C3664", "#016699", "#3499CC", "#6DC497",
        "#99CA3C", "#F79868", "#FFCD34", "#F5EC37"
      )
    )

  ggsave(
    "figures/linkage_map_population_sexavg.png",
    p_linkage_pop,
    width = 8,
    height = 6,
    dpi = 300
  )
}

if (file.exists("data/eco_sexdiff_map.txt")) {
  eco_LM <- read.table(
    "data/eco_sexdiff_map.txt",
    header = FALSE,
    sep = "\t"
  )

  colnames(eco_LM) <- c(
    "Ecotype", "Chr", "Physical",
    "Genetic_male", "Genetic_female"
  )

  eco_LM$Physical <- as.numeric(eco_LM$Physical)
  eco_LM$Chr <- factor(eco_LM$Chr, levels = chr_order)
  eco_LM$Ecotype <- recode(eco_LM$Ecotype, "Pond" = "Freshwater")

  eco_LM <- eco_LM %>%
    filter(!(Chr == "chrX" & Physical >= 18764539 & Physical <= 18858874)) %>%
    mutate(
      Genetic_female = if_else(
        Chr == "chrX" & Physical > 18875000 & Ecotype == "Freshwater",
        Genetic_female + 52.569,
        Genetic_female
      ),
      Genetic_female = if_else(
        Chr == "chrX" & Physical > 18875000 & Ecotype == "Marine",
        Genetic_female + 119.716,
        Genetic_female
      )
    )

  p_linkage_eco <- eco_LM %>%
    filter(Ecotype == "Marine") %>%
    ggplot(aes(group = interaction(Ecotype, Chr))) +
    geom_line(
      aes(x = Physical / 1e06, y = Genetic_female),
      color = "#e09351",
      linewidth = 0.4
    ) +
    geom_line(
      aes(x = Physical / 1e06, y = Genetic_male),
      color = "#94b594",
      linewidth = 0.4
    ) +
    facet_wrap(~ Chr, nrow = 3, ncol = 7, scales = "free_x") +
    theme_bw() +
    labs(
      x = "Physical position (Mb)",
      y = "Genetic position (cM)"
    )

  ggsave(
    "figures/linkage_map_ecotype_sexdiff_marine.png",
    p_linkage_eco,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# -----------------------------
# Path analysis
# -----------------------------

rrate.2.fw <- rrate.2.pop %>%
  filter(Ecotype == "Freshwater") %>%
  mutate(
    cpg_scale = as.numeric(scale(cpg)),
    r_scale = as.numeric(scale(r)),
    pi_scale = as.numeric(scale(pi))
  )

rrate.2.ma <- rrate.2.pop %>%
  filter(Ecotype == "Marine") %>%
  mutate(
    cpg_scale = as.numeric(scale(cpg)),
    r_scale = as.numeric(scale(r)),
    pi_scale = as.numeric(scale(pi))
  )

m.path.fw <- '
r_scale ~ a*cpg_scale
pi_scale ~ b*r_scale
pi_scale ~ i*cpg_scale
ab := a*b
total := i + (a*b)
'

m.path.ma <- '
cpg_scale ~ a*r_scale
pi_scale ~ b*cpg_scale
pi_scale ~ i*r_scale
ab := a*b
total := i + (a*b)
'

m.path.fit.fw <- sem(m.path.fw, data = rrate.2.fw)
m.path.fit.ma <- sem(m.path.ma, data = rrate.2.ma)

sink("results/path_analysis_summary.txt")

cat("Freshwater path model\n")
print(summary(m.path.fit.fw, standardized = TRUE, rsq = TRUE))
print(fitmeasures(m.path.fit.fw))

cat("\nMarine path model\n")
print(summary(m.path.fit.ma, standardized = TRUE, rsq = TRUE))
print(fitmeasures(m.path.fit.ma))

sink()

# -----------------------------
# Save session information
# -----------------------------

sink("results/sessionInfo_window_recombination_diversity.txt")
sessionInfo()
sink()
