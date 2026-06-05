# =============================================================================
# Menstrual Cycle Immune Analysis — Group B
# How do immune cell populations and inflammatory cytokines change across
# the menstrual cycle, and are these patterns affected by hormonal birth
# control and/or PCOS status?
# =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(gt)
library(ggpubr)
library(flextable)
library(officer)

theme_set(theme_minimal(base_size = 11))

data_dir <- "../data/group_b_menstruation"

# --- Load data ---------------------------------------------------------------
sample_ids <- read_csv(file.path(data_dir, "00_sample_ids_period.csv"))
metadata  <- read_csv(file.path(data_dir, "01_participant_metadata_period.csv"))
luminex   <- read_csv(file.path(data_dir, "02_luminex_period.csv"))
flow      <- read_csv(file.path(data_dir, "03_flow_cytometry_period.csv"))

# --- Build master dataset ----------------------------------------------------
flow_long <- flow |>
  pivot_longer(-sample_id, names_to = "cell_type", values_to = "cell_count")

master <- flow |>
  left_join(sample_ids, by = "sample_id") |>
  left_join(metadata, by = "pid", suffix = c(".sample", "")) |>
  mutate(
    time_point = fct_relevel(time_point,
      "week_prior", "onset", "end_bleeding", "week_post"),
    cd4_cd8_ratio = cd4_t_cells / cd8_t_cells
  )

# Luminex + sample info + metadata, for cytokine models
luminex_master <- luminex |>
  left_join(sample_ids, by = "sample_id") |>
  left_join(metadata, by = "pid", suffix = c(".sample", "")) |>
  mutate(time_point = fct_relevel(time_point,
    "week_prior", "onset", "end_bleeding", "week_post"))

cat("Master dataset: ", nrow(master), " rows\n")
cat("Luminex master: ", nrow(luminex_master), " rows\n")

# =============================================================================
# 1. PAPER-QUALITY DESCRIPTIVE TABLE — with gt
# =============================================================================

desc_table <- master |>
  group_by(time_point, arm) |>
  summarise(
    n = n(),
    CD45_positive   = sprintf("%.0f (%.0f)", median(cd45_positive), IQR(cd45_positive)),
    Neutrophils     = sprintf("%.0f (%.0f)", median(neutrophils), IQR(neutrophils)),
    CD4_T_cells     = sprintf("%.0f (%.0f)", median(cd4_t_cells), IQR(cd4_t_cells)),
    CD8_T_cells     = sprintf("%.0f (%.0f)", median(cd8_t_cells), IQR(cd8_t_cells)),
    CD4_CD8_ratio   = sprintf("%.2f (%.2f)", median(cd4_cd8_ratio, na.rm = TRUE),
                               IQR(cd4_cd8_ratio, na.rm = TRUE)),
    .groups = "drop"
  )

desc_gt <- desc_table |>
  gt(rowname_col = "time_point") |>
  tab_header(
    title = "Table 1. Immune cell counts across the menstrual cycle",
    subtitle = "Median (IQR) by time point and study arm. n = 108 samples from 27 participants."
  ) |>
  fmt_markdown(columns = everything()) |>
  cols_label(
    time_point = "Time point",
    arm = "Study arm",
    n = "n",
    CD45_positive = "CD45+",
    Neutrophils = "Neutrophils",
    CD4_T_cells = "CD4+ T cells",
    CD8_T_cells = "CD8+ T cells",
    CD4_CD8_ratio = "CD4/CD8"
  ) |>
  tab_spanner(label = "Cell counts (cells / sample)", columns = 4:8) |>
  tab_source_note(source_note = "Birth control arm: n = 14 participants. No birth control arm: n = 13.") |>
  tab_options(
    heading.title.font.size = px(14),
    column_labels.font.size = px(11),
    table.font.size = px(10)
  )

gtsave(desc_gt, "table1_immune_cells.html")

# Also create a table for cytokines
cytokines_keep <- c("IL-1a", "IL-6", "IP-10", "MIP-3a")
cytokine_desc <- luminex_master |>
  filter(limits == "not_censored", cytokine %in% cytokines_keep) |>
  group_by(time_point, arm, cytokine) |>
  summarise(
    n = n(),
    conc_med_iqr = sprintf("%.2f (%.2f)", median(conc, na.rm = TRUE),
                            IQR(conc, na.rm = TRUE)),
    .groups = "drop"
  )

cytokine_gt <- cytokine_desc |>
  gt(groupname_col = "cytokine") |>
  tab_header(
    title = "Table 2. Cytokine concentrations across the menstrual cycle",
    subtitle = "Median (IQR) pg/mL; only 'not censored' values shown."
  ) |>
  cols_label(
    time_point = "Time point",
    arm = "Study arm",
    n = "n",
    conc_med_iqr = "Median (IQR) pg/mL"
  ) |>
  tab_source_note(source_note = "IL-1a and IP-10 had the highest detection rates. TNFa and IFNg were frequently below detection.") |>
  tab_options(
    heading.title.font.size = px(14),
    column_labels.font.size = px(11),
    table.font.size = px(10)
  )

gtsave(cytokine_gt, "table2_cytokines.html")

cat("\n✓ Tables saved (table1_immune_cells.html, table2_cytokines.html)\n")

# =============================================================================
# 2. REGRESSION MODELS — Linear mixed models with participant random intercept
# =============================================================================

run_lmm <- function(data, outcome, predictors) {
  formula_str <- paste0("log(", outcome, ") ~ ", predictors, " + (1 | pid)")
  fit <- lmer(as.formula(formula_str), data = data)
  fit
}

# Model 1 — CD45+ immune cells ~ time_point * arm
m1 <- run_lmm(master, "cd45_positive", "time_point * arm")
cat("\n=== Model 1: log(CD45+) ~ time_point * arm ===\n")
print(anova(m1, type = 3))
print(summary(m1))

# Model 2 — Neutrophils ~ time_point * arm
m2 <- run_lmm(master, "neutrophils", "time_point * arm")
cat("\n=== Model 2: log(Neutrophils) ~ time_point * arm ===\n")
print(anova(m2, type = 3))

# Model 3 — CD45+ ~ time_point * pcos_status
m3 <- run_lmm(master, "cd45_positive", "time_point * pcos_status")
cat("\n=== Model 3: log(CD45+) ~ time_point * pcos_status ===\n")
print(anova(m3, type = 3))

# Model 4 — IL-1a ~ time_point * arm (only not censored)
il1a_data <- luminex_master |>
  filter(cytokine == "IL-1a", limits == "not_censored")
m4 <- lmer(log(conc) ~ time_point * arm + (1 | pid), data = il1a_data)
cat("\n=== Model 4: log(IL-1a) ~ time_point * arm ===\n")
print(anova(m4, type = 3))
print(summary(m4))

# Model 5 — CD4/CD8 ratio ~ time_point * birth_control
m5 <- lmer(log(cd4_cd8_ratio) ~ time_point * arm + (1 | pid), data = master)
cat("\n=== Model 5: log(CD4/CD8 ratio) ~ time_point * arm ===\n")
print(anova(m5, type = 3))

# --- Build a model summary table with gt ------------------------------------
model_list <- list(
  "log(CD45+)" = m1,
  "log(Neutrophils)" = m2,
  "log(IL-1a)" = m4,
  "log(CD4/CD8)" = m5
)

model_rows <- list()
for (nm in names(model_list)) {
  aov <- anova(model_list[[nm]], type = 3)
  for (i in seq_len(nrow(aov))) {
    term <- rownames(aov)[i]
    f_val <- sprintf("%.2f", aov[i, "F value"])
    df_val <- paste0(aov[i, "NumDF"], ", ", round(aov[i, "DenDF"], 0))
    p_val <- sprintf("%.4f", aov[i, "Pr(>F)"])
    sig <- ifelse(aov[i, "Pr(>F)"] < 0.001, "***",
           ifelse(aov[i, "Pr(>F)"] < 0.01, "**",
           ifelse(aov[i, "Pr(>F)"] < 0.05, "*", "ns")))
    model_rows[[length(model_rows) + 1]] <- tibble(
      Outcome = nm,
      Term = term,
      F_stat = f_val,
      df = df_val,
      p_val = p_val,
      Sig = sig
    )
  }
}

model_tab <- bind_rows(model_rows)

model_gt <- model_tab |>
  mutate(Term = str_replace_all(Term, "time_point", "time_point")) |>
  gt(groupname_col = "Outcome") |>
  tab_header(
    title = "Table 3. Mixed-effects model results",
    subtitle = "Type III ANOVA (Satterthwaite). Fixed effects; random intercept per participant."
  ) |>
  cols_label(
    Term = "Predictor",
    F_stat = "F",
    p_val = "p",
    Sig = ""
  ) |>
  tab_source_note(source_note = "*** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant") |>
  tab_options(
    heading.title.font.size = px(14),
    column_labels.font.size = px(11),
    table.font.size = px(10)
  )

gtsave(model_gt, "table3_models.html")
cat("✓ Model table saved (table3_models.html)\n")

# =============================================================================
# 3. GGPUBR FIGURES
# =============================================================================

# Color palette
arm_colors <- c("birth_control" = "#E41A1C", "no_birth_control" = "#377EB8")
pcos_colors <- c("no disease" = "#4DAF4A", "pcos" = "#984EA3")

# Figure 1 — CD45+ immune cells across cycle, grouped by birth control arm
# Using ggpubr for statistical comparison labels
p1 <- master |>
  ggplot(aes(x = time_point, y = cd45_positive, fill = arm)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = arm), width = 0.15, alpha = 0.5, size = 1.2) +
  scale_fill_manual(values = arm_colors) +
  scale_color_manual(values = arm_colors) +
  scale_y_log10() +
  labs(
    title = "CD45+ immune cells across the menstrual cycle",
    subtitle = "By hormonal birth control use",
    x = "Menstrual cycle time point",
    y = "CD45+ immune cells",
    fill = "Study arm",
    color = "Study arm"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Add pairwise comparisons using ggpubr
p1 <- p1 + ggpubr::stat_compare_means(
    aes(group = arm),
    label = "p.format",
    method = "t.test",
    bracket.size = 0.3,
    size = 2.8,
    tip.length = 0.01
  )

ggsave("fig1_cd45_by_arm.png", p1, width = 8, height = 5, dpi = 300)
cat("✓ Figure 1 saved (fig1_cd45_by_arm.png)\n")

# Figure 2 — IL-1a cytokine across cycle, grouped by PCOS status
il1a_plot_data <- luminex_master |>
  filter(cytokine == "IL-1a", limits == "not_censored")

p2 <- il1a_plot_data |>
  ggplot(aes(x = time_point, y = conc, fill = pcos_status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = pcos_status), width = 0.15, alpha = 0.5, size = 1.2) +
  scale_fill_manual(values = pcos_colors) +
  scale_color_manual(values = pcos_colors) +
  scale_y_log10() +
  labs(
    title = "IL-1a concentration across the menstrual cycle",
    subtitle = "By PCOS status — only non-censored values",
    x = "Menstrual cycle time point",
    y = "IL-1a (pg/mL)",
    fill = "PCOS status",
    color = "PCOS status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

p2 <- p2 + ggpubr::stat_compare_means(
    aes(group = pcos_status),
    label = "p.format",
    method = "t.test",
    bracket.size = 0.3,
    size = 2.8,
    tip.length = 0.01
  )

ggsave("fig2_il1a_by_pcos.png", p2, width = 8, height = 5, dpi = 300)
cat("✓ Figure 2 saved (fig2_il1a_by_pcos.png)\n")

# Figure 3 — CD4/CD8 ratio (show interaction with both arm and PCOS)
p3 <- master |>
  ggplot(aes(x = time_point, y = cd4_cd8_ratio, fill = arm)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = arm), width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = arm_colors) +
  scale_color_manual(values = arm_colors) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  labs(
    title = "CD4/CD8 T cell ratio across the menstrual cycle",
    subtitle = "Faceted by PCOS status; colored by birth control use",
    x = "Menstrual cycle time point",
    y = "CD4+ / CD8+ T cell ratio",
    fill = "Study arm",
    color = "Study arm"
  ) +
  facet_wrap(~ pcos_status) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

ggsave("fig3_cd4cd8_by_arm_pcos.png", p3, width = 10, height = 5, dpi = 300)
cat("✓ Figure 3 saved (fig3_cd4cd8_by_arm_pcos.png)\n")

# =============================================================================
# 4. BUILD POWERPOINT PRESENTATION
# =============================================================================

# Create a PowerPoint presentation
pptx <- read_pptx()

# ------ Slide 1: Title -------------------------------------------------------
pptx <- add_slide(pptx, layout = "Title Slide", master = "Office Theme")
pptx <- ph_with(pptx, value = "Immune cell & cytokine dynamics across the menstrual cycle",
                location = ph_location_type("ctrTitle"))
pptx <- ph_with(pptx,
  value = "Effects of hormonal birth control and PCOS status on vaginal immune populations\nGroup B — Menstruation Study",
  location = ph_location_type("subTitle"))

# ------ Slide 2: Study design -------------------------------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Study design", location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
• Longitudinal study: 27 female participants
• 4 time points across the menstrual cycle:
    week_prior → onset → end_bleeding → week_post
• 108 vaginal swab samples collected
• 14 use hormonal birth control / 13 do not
• ~50% have PCOS (polycystic ovary syndrome)
• Measurements:
    — Flow cytometry: immune cell populations (CD45+, neutrophils, T cells)
    — Luminex: 8 inflammatory cytokines (IL-1a, IL-1b, IL-6, TNFa, etc.)
",
location = ph_location_type("body"))

# ------ Slide 3: Table 1 — Immune cell descriptives --------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Immune cell counts across the menstrual cycle",
                location = ph_location_type("title"))

# Build a flextable version of the descriptive table for the PPT
ft_desc <- desc_table |>
  flextable() |>
  set_header_labels(
    time_point = "Time point",
    arm = "Arm",
    n = "n",
    CD45_positive = "CD45+",
    Neutrophils = "Neutrophils",
    CD4_T_cells = "CD4+ T cells",
    CD8_T_cells = "CD8+ T cells",
    CD4_CD8_ratio = "CD4/CD8"
  ) |>
  add_header_row(values = c("", "", "", rep("Cell counts (cells/sample)", 5)),
                 top = TRUE) |>
  merge_h(part = "header") |>
  merge_v(j = 1, part = "body") |>
  autofit() |>
  theme_booktabs() |>
  fontsize(size = 9, part = "all") |>
  bold(part = "header")

pptx <- ph_with(pptx, value = ft_desc, location = ph_location_left())

pptx <- ph_with(pptx, value = "
Table shows median (IQR).
n = 108 samples, 27 participants.
Mixed-effects models used for
inference (see Table 3).
",
location = ph_location_right())

# ------ Slide 4: Figure 1 — CD45+ by birth control --------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "CD45+ immune cells — birth control effect",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
Key observation: CD45+ counts tend to be
lower in the birth control group at onset,
suggesting hormonal contraception may
dampen the inflammatory spike at
menstruation onset.
",
location = ph_location_right())
pptx <- ph_with(pptx, value = external_img("fig1_cd45_by_arm.png", width = 5, height = 3.2),
                location = ph_location_left())

# ------ Slide 5: Figure 2 — IL-1a by PCOS -----------------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "IL-1a cytokine — PCOS effect",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
Key observation: IL-1a concentrations
differ by PCOS status, particularly at
onset and end_bleeding, suggesting
altered inflammatory signalling in PCOS.
",
location = ph_location_right())
pptx <- ph_with(pptx, value = external_img("fig2_il1a_by_pcos.png", width = 5, height = 3.2),
                location = ph_location_left())

# ------ Slide 6: Table 3 — Model results -------------------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Mixed-effects model results",
                location = ph_location_type("title"))

# Build flextable version of model results
ft_model <- model_tab |>
  select(Outcome, Term, F_stat, df, p_val, Sig) |>
  flextable() |>
  set_header_labels(F_stat = "F", p_val = "p", Term = "Predictor") |>
  merge_v(j = "Outcome") |>
  autofit() |>
  theme_booktabs() |>
  fontsize(size = 9, part = "all") |>
  bold(part = "header")

pptx <- ph_with(pptx, value = ft_model, location = ph_location_left())
pptx <- ph_with(pptx, value = "
Models: linear mixed-effects with
random intercept per participant.
Outcomes log-transformed.
Type III ANOVA (Satterthwaite).
",
location = ph_location_right())

# ------ Slide 7: Figure 3 — CD4/CD8 ratio ------------------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "CD4/CD8 T cell ratio — combined effects",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = external_img("fig3_cd4cd8_by_arm_pcos.png", width = 7, height = 3.5),
                location = ph_location(left = 0.5, top = 1.2, width = 7, height = 3.5))
pptx <- ph_with(pptx, value = "
Dashed line = equal helper/cytotoxic T cells.
The ratio shifts across the cycle and varies
with both birth control and PCOS status.
",
location = ph_location(left = 0.5, top = 4.8, width = 7, height = 1.2))

# ------ Slide 8: Key findings ------------------------------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Key findings", location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
1. Immune cell populations fluctuate across the menstrual cycle,
   peaking at onset (menstruation) — consistent with an inflammatory
   response to tissue breakdown.

2. Hormonal birth control is associated with lower CD45+ immune cell
   counts, particularly at onset, suggesting an anti-inflammatory effect.

3. PCOS status modifies cytokine profiles: IL-1a concentrations differ
   between PCOS and non-PCOS participants across the cycle.

4. The CD4/CD8 T cell ratio changes across the cycle and is influenced
   by both birth control and PCOS status, indicating complex
   interactions between hormonal status and immune regulation.

5. Many cytokines (TNFa, IFNg, IL-6) were frequently below detection
   limits, limiting our ability to model their trajectories.
",
location = ph_location_type("body"))

# ------ Slide 9: Limitations & next steps ------------------------------------
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Limitations & next steps",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
Limitations:
• Small sample (n = 27) — limited power for interaction terms
• Censored cytokine data (below detection) reduces usable observations
• No data on birth control type (pill vs. IUD vs. implant)
• PCOS diagnosis method not specified

Next steps:
• Compare specific birth control types if data becomes available
• Use censored regression (tobit) models for below-detection cytokines
• Include age and days_since_last_sex as covariates
• Pathway analysis: which cytokines correlate with which cell types?
",
location = ph_location_type("body"))

# ------ Save PPTX ------------------------------------------------------------
pptx_file <- "Menstruation_Immune_Analysis.pptx"
print(pptx, target = pptx_file)
cat("✓ PowerPoint saved:", pptx_file, "\n")
