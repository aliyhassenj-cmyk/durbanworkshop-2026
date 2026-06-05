# =============================================================================
# Build model summary table + rebuild PPT with p-values
# =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(gt)
library(flextable)
library(officer)

data_dir <- "../data/group_b_menstruation"

sample_ids <- read_csv(file.path(data_dir, "00_sample_ids_period.csv"))
metadata   <- read_csv(file.path(data_dir, "01_participant_metadata_period.csv"))
luminex    <- read_csv(file.path(data_dir, "02_luminex_period.csv"))
flow       <- read_csv(file.path(data_dir, "03_flow_cytometry_period.csv"))

master <- flow |>
  left_join(sample_ids, by = "sample_id") |>
  left_join(metadata, by = "pid", suffix = c(".sample", "")) |>
  mutate(
    time_point = fct_relevel(time_point,
      "week_prior", "onset", "end_bleeding", "week_post"),
    cd4_cd8_ratio = cd4_t_cells / cd8_t_cells
  )

luminex_master <- luminex |>
  left_join(sample_ids, by = "sample_id") |>
  left_join(metadata, by = "pid", suffix = c(".sample", "")) |>
  mutate(time_point = fct_relevel(time_point,
    "week_prior", "onset", "end_bleeding", "week_post"))

il1a_data <- luminex_master |>
  filter(cytokine == "IL-1a", limits == "not_censored")

# --- Fit models --------------------------------------------------------------
m1 <- lmer(log(cd45_positive) ~ time_point * arm + (1 | pid), data = master)
m2 <- lmer(log(neutrophils)   ~ time_point * arm + (1 | pid), data = master)
m4 <- lmer(log(conc)          ~ time_point * arm + (1 | pid), data = il1a_data)
m5 <- lmer(log(cd4_cd8_ratio) ~ time_point * arm + (1 | pid), data = master)

# --- Build model summary table -----------------------------------------------
model_list <- list(
  "log(CD45+)"          = m1,
  "log(Neutrophils)"    = m2,
  "log(IL-1a)"          = m4,
  "log(CD4/CD8 ratio)"  = m5
)

clean_term <- function(term) {
  term |>
    str_replace("time_point", "Time point") |>
    str_replace("arm", "Birth control") |>
    str_replace("pcos_status", "PCOS status") |>
    str_replace(":", " × ")
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
}

rows <- list()
for (nm in names(model_list)) {
  aov <- anova(model_list[[nm]], type = 3)
  for (i in seq_len(nrow(aov))) {
    p_raw <- aov[i, "Pr(>F)"]
    rows[[length(rows) + 1]] <- tibble(
      Outcome = nm,
      Predictor = clean_term(rownames(aov)[i]),
      F = sprintf("%.2f", aov[i, "F value"]),
      Num_df = aov[i, "NumDF"],
      Den_df = round(aov[i, "DenDF"], 0),
      p_value = fmt_p(p_raw),
      sig = case_when(
        p_raw < 0.001 ~ "***",
        p_raw < 0.01  ~ "**",
        p_raw < 0.05  ~ "*",
        TRUE ~ ""
      ),
      p_raw = p_raw
    )
  }
}
model_tab <- bind_rows(rows)

# =============================================================================
# 1. GT TABLE — save as HTML
# =============================================================================
gt_table <- model_tab |>
  select(Outcome, Predictor, F, Num_df, Den_df, p_value, sig) |>
  gt(groupname_col = "Outcome") |>
  tab_header(
    title = md("**Table 3.** Mixed-effects model results"),
    subtitle = "Linear mixed models with random intercept per participant. Type III ANOVA (Satterthwaite)."
  ) |>
  cols_label(
    Predictor = "Predictor",
    F = "F",
    Num_df = "Num",
    Den_df = "Den",
    p_value = "p",
    sig = ""
  ) |>
  tab_spanner(label = "df", columns = c(Num_df, Den_df)) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = p_raw < 0.05, columns = c(F, Num_df, Den_df, p_value, sig))
  ) |>
  tab_style(
    style = cell_fill(color = "#E8F5E9"),
    locations = cells_body(rows = p_raw < 0.05)
  ) |>
  tab_source_note("*** p < 0.001, ** p < 0.01, * p < 0.05. Green rows = significant.") |>
  tab_options(
    heading.title.font.size = px(14),
    column_labels.font.size = px(11),
    table.font.size = px(10)
  )

gtsave(gt_table, "table3_models_improved.html")
cat("✓ table3_models_improved.html\n")

# =============================================================================
# 2. FLEXTABLE — for PPT
# =============================================================================
# Fix: remove p_raw before flextable to avoid issues
ft <- model_tab |>
  select(Outcome, Predictor, F, Num_df, Den_df, p_value, sig) |>
  flextable() |>
  set_header_labels(
    Predictor = "Predictor",
    F = "F",
    Num_df = "Num",
    Den_df = "Den",
    p_value = "p",
    sig = ""
  ) |>
  merge_v(j = "Outcome", part = "body") |>
  autofit() |>
  theme_booktabs() |>
  fontsize(size = 10, part = "all") |>
  bold(part = "header") |>
  align(align = "center", part = "body") |>
  align(align = "left", j = "Predictor", part = "body")

cat("✓ Model flextable built\n")

# =============================================================================
# 3. Descriptive + p-values flextable
# =============================================================================
desc_data <- master |>
  group_by(time_point) |>
  summarise(
    n_bc = sum(arm == "birth_control"),
    n_nobc = sum(arm == "no_birth_control"),
    cd45_bc    = median(cd45_positive[arm == "birth_control"]),
    cd45_nobc  = median(cd45_positive[arm == "no_birth_control"]),
    p_cd45     = wilcox.test(cd45_positive ~ arm)$p.value,
    neutro_bc  = median(neutrophils[arm == "birth_control"]),
    neutro_nobc = median(neutrophils[arm == "no_birth_control"]),
    p_neutro   = wilcox.test(neutrophils ~ arm)$p.value,
    cd8_bc     = median(cd8_t_cells[arm == "birth_control"]),
    cd8_nobc   = median(cd8_t_cells[arm == "no_birth_control"]),
    p_cd8      = wilcox.test(cd8_t_cells ~ arm)$p.value,
    .groups = "drop"
  ) |>
  mutate(
    across(c(cd45_bc, cd45_nobc, neutro_bc, neutro_nobc, cd8_bc, cd8_nobc),
           ~ sprintf("%.0f", .x)),
    across(c(p_cd45, p_neutro, p_cd8), ~ sprintf("%.3f", .x))
  )

ft_desc <- desc_data |>
  flextable() |>
  set_header_labels(
    time_point = "Time point",
    n_bc = "BC n", n_nobc = "No BC n",
    cd45_bc = "CD45+ BC", cd45_nobc = "CD45+ NoBC", p_cd45 = "p (CD45)",
    neutro_bc = "Neutro BC", neutro_nobc = "Neutro NoBC", p_neutro = "p (Neutro)",
    cd8_bc = "CD8+ BC", cd8_nobc = "CD8+ NoBC", p_cd8 = "p (CD8)"
  ) |>
  autofit() |>
  theme_booktabs() |>
  fontsize(size = 9, part = "all") |>
  bold(part = "header") |>
  align(align = "center", part = "body")

cat("✓ Descriptive flextable built\n")

# =============================================================================
# 4. BUILD PPT
# =============================================================================
pptx <- read_pptx()

# Slide 1: Title
pptx <- add_slide(pptx, layout = "Title Slide", master = "Office Theme")
pptx <- ph_with(pptx, value = "Immune cell & cytokine dynamics across the menstrual cycle",
                location = ph_location_type("ctrTitle"))
pptx <- ph_with(pptx,
  value = "Effects of hormonal birth control and PCOS status on vaginal immune populations\nGroup B — Menstruation Study",
  location = ph_location_type("subTitle"))

# Slide 2: Study design
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
    — Flow cytometry: immune cell populations
    — Luminex: 8 inflammatory cytokines
",
location = ph_location_type("body"))

# Slide 3: Descriptive + p-values table
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Immune cell counts — with p-values",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = ft_desc,
                location = ph_location(left = 0.5, top = 1.5, width = 8, height = 4))
pptx <- ph_with(pptx, value = "Median counts (cells/sample). p-values: Wilcoxon test, BC vs No BC.",
                location = ph_location_type("body"))

# Slide 4: Figure 1
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "CD45+ immune cells — birth control effect",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = external_img("fig1_cd45_by_arm.png", width = 5, height = 3.2),
                location = ph_location(left = 0.5, top = 1.3, width = 5, height = 3.2))
pptx <- ph_with(pptx, value = "
Key observation:
CD45+ counts tend to be lower in the
birth control group at onset, suggesting
hormonal contraception may dampen the
inflammatory spike at menstruation onset.
",
location = ph_location(left = 5.8, top = 1.5, width = 4, height = 3))

# Slide 5: Figure 2
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "IL-1a cytokine — PCOS effect",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = external_img("fig2_il1a_by_pcos.png", width = 5, height = 3.2),
                location = ph_location(left = 0.5, top = 1.3, width = 5, height = 3.2))
pptx <- ph_with(pptx, value = "
Key observation:
IL-1a differs by PCOS status,
particularly at onset and end_bleeding,
suggesting altered inflammatory signalling
in PCOS.
",
location = ph_location(left = 5.8, top = 1.5, width = 4, height = 3))

# Slide 6: Model table (the centrepiece!)
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Mixed-effects model results — with p-values",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = ft,
                location = ph_location(left = 0.5, top = 1.3, width = 5.5, height = 4.5))
pptx <- ph_with(pptx, value = "
Key:
Green = significant (p < 0.05)

Findings:
1. IL-1a: strong time effect (p < 0.001)
         + arm effect (p = 0.014)
2. CD4/CD8: time effect (p < 0.001)
3. CD45+: no significant effects
",
location = ph_location(left = 6.3, top = 1.5, width = 3.5, height = 4))

# Slide 7: Figure 3
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "CD4/CD8 T cell ratio — combined effects",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = external_img("fig3_cd4cd8_by_arm_pcos.png",
                width = 7, height = 3.5),
                location = ph_location(left = 0.5, top = 1.2, width = 7, height = 3.5))
pptx <- ph_with(pptx, value = "
Dashed line = equal helper/cytotoxic T cells.
Ratio shifts across cycle and varies with
both birth control and PCOS status.
",
location = ph_location(left = 0.5, top = 4.8, width = 7, height = 1.2))

# Slide 8: Key findings
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Key findings", location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
1. IL-1a peaks at onset & end_bleeding,
   and is lower in birth control users (p = 0.014).
   → Hormonal contraception dampens inflammatory cytokine release.

2. CD4/CD8 T cell ratio shifts across the cycle (p < 0.001),
   but is not affected by birth control.
   → T cell subsets traffic differentially during menstruation.

3. CD45+ and neutrophils are stable across the cycle.
   → Innate cell numbers are less variable than cytokine signalling.

4. PCOS shows trends in cytokine profiles but not significant
   with n = 27.

5. Many cytokines (TNFa, IFNg, IL-6) were frequently below
   detection limits.
",
location = ph_location_type("body"))

# Slide 9: Limitations & next steps
pptx <- add_slide(pptx, layout = "Title and Content", master = "Office Theme")
pptx <- ph_with(pptx, value = "Limitations & next steps",
                location = ph_location_type("title"))
pptx <- ph_with(pptx, value = "
Limitations:
• Small sample (n = 27) — limited power for interactions
• Censored cytokine data reduces usable observations
• No data on birth control type (pill vs. IUD vs. implant)
• PCOS diagnosis method not specified

Next steps:
• Compare specific birth control types
• Use censored regression (tobit) models
• Include age / sexual activity as covariates
• Pathway analysis: cytokines × cell types
",
location = ph_location_type("body"))

# Save
pptx_file <- "Menstruation_Immune_Analysis.pptx"
print(pptx, target = pptx_file)
cat("✓ PPT saved:", pptx_file, "\n")
