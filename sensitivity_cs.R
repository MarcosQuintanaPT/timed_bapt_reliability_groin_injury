# ==============================================================================
# Title: Supplementary Sensitivity Analysis (4-Level Limb Injury Status)
# Project: Reliability of the timed BAPT and comparison in male football players
# Author: Marcos Quintana
# ==============================================================================
pacman::p_load(readxl, tidyverse, lme4, lmerTest, emmeans, performance)

round_df <- function(d, digits = 3) {
  d %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = digits)))
}

# ==============================================================================
# 0. DATA IMPORT & 4-LEVEL LIMB STATUS CODING
# ==============================================================================
df_base <- read_excel("basic_ds.xlsx") %>%
  mutate(
    id          = factor(id),
    limb        = factor(limb),
    prior_groin = factor(prior_groin, levels = c("no", "yes"))
  ) %>%
  group_by(id) %>%
  mutate(
    n_groin_limbs = sum(prior_groin == "yes", na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    limb_status = case_when(
      prior_groin == "no"  & n_groin_limbs == 0 ~ "control",
      prior_groin == "no"  & n_groin_limbs == 1 ~ "contralateral",
      prior_groin == "yes" & n_groin_limbs == 1 ~ "injured_uni",
      prior_groin == "yes" & n_groin_limbs >= 2 ~ "injured_bi",
      TRUE ~ NA_character_
    ),
    limb_status = factor(
      limb_status,
      levels = c("control", "contralateral", "injured_uni", "injured_bi")
    )
  )

cat("\n--- Sensitivity Limb-Status Counts (Limbs & Unique Players) ---\n")
df_base %>%
  filter(!is.na(limb_status)) %>%
  group_by(limb_status) %>%
  summarise(
    n_limbs   = n(),
    n_players = n_distinct(id),
    .groups   = "drop"
  ) %>%
  print()

cat("\n--- Unilateral Players Not Contributing 2 Tested Limbs ---\n")
df_base %>%
  filter(id %in% id[limb_status %in% c("injured_uni", "contralateral")]) %>%
  count(id, name = "n_tested_limbs") %>%
  filter(n_tested_limbs != 2) %>%
  print()


# ==============================================================================
# 1. SENSITIVITY ANALYSIS: TIME UNDER TENSION (BAPT_TIME)
# ==============================================================================
cat("\n========================================================\n")
cat("=== SENSITIVITY ANALYSIS: BAPT_TIME (LIMB_STATUS)    ===\n")
cat("========================================================\n")

df_time <- df_base %>%
  select(id, limb, prior_groin, limb_status, bapt_time) %>%
  drop_na()

cat("Total observations:", nrow(df_time), "| Unique subjects:", n_distinct(df_time$id), "\n\n")

tabla_sens_time <- df_time %>%
  group_by(limb_status) %>%
  summarise(
    n_limbs     = n(),
    n_players   = n_distinct(id),
    `Mean ± SD` = paste0(round(mean(bapt_time), 1), " ± ", round(sd(bapt_time), 1), " s"),
    .groups     = "drop"
  )
print(tabla_sens_time)

fit_sens_time <- lmer(bapt_time ~ limb_status + (1 | id), data = df_time)
cat("\n--- Model Summary (BAPT_TIME) ---\n")
print(summary(fit_sens_time))

# Optional likelihood-ratio comparison vs. primary binary limb model (ML)
fit_base_time_ml <- lmer(bapt_time ~ prior_groin + (1 | id), data = df_time, REML = FALSE)
fit_sens_time_ml <- update(fit_sens_time, REML = FALSE)
cat("\n--- LRT vs. Primary Binary Model (BAPT_TIME) ---\n")
print(anova(fit_base_time_ml, fit_sens_time_ml))

cat("\n--- Estimated Marginal Means & 95% CIs (BAPT_TIME) ---\n")
emm_sens_time <- emmeans(fit_sens_time, ~ limb_status)
print(confint(emm_sens_time))

cat("\n--- Tukey-Adjusted Pairwise Contrasts & 95% CIs (BAPT_TIME) ---\n")
pairs_sens_time <- pairs(emm_sens_time, adjust = "tukey")
print(pairs_sens_time)
print(confint(pairs_sens_time))

cat("\n--- Effect Sizes: Cohen's d (BAPT_TIME) ---\n")
eff_sens_time <- eff_size(emm_sens_time, sigma = sigma(fit_sens_time), edf = df.residual(fit_sens_time))
print(eff_sens_time)

check_model(fit_sens_time, check = c("linearity", "qq", "normality"))


# ==============================================================================
# 2. SENSITIVITY ANALYSIS: REPETITIONS (BAPT_REPS)
# ==============================================================================
cat("\n========================================================\n")
cat("=== SENSITIVITY ANALYSIS: BAPT_REPS (LIMB_STATUS)    ===\n")
cat("========================================================\n")

df_reps <- df_base %>%
  select(id, limb, prior_groin, limb_status, bapt_reps) %>%
  drop_na()

cat("Total observations:", nrow(df_reps), "| Unique subjects:", n_distinct(df_reps$id), "\n\n")

tabla_sens_reps <- df_reps %>%
  group_by(limb_status) %>%
  summarise(
    n_limbs     = n(),
    n_players   = n_distinct(id),
    `Mean ± SD` = paste0(round(mean(bapt_reps), 1), " ± ", round(sd(bapt_reps), 1), " reps"),
    .groups     = "drop"
  )
print(tabla_sens_reps)

fit_sens_reps <- lmer(bapt_reps ~ limb_status + (1 | id), data = df_reps)
cat("\n--- Model Summary (BAPT_REPS) ---\n")
print(summary(fit_sens_reps))

# Optional likelihood-ratio comparison vs. primary binary limb model (ML)
fit_base_reps_ml <- lmer(bapt_reps ~ prior_groin + (1 | id), data = df_reps, REML = FALSE)
fit_sens_reps_ml <- update(fit_sens_reps, REML = FALSE)
cat("\n--- LRT vs. Primary Binary Model (BAPT_REPS) ---\n")
print(anova(fit_base_reps_ml, fit_sens_reps_ml))

cat("\n--- Estimated Marginal Means & 95% CIs (BAPT_REPS) ---\n")
emm_sens_reps <- emmeans(fit_sens_reps, ~ limb_status)
print(confint(emm_sens_reps))

cat("\n--- Tukey-Adjusted Pairwise Contrasts & 95% CIs (BAPT_REPS) ---\n")
pairs_sens_reps <- pairs(emm_sens_reps, adjust = "tukey")
print(pairs_sens_reps)
print(confint(pairs_sens_reps))

cat("\n--- Effect Sizes: Cohen's d (BAPT_REPS) ---\n")
eff_sens_reps <- eff_size(emm_sens_reps, sigma = sigma(fit_sens_reps), edf = df.residual(fit_sens_reps))
print(eff_sens_reps)

check_model(fit_sens_reps, check = c("linearity", "qq", "normality"))