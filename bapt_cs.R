# ==============================================================================
# Title: Descriptive Statistics and Linear Mixed Models for BAPT Analysis
# Project: Reliability of the timed BAPT and comparison in male football players
# Author: Marcos Quintana
# ==============================================================================
pacman::p_load(readxl, tidyverse, lmtest, lme4, emmeans, performance)

round_df <- function(d, digits = 3) {
  d %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = digits)))
}

df_base <- read_excel("basic_ds.xlsx") %>%
  mutate(
    id = factor(id)
  )

# Inspect columns and raw row-level counts
table(df_base$prior_groin)

# True player-level count
df_base %>%
  group_by(id) %>%
  summarise(
    has_prior_groin = any(prior_groin == "yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  count(has_prior_groin, name = "n_players")


# ==============================================================================
# 1. SUBJECT-LEVEL DEMOGRAPHICS (INJURED VS. UNINJURED)
# ==============================================================================
cat("\n========================================================\n")
cat("=== 1. SUBJECT-LEVEL DEMOGRAPHICS                    ===\n")
cat("========================================================\n")

# 1 row per unique subject (id = 236: 197 Uninjured, 39 Injured)
df_subj <- df_base %>%
  group_by(id) %>%
  summarise(
    group = factor(
      if_else(any(prior_groin == "yes", na.rm = TRUE), "Injured", "Uninjured"),
      levels = c("Uninjured", "Injured")
    ),
    age            = first(na.omit(age)),
    height         = first(na.omit(height)),
    body_mass      = first(na.omit(body_mass)),
    pitch_position = first(na.omit(pitch_position)),
    .groups = "drop"
  )

# Age, height, and body mass by group (Mean ± SD only, subject level)
cat("\n--- Age, Height, and Body Mass by Group (Subject Level) ---\n")
subj_cont <- df_subj %>%
  pivot_longer(cols = c(age, height, body_mass), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(group, variable) %>%
  summarise(
    n           = n(),
    `Mean ± SD` = paste0(round(mean(value), 1), " ± ", round(sd(value), 1)),
    .groups     = "drop"
  )
print(subj_cont)

# Pitch position by group (subject level)
cat("\n--- Pitch Position by Group (Subject Level) ---\n")
subj_pos <- df_subj %>%
  filter(!is.na(pitch_position)) %>%
  count(group, pitch_position, name = "n") %>%
  group_by(group) %>%
  mutate(
    pct     = round(100 * n / sum(n), 1)
  ) %>%
  ungroup()
print(subj_pos)


# Subject-level counts for knee arthroscopy
kas_by_subject <- df_base %>%
  group_by(id) %>%
  summarise(
    n_rows  = n(),
    n_kas   = sum(kas == "yes", na.rm = TRUE),
    has_kas = n_kas > 0,
    .groups = "drop"
  )

sum(kas_by_subject$has_kas)

count(kas_by_subject, n_kas, name = "n_subjects")

df_base %>%
  group_by(id) %>%
  summarise(
    has_groin = any(prior_groin == "yes"),
    has_kas   = any(kas == "yes"),
    .groups   = "drop"
  ) %>%
  count(has_groin, has_kas)

# Dominant foot
subj_dom <- df_base %>%
  filter(dominant == "yes") %>%
  left_join(df_subj %>% select(id, group), by = "id") %>%
  count(group, limb, name = "n") %>%
  group_by(group) %>%
  mutate(
    pct     = round(100 * n / sum(n), 1)
  ) %>%
  ungroup()
print(subj_dom)


# ==============================================================================
# 2. ANALYSIS 1: TIME UNDER TENSION (BAPT_TIME)
# ==============================================================================
cat("\n========================================================\n")
cat("=== ANALYSIS 1: BAPT_TIME (DURATION)                 ===\n")
cat("========================================================\n")

df_time <- df_base %>% select(id, prior_groin, bapt_time) %>% drop_na()

cat("Total observations:", nrow(df_time), "| Unique subjects:", n_distinct(df_time$id), "\n\n")

# Manuscript descriptive table (BAPT_TIME)
tabla_paper_time <- df_time %>%
  group_by(prior_groin) %>%
  summarise(
    n               = n(),
    `Mean ± SD`     = paste0(round(mean(bapt_time), 1), " ± ", round(sd(bapt_time), 1), " s"),
    `Median [IQR]`  = paste0(round(median(bapt_time), 1), " [", 
                            round(as.numeric(quantile(bapt_time, 0.25)), 1), " to ", 
                            round(as.numeric(quantile(bapt_time, 0.75)), 1), "] s")
  )
print(tabla_paper_time)

# Linear Mixed-Effects Model (Estimated Marginal Means)
fit_lmm_time <- lmer(bapt_time ~ prior_groin + (1 | id), data = df_time)
cat("\n--- Estimated Marginal Means (BAPT_TIME) ---\n")

emm_time <- emmeans(fit_lmm_time, ~ prior_groin)
print(emm_time)
print(confint(emm_time)) # Confidence intervals for marginal means

cat("\n--- Pairwise Comparisons (Contrasts) ---\n")
contrast_pairs <- pairs(emm_time)
print(contrast_pairs)

cat("\n--- 95% Confidence Intervals for Contrasts ---\n")
print(confint(contrast_pairs))

# Effect Size (Cohen's d) based on the linear mixed model
cat("\n--- Effect Size (Cohen's d) ---\n")
eff_time <- eff_size(emm_time, sigma = sigma(fit_lmm_time), edf = df.residual(fit_lmm_time))
print(eff_time)

# Model assumptions check
check_model(fit_lmm_time, check = c("linearity", "qq", "normality"))


# ==============================================================================
# 3. ANALYSIS 2: REPETITIONS (BAPT_REPS)
# ==============================================================================
cat("\n========================================================\n")
cat("=== ANALYSIS 2: BAPT_REPS (REPETITIONS)              ===\n")
cat("========================================================\n")

df_reps <- df_base %>% select(id, prior_groin, bapt_reps) %>% drop_na()

cat("Total observations:", nrow(df_reps), "| Unique subjects:", n_distinct(df_reps$id), "\n\n")

# Manuscript descriptive table (BAPT_REPS)
tabla_paper_reps <- df_reps %>%
  group_by(prior_groin) %>%
  summarise(
    n               = n(),
    `Mean ± SD`     = paste0(round(mean(bapt_reps), 1), " ± ", round(sd(bapt_reps), 1), " reps"),
    `Median [IQR]`  = paste0(round(median(bapt_reps), 1), " [", 
                            round(as.numeric(quantile(bapt_reps, 0.25)), 1), " to ", 
                            round(as.numeric(quantile(bapt_reps, 0.75)), 1), "] reps")
  )
print(tabla_paper_reps)

# Linear Mixed-Effects Model for bapt_reps
fit_lmm_reps <- lmer(bapt_reps ~ prior_groin + (1 | id), data = df_reps)
cat("\n--- Estimated Marginal Means (BAPT_REPS) ---\n")

emm_reps <- emmeans(fit_lmm_reps, ~ prior_groin)
print(emm_reps)
print(confint(emm_reps)) # Confidence intervals for marginal means

cat("\n--- Pairwise Comparisons (Contrasts) ---\n")
contrast_pairs_reps <- pairs(emm_reps)
print(contrast_pairs_reps)

cat("\n--- 95% Confidence Intervals for Contrasts ---\n")
print(confint(contrast_pairs_reps))

# Effect Size (Cohen's d) based on the linear mixed model
cat("\n--- Effect Size (Cohen's d) ---\n")
eff_reps <- eff_size(emm_reps, sigma = sigma(fit_lmm_reps), edf = df.residual(fit_lmm_reps))
print(eff_reps)

# Model assumptions check
check_model(fit_lmm_reps, check = c("linearity", "qq", "normality"))