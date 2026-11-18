# ==============================================================================
# Title: Descriptive Statistics and Linear Mixed Models for BAPT Analysis
# Project: Reliability of the timed BAPT and comparison in male football players
# Author: Marcos Quintana
# ==============================================================================

# Load required packages
pacman::p_load(readxl, tidyverse, lmtest, lme4, emmeans, performance)

# Utility function to round numeric columns in a data frame
round_df <- function(d, digits = 3) {
  nums <- vapply(d, is.numeric, FUN.VALUE = logical(1))
  d[, nums] <- round(d[, nums], digits = digits)
  return(d)
}

# Load the cleaned dataset from the data directory
df_base <- read_excel("Cross_sectional_cleaned_base.xlsx")

# Ensure 'id' and 'grupo' (group) are treated as factors
df_base <- df_base %>%
  mutate(
    id = factor(id),
    grupo = factor(grupo)
  )


# ==============================================================================
# 1. SUBJECT-LEVEL DESCRIPTIVE STATISTICS (Table 1 / Demographics)
# ==============================================================================
cat("\n========================================================\n")
cat("=== SUBJECT-LEVEL DESCRIPTIONS (PATIENT LEVEL)        ===\n")
cat("========================================================\n")

# Filter to one row per unique subject for demographic and baseline variables
df_subjects <- df_base %>%
  distinct(id, .keep_all = TRUE)

cat("Total unique subjects:", nrow(df_subjects), "\n\n")

# A. Continuous variables (Age, Height, Body Mass)
tabla_continuous_subjects <- df_subjects %>%
  group_by(grupo) %>%
  summarise(
    n                    = n(),
    `Age (Mean ± SD)`    = paste0(round(mean(age, na.rm = TRUE), 1), " ± ", round(sd(age, na.rm = TRUE), 1)),
    `Height (Mean ± SD)` = paste0(round(mean(height, na.rm = TRUE), 1), " ± ", round(sd(height, na.rm = TRUE), 1)),
    `Mass (Mean ± SD)`   = paste0(round(mean(body_mass, na.rm = TRUE), 1), " ± ", round(sd(body_mass, na.rm = TRUE), 1))
  )
print(tabla_continuous_subjects)

# B. Categorical variables (Pitch Position, KAS)
cat("\n--- Categorical Variables Distribution (Subjects) ---\n")
print(table(df_subjects$pitch_position, df_subjects$grupo, useNA = "ifany"))
print(table(df_subjects$kas, df_subjects$grupo, useNA = "ifany"))


# ==============================================================================
# 2. ANALYSIS 1: TIME UNDER TENSION (BAPT_TIME)
# ==============================================================================
cat("\n========================================================\n")
cat("=== ANALYSIS 1: BAPT_TIME (DURATION)                 ===\n")
cat("========================================================\n")

df_time <- df_base %>% select(id, prior_groin_num, grupo, bapt_time) %>% drop_na()

cat("Total observations:", nrow(df_time), "| Unique subjects:", n_distinct(df_time$id), "\n\n")

# Manuscript descriptive table (BAPT_TIME)
tabla_paper_time <- df_time %>%
  group_by(grupo) %>%
  summarise(
    n               = n(),
    `Mean ± SD`     = paste0(round(mean(bapt_time), 1), " ± ", round(sd(bapt_time), 1), " s"),
    `Median [IQR]`  = paste0(round(median(bapt_time), 1), " [", 
                            round(as.numeric(quantile(bapt_time, 0.25)), 1), " to ", 
                            round(as.numeric(quantile(bapt_time, 0.75)), 1), "] s")
  )
print(tabla_paper_time)

# Linear Mixed-Effects Model (Estimated Marginal Means)
fit_lmm_time <- lmer(bapt_time ~ grupo + (1 | id), data = df_time)
cat("\n--- Estimated Marginal Means (BAPT_TIME) ---\n")

emm_time <- emmeans(fit_lmm_time, ~ grupo)
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

df_reps <- df_base %>% select(id, prior_groin_num, grupo, bapt_reps) %>% drop_na()

cat("Total observations:", nrow(df_reps), "| Unique subjects:", n_distinct(df_reps$id), "\n\n")

# Manuscript descriptive table (BAPT_REPS)
tabla_paper_reps <- df_reps %>%
  group_by(grupo) %>%
  summarise(
    n               = n(),
    `Mean ± SD`     = paste0(round(mean(bapt_reps), 1), " ± ", round(sd(bapt_reps), 1), " reps"),
    `Median [IQR]`  = paste0(round(median(bapt_reps), 1), " [", 
                            round(as.numeric(quantile(bapt_reps, 0.25)), 1), " to ", 
                            round(as.numeric(quantile(bapt_reps, 0.75)), 1), "] reps")
  )
print(tabla_paper_reps)

# Linear Mixed-Effects Model for bapt_reps
fit_lmm_reps <- lmer(bapt_reps ~ grupo + (1 | id), data = df_reps)
cat("\n--- Estimated Marginal Means (BAPT_REPS) ---\n")

emm_reps <- emmeans(fit_lmm_reps, ~ grupo)
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