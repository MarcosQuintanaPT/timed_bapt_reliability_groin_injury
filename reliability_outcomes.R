# ==============================================================================
# Title: Test-Retest Reliability and Bland-Altman Analysis for BAPT
# Project: Reliability of the timed BAPT and comparison in male football players
# Author: Marcos Quintana
# ==============================================================================

# Load required packages
pacman::p_load(tidyverse, psych, readxl, ICC.Sample.Size, rmcorr, BlandAltmanLeh, skimr, cowplot)

# ==============================================================================
# 1. SAMPLE SIZE & DATA LOADING
# ==============================================================================

# Sample size calculation based on ICC expected values
calc_size <- calculateIccSampleSize(p = 0.9, p0 = 0.5, k = 2, 
                                   alpha = 0.05, tails = 2, 
                                   power = 0.8)
print(calc_size)

# Load the reliability dataset from the data directory
df <- read_excel("reliability_bapt.xlsx") %>% 
  na.omit()

# Inspect descriptive statistics
df |> 
  skim(r1, r2, tut1, tut2, age, height, body_mass)


# ==============================================================================
# 2. RELIABILITY METRICS FUNCTION
# ==============================================================================
calculate_metrics <- function(data, col1, col2) {
  res_icc <- ICC(data[, c(col1, col2)], lmer = FALSE)
  
  val_icc <- res_icc$results$ICC[3] 
  low_icc <- res_icc$results$`lower bound`[3]
  up_icc  <- res_icc$results$`upper bound`[3]
  
  sd_p    <- sqrt((sd(data[[col1]], na.rm = TRUE)^2 + sd(data[[col2]], na.rm = TRUE)^2) / 2)
  val_sem <- sd_p * sqrt(1 - val_icc)
  val_mdc <- val_sem * 1.96 * sqrt(2)
  cv_v    <- mean((abs(data[[col1]] - data[[col2]]) / ((data[[col1]] + data[[col2]])/2)) * 100)
  
  return(list(icc = val_icc, icc_l = low_icc, icc_u = up_icc, sem = val_sem, mdc = val_mdc, cv = cv_v))
}


# ==============================================================================
# 3. RESTRUCTURING FOR REPEATED MEASURES CORRELATION (RMCORR)
# ==============================================================================
df_long <- df %>%
  pivot_longer(
    cols = -sujeto, 
    names_to = c(".value", "session"), 
    names_pattern = "([a-z]+)(\\d)"
  )


# ==============================================================================
# 4. REPEATED MEASURES CORRELATION
# ==============================================================================
rmcorr(participant = sujeto, measure1 = tut, measure2 = r, dataset = df_long)


# ==============================================================================
# 5. FINAL RELIABILITY SUMMARY TABLE (ICC, SEM, MDC, CV)
# ==============================================================================
res_tut <- calculate_metrics(df, "tut1", "tut2")
res_r   <- calculate_metrics(df, "r1", "r2")

cat("\nRELIABILITY RESULTS SUMMARY\n")
resumen_final <- data.frame(
  Variable  = c("Time under tension (s)", "Repetitions (n)"),
  ICC       = c(round(res_tut$icc, 3), round(res_r$icc, 3)),
  CI_Lower  = c(round(res_tut$icc_l, 3), round(res_r$icc_l, 3)),
  CI_Upper  = c(round(res_tut$icc_u, 3), round(res_r$icc_u, 3)),
  SEM       = c(round(res_tut$sem, 2),   round(res_r$sem, 2)),
  MDC95     = c(round(res_tut$mdc, 2),   round(res_r$mdc, 2)),
  CV_perc   = c(round(res_tut$cv, 2),    round(res_r$cv, 2))
)

print(resumen_final, row.names = FALSE)


# ==============================================================================
# 6. BLAND-ALTMAN ANALYSIS & PLOTS
# ==============================================================================

# Helper function to compute Bland-Altman differences and means
bland_altman_data <- function(x, y) {
  mean_val <- (x + y) / 2
  diff_val <- x - y
  data.frame(mean = mean_val, diff = diff_val)
}

# --- A. Time Under Tension (tut1 vs tut2) ---
t1 <- as.numeric(df$tut1)
t2 <- as.numeric(df$tut2)
ba_tut_df <- bland_altman_data(t1, t2)

bias_tut  <- mean(ba_tut_df$diff, na.rm = TRUE)
upper_tut <- bias_tut + 1.96 * sd(ba_tut_df$diff, na.rm = TRUE)
lower_tut <- bias_tut - 1.96 * sd(ba_tut_df$diff, na.rm = TRUE)

cat("\n=== BLAND-ALTMAN: Time Under Tension (Session 1 vs Session 2) ===\n")
cat("Mean Bias: ", round(bias_tut, 2), "\n")
cat("95% Limits of Agreement: [", round(lower_tut, 2), ", ", round(upper_tut, 2), "]\n", sep = "")

p_ba_tut <- ggplot(ba_tut_df, aes(x = mean, y = diff)) +
  geom_point(color = "skyblue3", size = 2.2, alpha = 0.7) +
  geom_hline(yintercept = bias_tut, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = upper_tut, color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = lower_tut, color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  theme_classic() +
  labs(
    x = "Mean Time Under Tension (Session 1 & 2) [s]",
    y = "Difference (Session 1 - Session 2) [s]"
  )

print(p_ba_tut)

# --- B. Repetitions (r1 vs r2) ---
r1_val <- as.numeric(df$r1)
r2_val <- as.numeric(df$r2)
ba_r_df <- bland_altman_data(r1_val, r2_val)

bias_r  <- mean(ba_r_df$diff, na.rm = TRUE)
upper_r <- bias_r + 1.96 * sd(ba_r_df$diff, na.rm = TRUE)
lower_r <- bias_r - 1.96 * sd(ba_r_df$diff, na.rm = TRUE)

cat("\n=== BLAND-ALTMAN: Repetitions (Session 1 vs Session 2) ===\n")
cat("Mean Bias: ", round(bias_r, 2), "\n")
cat("95% Limits of Agreement: [", round(lower_r, 2), ", ", round(upper_r, 2), "]\n", sep = "")

p_ba_r <- ggplot(ba_r_df, aes(x = mean, y = diff)) +
  geom_point(color = "skyblue3", size = 2.2, alpha = 0.7) +
  geom_hline(yintercept = bias_r, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = upper_r, color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = lower_r, color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  theme_classic() +
  labs(
    x = "Mean Repetitions (Session 1 & 2) [n]",
    y = "Difference (Session 1 - Session 2) [n]"
  )

print(p_ba_r)

# Merge Bland-Altman plots into a side-by-side layout
baltman_p <- plot_grid(p_ba_tut, p_ba_r, nrow = 1, align = "h", labels = "auto")
print(baltman_p)