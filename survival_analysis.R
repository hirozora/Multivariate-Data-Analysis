# Survival Analysis for Cancer Dataset (Simplified Version)

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install and load survival package
if (!require("survival", quietly = TRUE)) {
  cat("Installing survival package...\n")
  install.packages("survival", dependencies = TRUE)
}
library(survival)

# Load data
cancer_data <- read.csv("cancer.csv", row.names = 1)

# Data preprocessing - Recode status variable
# Original: 1 = censored, 2 = dead
# Recoded: 0 = censored, 1 = event (dead)
cancer_data$status_recoded <- ifelse(cancer_data$status == 2, 1, 0)

# Display data summary
cat("=== Data Summary ===\n")
summary(cancer_data[, c("time", "status_recoded", "age", "sex", "ph.ecog")])

# 1. Survival Curve, Median Survival Time, and Survival Rate

surv_obj <- Surv(time = cancer_data$time, event = cancer_data$status_recoded)
km_fit <- survfit(surv_obj ~ 1, data = cancer_data)
cat("\n=== Survival Analysis Summary ===\n")
print(summary(km_fit))
cat("\n=== Median Survival Time ===\n")
print(km_fit)
cat("\n=== Survival Rate at 365 Days (1 Year) ===\n")
surv_365 <- summary(km_fit, times = 365)
cat(sprintf("Time: %d days\n", surv_365$time))
cat(sprintf("Survival Rate: %.4f\n", surv_365$surv))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", surv_365$lower, surv_365$upper))

# Plot survival curve using base graphics
png("survival_curve.png", width = 800, height = 600)
plot(km_fit, 
     xlab = "Time (days)", 
     ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curve",
     conf.int = TRUE,
     col = "blue",
     lwd = 2)
grid()
legend("topright", 
       legend = c("Survival Curve", "95% CI"), 
       col = c("blue", "lightblue"), 
       lwd = c(2, 1),
       bty = "n")
dev.off()
cat("\nSurvival curve saved as 'survival_curve.png'\n")

# 2. Cox Proportional Hazards Model
# Remove rows with missing values in covariates
cancer_clean <- na.omit(cancer_data[, c("time", "status_recoded", "sex", "age", "ph.ecog")])

# Fit Cox model
cox_model <- coxph(Surv(time, status_recoded) ~ sex + age + ph.ecog, 
                   data = cancer_clean)
cat("\n=== Cox Proportional Hazards Model ===\n")
print(summary(cox_model))
cat("\n=== Cox Model Coefficients ===\n")
coef_table <- summary(cox_model)$coefficients
print(coef_table)

cat("\n=== Hazard Ratios (with 95% CI) ===\n")
hr_ci <- summary(cox_model)$conf.int
print(hr_ci)

# 3. Test Proportional Hazards Assumption

# Schoenfeld residuals test
cat("\n=== Proportional Hazards Assumption Test (Schoenfeld Residuals) ===\n")
ph_test <- cox.zph(cox_model)
print(ph_test)

cat("\n=== Interpretation ===\n")
cat("If p-value > 0.05, the proportional hazards assumption holds.\n")
cat("If p-value <= 0.05, the assumption is violated.\n\n")

# Check each covariate
for (i in 1:nrow(ph_test$table)) {
  var_name <- rownames(ph_test$table)[i]
  p_val <- ph_test$table[i, "p"]
  if (p_val > 0.05) {
    cat(sprintf("%s: p = %.4f (assumption holds)\n", var_name, p_val))
  } else {
    cat(sprintf("%s: p = %.4f (assumption VIOLATED)\n", var_name, p_val))
  }
}

# Plot Schoenfeld residuals
png("schoenfeld_residuals.png", width = 1200, height = 400)
par(mfrow = c(1, 3))
plot(ph_test)
dev.off()
cat("\nSchoenfeld residuals plot saved as 'schoenfeld_residuals.png'\n")

# Summary of Results
cat("\n=== Analysis Complete ===\n")
cat("Files generated:\n")
cat("1. survival_curve.png - Kaplan-Meier survival curve\n")
cat("2. schoenfeld_residuals.png - Proportional hazards test plot\n")

