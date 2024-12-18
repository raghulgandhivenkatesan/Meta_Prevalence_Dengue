####################Load necessary libraries#######################
library(metafor)
library(meta)

####################Data for meta-analysis########################
dengue_data <- data.frame(
  Author = c("Buntubatu et al., (2019)", "Baqi et al., (2022)", "Bhatt et al., (2020)", 
             "Gupta et al., (2022)", "Hussain et al., (2017)", "Li et al., (2016)", 
             "Mansanguan et al., (2021)", "Salgado et al., (2010)", "Satarasinghe et al., (2007)", 
             "Soneja et al., (2019)", "Weerakoon et al., (2011)", "Yadav et al., (2013)"),
  Events = c(39, 42, 13, 40, 24, 201, 2, 11, 42, 11, 45, 32), # Myocarditis cases
  N = c(50, 1008, 182, 150, 128, 351, 81, 102, 174, 183, 319, 67) # Total sample size
)

####################Step 1: Calculate observed proportions######################
ies <- escalc(measure = "PR", xi = Events, ni = N, data = dengue_data)
print(ies)

####################Step 2: Perform random-effects meta-analysis##################
pes <- rma(yi, vi, data = ies, method = "DL", weighted = TRUE)
print(pes, digits = 2)
confint(pes)

####################Step 3: Identifying outliers with residuals######################
stud_res <- rstudent(pes)
abs_z <- abs(stud_res$z)
outliers <- stud_res[order(-abs_z), ]
print(outliers)

####################Step 4: Leave-One-Out Analysis####################################
l1o <- leave1out(pes)
# Extract estimates and confidence intervals
yi <- l1o$estimate             # Leave-one-out effect sizes
sei <- l1o$se                  # Standard errors
ci.lb <- yi - 1.96 * sei       # Lower bounds of 95% CI
ci.ub <- yi + 1.96 * sei       # Upper bounds of 95% CI
# Create Forest Plot for Leave-One-Out Analysis
forest(yi, sei = sei,
       slab = dengue_data$Author, 
       xlab = "Summary Proportions Leaving Out Each Study",
       refline = pes$b,         # Overall summary proportion
       digits = 4,              # Number of digits
       alim = c(min(ci.lb), max(ci.ub)),  # Adjust x-axis limits
       cex = 0.8,               # Text size
       col = "blue",            # Point color
       main = "Leave-One-Out Sensitivity Analysis")

####################Step 5: Baujat Plot (Heterogeneity Diagnostics)####################
baujat(pes, main = "Baujat Plot for Heterogeneity Diagnostics")

####################Step 6: Influence Diagnostics######################################
inf <- influence(pes)
print(inf)
# Plot the Influence Diagnostics
plot(inf)

####################Step 7: Remove Outliers (if needed, for example 1 and 6)############
outlier_removed <- ies[-c(1, 6), ]
pes_updated <- rma(yi, vi, data = outlier_removed, method = "DL", weighted = TRUE)
print(pes_updated, digits = 2)

####################Step 8: Forest Plot for Overall Meta-Analysis#####################
meta_summary <- metaprop(
  event = Events, 
  n = N, 
  studlab = Author, 
  data = dengue_data, 
  sm = "PRAW", method.ci = "NAsm", method.tau = "DL", 
  incr = 0.5, tau.common = TRUE, title = "Forest Plot for Dengue Data"
)

# Save forest plot as an image
png("forestplot_dengue.png", width = 1000, height = 1000)
forest(meta_summary,
       xlim = c(0, 1), pscale = 1, 
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95% C.I.", "Weights"),
       leftcols = c("studlab", "event", "n"), 
       leftlabs = c("Study", "Cases", "Total"),
       squaresize = 0.5, col.square = "navy",
       col.diamond = "red", col.diamond.lines = "navy",
       type.random = "diamond", print.Q = TRUE, print.I2 = TRUE)
dev.off()

#####################Step 9: Funnel Plot for Publication Bias###########################
funnel(pes, main = "Funnel Plot for Dengue Data")

#####################Step 10: Trim-and-Fill Analysis for Publication Bias###################
pes_trimfill <- trimfill(pes)
print(pes_trimfill)
funnel(pes_trimfill, main = "Trim-and-Fill Funnel Plot")

#####################Step 11: Egger's Regression Test for Small-Study Effects##############
egger_test <- regtest(pes, model = "lm", predictor = "sei")
print(egger_test)

#####################Step 12: Rank Correlation Test######################################
rank_corr <- ranktest(pes)
print(rank_corr)


####################Step 13: Perform meta-regression with sample size as a moderator#################
meta_reg_sample <- rma(yi, vi, mods = ~ N, data = ies, method = "DL")
print(meta_reg_sample, digits = 2)
# Summary of meta-regression
summary(meta_reg_sample)

####################Step 14: Bubble Plot for Meta-Regression with Sample Size#############################
bubble <- predict(meta_reg_sample, newmods = ies$N)
plot(ies$N, ies$yi,
     xlab = "Sample Size (N)", 
     ylab = "Observed Effect Size (Log Proportions)", 
     main = "Bubble Plot for Meta-Regression (Sample Size)",
     pch = 21, bg = "lightblue", cex = 1.2)
# Add meta-regression line
lines(ies$N, bubble$pred, col = "red", lwd = 2)
# Add legend
legend("topright", legend = "Meta-Regression Line", col = "red", lwd = 2)


######################Step 15: Generate Distribution of the True Effect###########################
# Create a sequence of possible true effects based on the random-effects model
true_effects <- seq(pes$b - 4 * sqrt(pes$tau2), pes$b + 4 * sqrt(pes$tau2), length.out = 1000)
# Density of the distribution
density_values <- dnorm(true_effects, mean = pes$b, sd = sqrt(pes$tau2))
# Plot the distribution of the true effect
plot(true_effects, density_values, type = "l", lwd = 2, col = "blue",
     xlab = "True Effect Size (Log Proportions)",
     ylab = "Density",
     main = "Distribution of the True Effect")
# Add vertical line for the estimated overall effect size
abline(v = pes$b, col = "red", lwd = 2, lty = 2)
# Add legend
legend("topright", legend = c("True Effect Distribution", "Overall Effect Size"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)



