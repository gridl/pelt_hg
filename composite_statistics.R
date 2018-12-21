#############################################################################
# Statistics for Pelt 4
# Written in R Version 3.5.0
#############################################################################
# Load data
compiled = read.csv("pelts_compiled.csv")

# Load Libraries
library(ggplot2)
library(boot)
library(lmtest)
library(car)
library(multcomp)

#########################################################
#### Pelt 1 Statistical Analysis ####
#### T-test ####
# Check assumptions of raw variables
shapiro.test(compiled$THg_TC) #Not normal 
hist(compiled$THg_TC, freq = FALSE, breaks = 30, xlab = "[THg] (ppm)") # Approximately normal

shapiro.test(compiled$THg_UC) # Not normal
hist(compiled$THg_UC, freq = FALSE, breaks = 30, xlab = "[THg] (ppm)") # Approximately normal

# Equal variances Levene's test
var.test(compiled$THg_TC, compiled$THg_UC) #variances are not equal

# Differene between TC and UC (paired)- Welch's paired t-test
t.test(compiled$THg_TC, compiled$THg_UC, paired = TRUE, var.equal = FALSE, 
       alternative = "two.sided")
shapiro.test((compiled$THg_TC-compiled$THg_UC)) # Differences Not normal
hist(compiled$THg_TC-compiled$THg_UC) # Differences Approximately normal

#### Variance between TC and UC ####
# Calculate absolute ranges
comp_TC_range = max(compiled$THg_TC) - min(compiled$THg_TC)
comp_UC_range = max(compiled$THg_UC) - min(compiled$THg_UC)

# Calculate ratio of ranges relative to TC THg range
UC_TC_ratio = (comp_UC_range)/(comp_TC_range)

# Calculate percent diff of ranges relative to TC THg range
percent_UC_TC_diff = (((comp_UC_range)/(comp_TC_range))-1)*100

# Calculate mean, SD and 95% CI (t(3, 0.05) = 3.1824) for percent diffs
paste(mean(percent_UC_TC_diff), "+/-",(3.1824*(260/sqrt(4))))

#### Anatomical Region Differences ####
# One-way ANOVA for Anatomical region 
# TOP COAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
comp_TC_anova1 = lm(compiled$THg_TC ~ compiled$Anatomical_region)
summary(comp_TC_anova1)
plot(comp_TC_anova1, main = "TC [THg] by anatomical region")

# Normality test of residuals
shapiro.test(residuals(comp_TC_anova1)) #not normal
hist(residuals(comp_TC_anova1), breaks = 30, freq = FALSE, 
     main = "Histogram of PELT 1 TC [THg] one-way ANOVA residuals by anatomical region")
curve(dnorm(x, mean(residuals(comp_TC_anova1)), sd(residuals(comp_TC_anova1))), 
      add=TRUE, col="darkblue", lwd=2)

# Homoscedasticity test of residuals
bptest(compiled$THg_TC ~ compiled$Anatomical_region) #homoscedastic

# Post-hoc Tukey's multiple comparison
# ANOVA 
comp_TC_posthoc = aov(compiled$THg_TC ~ compiled$Anatomical_region)
comp_TC_tukey = TukeyHSD(comp_TC_posthoc)
plot(comp_TC_tukey)

# UNDERCOAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
comp_UC_anova1 = lm(compiled$THg_UC ~ compiled$Anatomical_region)
summary(comp_UC_anova1)
plot(comp_UC_anova1, main = "UC [THg] by anatomical region")

# Normality test of residuals
shapiro.test(residuals(comp_UC_anova1)) #not normal
par(mfrow = c(1,1), las = 1, mar = c(5,5,5,5))
hist(residuals(comp_UC_anova1), breaks = 30, freq = FALSE)
curve(dnorm(x, mean(residuals(comp_UC_anova1)), sd(residuals(comp_UC_anova1))), 
      add=TRUE, col="darkblue", lwd=2) # Approximately normal

# Homoscedasticity test of residuals
leveneTest(compiled$THg_UC ~ compiled$Anatomical_region) # homoscedastic

#Tukey's multiple comparison
comp_UC_posthoc = aov(compiled$THg_UC ~ compiled$Anatomical_region)
comp_UC_tukey = TukeyHSD(comp_UC_posthoc)
plot(comp_UC_tukey)

#### Fur Region Differences ####
#TOP COAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
comp_TC_anova2 = lm(compiled$THg_TC ~ compiled$Fur_region)
comp_TC_anova2
plot(comp_TC_anova2, main = "TC [THg] by fur region")

#normality test of residuals
shapiro.test(residuals(comp_TC_anova2)) # Not normal

par(mfrow = c(1,1), las = 1, mar = c(5,5,5,5))
hist(residuals(comp_TC_anova2), breaks = 30, freq = FALSE, 
     main = "Histogram of PELT 1 TC [THg] one-way ANOVA residuals by fur region")
curve(dnorm(x, mean(residuals(comp_TC_anova2)), sd(residuals(comp_TC_anova2))), 
      add=TRUE, col="darkblue", lwd=2) # Approx normal


# Homoscedasticity test of residuals
leveneTest(compiled$THg_TC ~ compiled$Fur_region) # Homoscedastic

#Tukey's multiple comparison
comp_TC_posthoc2 = aov(compiled$THg_TC ~ compiled$Fur_region)
comp_TC_tukey2 = TukeyHSD(comp_TC_posthoc2)
plot(comp_TC_tukey2)

# UNDERCOAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
comp_UC_anova2 = lm(compiled$THg_UC ~ compiled$Fur_region)
comp_UC_anova2
plot(comp_UC_anova2, main = "UC [THg] by fur region")

#normality test of residuals
shapiro.test(residuals(comp_UC_anova2)) # Not normal

par(mfrow = c(1,1), las = 1, mar = c(5,5,5,5))
hist(residuals(comp_UC_anova2), breaks = 30, freq = FALSE, 
     main = "Histogram of PELT 1 UC [THg] one-way ANOVA residuals by fur region")
curve(dnorm(x, mean(residuals(comp_UC_anova2)), sd(residuals(comp_UC_anova2))), 
      add=TRUE, col="darkblue", lwd=2) # Approximately normal

# Homoscedasticity test of residuals
leveneTest(compiled$THg_UC ~ compiled$Fur_region) # Homoscedastic

#Tukey's multiple comparison
comp_UC_posthoc2 = aov(compiled$THg_UC ~ compiled$Fur_region)
comp_UC_tukey2 = TukeyHSD(comp_UC_posthoc2)
plot(comp_UC_tukey2)