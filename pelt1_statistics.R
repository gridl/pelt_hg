#############################################################################
# Statistics for Pelt 1
# Written in R Version 3.5.0
#############################################################################
# Load data
pelt1 = read.csv("PELT1.csv")

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
shapiro.test(pelt1$THg_TC) #Not normal 
hist(pelt1$THg_TC, freq = FALSE, breaks = 30, xlab = "[THg] (ppm)") # Approximately normal

shapiro.test(pelt1$THg_UC) # Not normal
hist(pelt1$THg_UC, freq = FALSE, breaks = 30, xlab = "[THg] (ppm)") # Approximately normal

# Equal variances Levene's test
var.test(pelt1$THg_TC, pelt1$THg_UC) #variances are not equal

# Differene between TC and UC (paired)- Welch's paired t-test
t.test(pelt1$THg_TC, pelt1$THg_UC, paired = TRUE, var.equal = FALSE, 
       alternative = "two.sided")
shapiro.test((pelt1$THg_TC-pelt1$THg_UC)) # Differences Not normal
hist(pelt1$THg_TC-pelt1$THg_UC) # Differences Approximately normal

#### Variance between TC and UC ####
# Calculate absolute ranges
p1_TC_range = max(pelt1$THg_TC) - min(pelt1$THg_TC)
p1_UC_range = max(pelt1$THg_UC) - min(pelt1$THg_UC)

# Calculate ratio of ranges relative to TC THg range
UC_TC_ratio = (p1_UC_range)/(p1_TC_range)

# Calculate percent diff of ranges relative to TC THg range
percent_UC_TC_diff = (((p1_UC_range)/(p1_TC_range))-1)*100

# Calculate mean, SD and 95% CI (t(3, 0.05) = 3.1824) for percent diffs
paste(mean(percent_UC_TC_diff), "+/-",(3.1824*(260/sqrt(4))))

#### Anatomical Region Differences ####
# One-way ANOVA for Anatomical region 
# TOP COAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
p1_TC_anova1 = lm(pelt1$THg_TC ~ pelt1$Anatomical_region)
summary(p1_TC_anova1)
plot(p1_TC_anova1, main = "TC [THg] by anatomical region")

# Normality test of residuals
shapiro.test(residuals(p1_TC_anova1)) #not normal
hist(residuals(p1_TC_anova1), breaks = 30, freq = FALSE, 
     main = "Histogram of PELT 1 TC [THg] one-way ANOVA residuals by anatomical region")
curve(dnorm(x, mean(residuals(p1_TC_anova1)), sd(residuals(p1_TC_anova1))), 
      add=TRUE, col="darkblue", lwd=2)

# Homoscedasticity test of residuals
bptest(pelt1$THg_TC ~ pelt1$Anatomical_region) #homoscedastic

# Post-hoc Tukey's multiple comparison
# ANOVA 
p1_TC_posthoc = aov(pelt1$THg_TC ~ pelt1$Anatomical_region)
p1_TC_tukey = TukeyHSD(p1_TC_posthoc)
plot(p1_TC_tukey)

# UNDERCOAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
p1_UC_anova1 = lm(pelt1$THg_UC ~ pelt1$Anatomical_region)
summary(p1_UC_anova1)
plot(p1_UC_anova1, main = "UC [THg] by anatomical region")

# Normality test of residuals
shapiro.test(residuals(p1_UC_anova1)) #not normal
par(mfrow = c(1,1), las = 1, mar = c(5,5,5,5))
hist(residuals(p1_UC_anova1), breaks = 30, freq = FALSE)
curve(dnorm(x, mean(residuals(p1_UC_anova1)), sd(residuals(p1_UC_anova1))), 
      add=TRUE, col="darkblue", lwd=2) # Approximately normal

# Homoscedasticity test of residuals
leveneTest(pelt1$THg_UC ~ pelt1$Anatomical_region) # homoscedastic

#Tukey's multiple comparison
p1_UC_posthoc = aov(pelt1$THg_UC ~ pelt1$Anatomical_region)
p1_UC_tukey = TukeyHSD(p1_UC_posthoc)
plot(p1_UC_tukey)

#### Fur Region Differences ####
#TOP COAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
p1_TC_anova2 = lm(pelt1$THg_TC ~ pelt1$Fur_region)
p1_TC_anova2
plot(p1_TC_anova2, main = "TC [THg] by fur region")

#normality test of residuals
shapiro.test(residuals(p1_TC_anova2)) # Not normal

par(mfrow = c(1,1), las = 1, mar = c(5,5,5,5))
hist(residuals(p1_TC_anova2), breaks = 30, freq = FALSE, 
     main = "Histogram of PELT 1 TC [THg] one-way ANOVA residuals by fur region")
curve(dnorm(x, mean(residuals(p1_TC_anova2)), sd(residuals(p1_TC_anova2))), 
      add=TRUE, col="darkblue", lwd=2) # Approx normal


# Homoscedasticity test of residuals
leveneTest(pelt1$THg_TC ~ pelt1$Fur_region) # Homoscedastic

#Tukey's multiple comparison
p1_TC_posthoc2 = aov(pelt1$THg_TC ~ pelt1$Fur_region)
p1_TC_tukey2 = TukeyHSD(p1_TC_posthoc2)
plot(p1_TC_tukey2)

# UNDERCOAT
par(mfrow = c(2,2), mar = c(5,5,5,5))
p1_UC_anova2 = lm(pelt1$THg_UC ~ pelt1$Fur_region)
p1_UC_anova2
plot(p1_UC_anova2, main = "UC [THg] by fur region")

#normality test of residuals
shapiro.test(residuals(p1_UC_anova2)) # Not normal

par(mfrow = c(1,1), las = 1, mar = c(5,5,5,5))
hist(residuals(p1_UC_anova2), breaks = 30, freq = FALSE, 
     main = "Histogram of PELT 1 UC [THg] one-way ANOVA residuals by fur region")
curve(dnorm(x, mean(residuals(p1_UC_anova2)), sd(residuals(p1_UC_anova2))), 
      add=TRUE, col="darkblue", lwd=2) # Approximately normal

# Homoscedasticity test of residuals
leveneTest(pelt1$THg_UC ~ pelt1$Fur_region) # Homoscedastic

#Tukey's multiple comparison
p1_UC_posthoc2 = aov(pelt1$THg_UC ~ pelt1$Fur_region)
p1_UC_tukey2 = TukeyHSD(p1_UC_posthoc2)
plot(p1_UC_tukey2)