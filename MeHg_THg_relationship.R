#############################################################################
# Analysis for MeHg and THg data
# Written in R Version 3.5.0
#############################################################################
# Load data
compiled = read.csv("pelts_compiled.csv")

# Load Libraries
library(ggplot2)
library(boot)
library(lmtest)
library(car)

############################################################
#### THg to MeHg relationship ####
comp_MeHg = na.omit(compiled)
#add column with MeHg/THg ratio in combined df
comp_MeHg$Percent_MeHg = (comp_MeHg$MeHg_TC/comp_MeHg$THg_TC)*100
subset_110 = comp_MeHg[which(comp_MeHg$Percent_MeHg <= 110),]

#Regression Coefficient for MeHg/THg < 110%
#Ensure normality
qqnorm(comp_MeHg$Percent_MeHg[which(comp_MeHg$Percent_MeHg <= 110)],
       ylab="Percent MeHg in THg (%)")
qqline(comp_MeHg$Percent_MeHg[which(comp_MeHg$Percent_MeHg <= 110)])

#Generate regression model with fixed intercept at 0
MeHg_THg_lm = lm(formula = MeHg_TC ~ 0 + THg_TC, 
                  data = comp_MeHg[which(comp_MeHg$Percent_MeHg <= 110),])
summary(MeHg_THg_lm)

par(mfrow = c(2, 2), las=1)
plot(MeHg_THg_lm)

# Test assumptions 
#Test of homoscedasticty
ncvTest(MeHg_THg_lm) #homoscedastic
# Test for autocorrelation
dwtest(MeHg_THg_lm) #no autocorrelation
# Test for linearity
resettest(MeHg_THg_lm) #linear
# Test for Normality
shapiro.test(resid(MeHg_THg_lm)) #violation
hist((resid(MeHg_THg_lm))) #visual assessment- approximately normal

# Bootstrap regression coefficient to estimate 95% CI
# function to obtain regression weights
bs = function(formula, data, indices) {
  d = data[indices,] # allows boot to select sample
  fit = lm(formula, data=d)
  return(coef(fit))
}
# bootstrapping with 1000 replications
results = boot(data = comp_MeHg[which(comp_MeHg$Percent_MeHg <= 110),], 
                statistic=bs,
                R=10000, formula = MeHg_TC ~ 0 + THg_TC)
results
plot(results)
#show 95% CI for bootstrapped regression coefficient
boot.ci(results, type = "all", index = 1)

# Plot regression
MeHg_THg_scatter = ggplot(data = subset_110, aes(x = THg_TC, y = MeHg_TC)) +
  stat_smooth(method = lm, se = TRUE) +
  geom_point() +
  ggtitle("[MeHg] vs [THg] relationship") +
  xlab(label = "THg concentration (ppm)") +
  ylab(label = "MeHg concentration (ppm)")
MeHg_THg_scatter



