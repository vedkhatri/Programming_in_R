###############################################################################
                                 #DATASET #1
###############################################################################

# To import the dataset,
arth <- read.csv("arthritis.csv", stringsAsFactors = TRUE)

# To check the normality of the dataset, we use the SHAPIRO-WILK TEST
shapiro.test(arth$tetanus_titer) #Tetanus_titer
shapiro.test(arth$pertussis_titer) #Pertussis_titer
shapiro.test(arth$age) #Age
shapiro.test(arth$age_at_vaccination) #Age_at_vaccination
#When p > 0.05, the datasets are normally-distributed.
#In this case, the normally-distributed datasets are,
#       Tetanus_titer (p-value = 0.2697)
#       Age (p-value = 0.7567)

#For normally-distributed datasets, we perform t-test to check whether
#there is a significant difference between the mean values
t.test(tetanus_titer ~ arthritis, data = arth)
#Because p-value < 0.05, there IS a significant difference between the mean
#values between tetanus_titer for group 0 and group 1.

#Similarly, t-test for age with respect to arthritis
t.test(age ~ arthritis, data = arth)
#Because p-value = 0.163 (> 0.05), there IS NOT a significant difference between
#the mean values between tetanus_titer for group 0 and group 1.

#For other two datasets, we use Wilcoxon rank sum test,
wilcox.test(pertussis_titer ~ arthritis, data = arth)
wilcox.test(age_at_vaccination ~ arthritis, data = arth)

#     CORRELATION TESTING
#For normally-distributed datasets, we use Pearson correlation testing
cor.test(arth$tetanus_titer, arth$age) # Negative weak correlation

#For non-normal distributions, we use the Spearman correlation testing
cor.test(arth$pertussis_titer, arth$age_at_vaccination, method = "spearman")
cor.test(arth$tetanus_titer, arth$age_at_vaccination, method = "spearman")
cor.test(arth$tetanus_titer, arth$age, method = "spearman")
# Negative weak correlations for all, expcept positive weak correlation for
#   Pertussis_titer ~ age_at_vaccination


# Considering arth$arthritis as a factor
arth$arthrts <- as.factor(arth$arthritis)

#Linear models
library(beeswarm)
a1 <- lm(arth$tetanus_titer ~ arth$rthrts) #Tetanus_titer ~ Arthrts
a2 <- lm(arth$pertussis_titer ~ arth$arthrts) #Pertussis_titer ~ Arthrts
a3 <- lm(arth$age ~ arth$arthrts) #Age ~ arthrts
a4 <- lm(arth$age_at_vaccination ~ arth$arthrts) #age_at_vaccination ~ arthrts


# To perform Levene Test
leveneTest(arth$tetanus_titer ~ arth$arthrts)
leveneTest(arth$pertussis_titer ~ arth$arthrts)
leveneTest(arth$age ~ arth$arthrts)
leveneTest(arth$age_at_vaccination ~ arth$arthrts)
# From the results, the p-value > 0.05 is for (differences of variances NOT significant)
#   Tetanus_titer ~ arthrts
#   Age ~ artrts
#   Age_by_vaccination ~ arthrts

#To confirm homogeneity of variances for Pertussis_titer ~ arthrts,
a5 <- lm(log(arth$pertussis_titer+1) ~ arth$arthrts)
leveneTest(log(arth$pertussis_titer+1) ~ arth$arthrts)

# To perform normality of residuals
shapiro.test(a1$residuals) # For Tetanus_titer ~ arthrts
shapiro.test(a5$residuals) # For Pertussis_titer ~ arthrts
shapiro.test(a3$residuals) # For Age ~ arthrts
shapiro.test(a4$residuals) # For Age_at_vaccination ~ arthrts
# The ones whose residuals ARE normally-distributed are
#   a3: Age ~ arthrts
#   a5: Pertussis_titer ~ arthrts
# The ones whose residuals ARE NOT normally-distributed are:
#   a1: Tetanus_titer ~ arthrts
#   a4: Age_at_vaccination ~ arthrts


########## Regression Analyses begins

# Linear models again, however, these are for IVs
ar1 <- lm(arth$tetanus_titer ~ arth$age)    # Tetanus_titer ~ Age
ar2 <- lm(arth$pertussis_titer ~ arth$age) # Pertussis_titer ~ Age
ar3 <- lm(arth$tetanus_titer ~ arth$age_at_vaccination)  #Tetanus_titer ~ Age_at_vaccination
ar4 <- lm(arth$pertussis_titer ~ arth$age_at_vaccination) #Pertussis_titer ~ Age_at_vaccination

# Plotting
plot(arth$tetanus_titer ~ arth$age) # For ar1
abline(ar1)
plot(arth$pertussis_titer ~ arth$age) # For ar2
abline(ar2)
plot(arth$tetanus_titer ~ arth$age_at_vaccination) # For ar3
abline(ar3)
plot(arth$pertussis_titer ~ arth$age_at_vaccination) # For ar4
abline(ar4)

# Getting the abline equation y = a + bx
summary(ar1)
summary(ar2)
summary(ar3)
summary(ar4)

## How to read the summary table
    # p-value < 0.05 indicates a significant relationship of the IV
    # Adjusted R-squared (when closer to 1), indicates a strong relationship
    # The overall regression is denoted by p-value in the bottom-most line


### Considering LOGISTIC PREGRESSION

library(lmtest)

# Creating generalised linear models
b1 <- glm(arthritis ~ tetanus_titer, data = arth, family = "binomial")
b2 <- glm(arthritis ~ pertussis_titer, data = arth, family = "binomial")
b3 <- glm(arthritis ~ age, data = arth, family = "binomial")
b4 <- glm(arthritis ~ age_at_vaccination, data = arth, family = "binomial")
b5 <- glm(arthritis ~ tetanus_titer+pertussis_titer, data = arth, family = "binomial")
b6 <- glm(arthritis ~ tetanus_titer+pertussis_titer+age+age_at_vaccination, data = arth, family = "binomial")


# Summary
summary(b1)
summary(b2)
summary(b3)
summary(b4)
summary(b5)
summary(b6)

# VIF
vif(b5)
vif(b6)

library(lmtest)
# Likelihood Ratio tests
lrtest(b1)
lrtest(b2)
lrtest(b3)
lrtest(b4)
lrtest(b5)
lrtest(b6)


###############################################################################
                                  #DATASET #2
###############################################################################

# To import the dataset,
bp <- read.csv("blood_pressure.csv", stringsAsFactors = TRUE)

# To check the normality of the dataset, we use the SHAPIRO-WILK TEST
shapiro.test(bp$syst_change) #Syst_change
shapiro.test(bp$PTH) #PTH
shapiro.test(bp$renin) #renin
shapiro.test(bp$Ca) #Ca
shapiro.test(bp$age) #age
# The datasets that ARE normally-distributed are
#   syst_change
#   Ca

# CORRELATION TESTING
  # For normally-distributed datasets
cor.test(bp$age, bp$syst_change)
plot(bp$syst_change ~ bp$Ca)
  # For the remaining datasets
cor.test(bp$PTH, bp$syst_change, method = "spearman")
plot(bp$syst_change ~ bp$PTH)
cor.test(bp$renin, bp$syst_change, method = "spearman")
plot(bp$syst_change ~ bp$renin)
cor.test(bp$age, bp$syst_change, method = "spearman")
plot(bp$syst_change ~ bp$age)

######### REGRESSION ANALYSES BEGINS

# Linear models
bp1 <- lm(bp$syst_change ~ bp$Ca)
bp2 <- lm(bp$syst_change ~ bp$PTH)
bp3 <- lm(bp$syst_change ~ bp$renin)
bp4 <- lm(bp$syst_change ~ bp$age)
bp5 <- lm(syst_change ~ Ca+PTH, data = bp)


# Plotting and ABLine verification
plot(bp$syst_change ~ bp$Ca)
abline(bp1)
plot(bp$syst_change ~ bp$PTH)
abline(bp2)
plot(bp$syst_change ~ bp$renin)
abline(bp3)
plot(bp$syst_change ~ bp$age)
abline(bp4)

# Check for Variances
plot(bp1, which = 1)
plot(bp2, which = 1)
plot(bp3, which = 1)
plot(bp4, which = 1)
plot(bp5, which = 1)

# Check for residuals to be normally-distributed by plotting
plot(bp1, which = 2)
plot(bp2, which = 2)
plot(bp3, which = 2)
plot(bp4, which = 2)
plot(bp5, which = 2)

# Check for residuals to be normally-distributed by Shapiro-test
shapiro.test(bp1$residuals)
shapiro.test(bp2$residuals)
shapiro.test(bp3$residuals)
shapiro.test(bp4$residuals)
shapiro.test(bp5$residuals)

# Summary
summary(bp1)
summary(bp2)
summary(bp3)
summary(bp4)
summary(bp5)

## How to read the summary table
  # p-value < 0.05 indicates a significant relationship of the IV
  # Adjusted R-squared (when closer to 1), indicates a strong relationship
  # The overall regression is denoted by p-value in the bottom-most line

# To assess the multicollinarity by using Variance Inflation Factors (vif)
vif(bp5)


###############################################################################
                                  #DATASET #3
###############################################################################

# To import the dataset,
gw <- read.csv("groundwater.csv", stringsAsFactors = TRUE)

# To check the normality of the dataset, we use the SHAPIRO-WILK TEST
shapiro.test(gw$coliform_CFU) #coliform_CFU

# To check the differences between mean by one-way Anova as the distribution is
# NOT NORMALLY-DISTRIBUTED

# The first assumption is to check the difference in variances is not significant
leveneTest(gw$coliform_CFU ~ gw$area)
# Pr(>F) is the probability value.
# Because p-value > 0.05, the difference of variances is not significant.

#The second assumption is that the residual values should be normally-distrbuted
g1 <- lm(gw$coliform_CFU ~ gw$area) # Linear model
shapiro.test(g1$residuals)
# Because the p-value is NOT > 0.05, the residuals are NOT normally-distributed
# Therefore, one-way Anova conditions are NOT met.

# Alternatively, these could be plotted to check the conditions are met or not
plot(g1, which = 1) # for VARIANCES
plot(g1, which = 2) # for NORMALITY of residuals

# So, we use Kruskal-Wallis test and post-hoc Dunn Test
kruskal.test(gw$coliform_CFU ~ gw$area)
# Because p-value < 0.05, the difference between the MEDIAN values IS significant

# Therefore, we perform the Dunn Test
library(dunn.test)
dunn.test(gw$coliform_CFU, gw$area, kw = FALSE, method = 'bh')
# The data with significant differences in their MEDIANS are:
#   Callenda ~ Fairmile
#   Esk sand ~ Torness




################################################################################
                               # DATASET 4
################################################################################

ser <- read.csv("sertoli.csv", stringsAsFactors = TRUE)

## Plotting for assumptions
s1 <- lm(OD ~ conc*temp, data = ser)
plot(s1, which = 1)
plot(s1, which = 2)

leveneTest(OD ~ temp*conc, data = ser)
# Anova
Anova(s1, type = 3)
TukeyHSD(aov(s1))