library(ggplot2)
library(dplyr)
library(HH)
library(knitr)

# Functions that implement theorems from linear model theory are
# defined in a separate file.
source("C:\\Users\\Todd\\R\\StatisticsProjects\\LinearModelFunctions.R")


# Two-Way ANOVA Analysis using R lm() function.
# Duplicate lm() results using linear model theory.
# See Christensen, "Plane Answers to Complex Questions", 2002

# Data from Christensen, table 7.4
obs = c(61.5, 68.2, 64.0 ,65.0 ,59.7, 55.0, 42.0, 60.2, 52.5, 61.8, 49.5, 52.7, 42.0, 54.0 ,61.0, 48.2, 39.6)
obs = c(obs, 60.3, 51.7, 49.3, 48.0, 50.8, 64.7, 61.7, 64.0, 62.0, 56.5, 59.0, 47.2, 53.0, 51.3, 40.5)
obs = c(obs, 37.0, 36.3, 68.0, 56.3, 69.8, 67.0, 39.7, 46.0, 61.3, 55.3, 55.7, 50.0, 43.8, 54.5)
obs = c(obs, 59.0, 57.4, 54.0, 47.0, 59.5, 52.8, 56.0, 45.2, 57.0, 61.4, 44.8, 51.5, 53.0, 42.0, 54.0)

# Rename matrix of observed weight gain
WeightGain = obs

Litter = c(rep("LA", 17), rep("LF", 15), rep("LI", 14), rep("LJ", 15))

Mother = c(rep("MA", 5), rep("MF", 3), rep("MI", 4), rep("MJ", 5))
Mother = c(Mother, rep("MA", 4), rep("MF", 5), rep("MI", 4), rep("MJ", 2))
Mother = c(Mother, rep("MA", 3), rep("MF", 3), rep("MI", 5), rep("MJ", 3))
Mother = c(Mother, rep("MA", 4), rep("MF", 3), rep("MI", 3), rep("MJ", 5))

# Put data into data frame.
dfRats = data.frame(WeightGain, Litter, Mother)

# ANOVA Type I
outI = lm(WeightGain ~ Litter * Mother, data = dfRats)
anova(outI)

# Display interaction plot.
dev.new() # New plot in separate window.
prettyInteractionPlot(dfRats, 1, 2, 3, "Litter:Mother Interaction")

# Duplicate anova(lm()) results using linear model theory.

# Put data into equivalent over-parameterized design matrix.
mRats = buildDesignAllMain(dfRats, 1, 2:3)
mRatsi = buildInteractionEffect(mRats, 3:6, 7:10)
mRats = cbind(mRats, mRatsi)

Y = mRats[, 1]
X = mRats[, -1]

# Duplicate anova(lm()) results using corresponding
# comparisons of full vs. reduced model.
outL = meanSquares(Y, X, 1:5, 1, 1:dim(X)[2], "R(L | u)")
outM = meanSquares(Y, X, 1:9, 1:5, 1:dim(X)[2], "R(M | u,L)")
outLM = meanSquares(Y, X, 1:dim(X)[2], 1:9, 1:dim(X)[2], "R(L:M | u,L,M)")
outE = meanSquares(Y, X, 1, 1, 1:dim(X)[2], "Residuals")
outWithI = rbind(outL, outM, outLM, outE)
outWithI

# Since interaction is not statistically significant,
# remove interaction from model.
out = lm(WeightGain ~ Litter + Mother, data = dfRats)
anova(out)

# We can use drop1() for this unbalanced model since
# there is no interaction term.
drop1(out, ~., test = "F")

# To be redundant, we can also us lm() and anova() to
# test the significance of the Mother effect via
# reduction in SSE.
outX = lm(WeightGain ~ Litter + Mother, data = dfRats)
out0 = lm(WeightGain ~ Litter, data = dfRats)
anova(out0, outX)

# Duplicate anova(lm()) results using corresponding
# comparisons of full vs. reduced model.
outL = meanSquares (Y, X, 1:9, c(1, 6:9), 1:9, "R(L | u,M)")
outM = meanSquares(Y, X, 1:9, 1:5, 1:9, "R(M | u,L)")
outE = meanSquares(Y, X, 1, 1, 1:9, "Residuals")
outNoI = rbind(outL, outM, outE)
outNoI

# Results above indicate Mother effect is statistically significant.

# Residual analysis to check ANOVA assumptions.
 
# Check assumption of normally distributed residuals.
shapiro.test(out$residuals)

dev.new() # New plot in separate window.
prettyQQ(out$residuals, "QQ Plot of Residuals")

# Brown-Forsyth test for homogeneous variances between categories of effect M
hov(WeightGain ~ Mother, data = dfRats)

# Graphical check of homogeneous variances assumption via
# categorized scatter plot of residuals vs. fitted values
dfOut = data.frame(out$fitted.values, out$residuals, rstudent(out), dfRats$Mother)
colnames(dfOut) = c("Fitted", "Residuals", "StudentizedRes", "Mother")

dev.new() # New plot in separate window.
prettyResVsFit(dfOut, 1, 2, FALSE, "Residuals vs. Fitted Values, Categorized by Mother", 4)

# Box plots to compare levels within Mother effect.
dev.new() # New plot in separate window.
prettyBoxPlot(dfRats, 1, 3, "Weight Gain Categorized by Mother")

# Post-hoc test for significant pair-wise differences
# between levels of Mother effect.
tukeyLevel = 0.95
tOutA = TukeyHSD(aov(out), "Mother", conf.level = tukeyLevel)
tOutB = round(tOutA$Mother, 3)
kable(tOutB)

# Check for outliers by reviewing studentized residuals.
dev.new() # New plot in separate window.
prettyResVsFit(dfOut, 1, 3, TRUE, "Studentized Residuals vs. Fitted Values")

