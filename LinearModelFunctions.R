library(stringr)
library(MASS)

# Implement functions from linear model theory.
# See Christensen, "Plane Answers to Complex Questions", 2002


# Create perpendicular project matrix with same column space as X.
# X: design matrix
# return: perpendicular project matrix onto C(X)
ppm = function(X)
{
  # Create orthonomal matrix with same column space as X.
  # Christensen, theorem B.44
  M = X %*% ginv(t(X) %*% X) %*% t(X)

  M
}

 
# Trace of perpendicular projection matrix.
# M: perpendicular project matrix.
# return: trace of M
tr = function(M)
{
  # Note that for a perpendicular project matrix, trace = rank.
  # Christensen definition B.27 and discussion of Colorllar B.49.
  sum(diag(M))
}


# Calc SS, df, MS, F, and p-value given full and reduced model
# using Christensen section 3.2.
# Y: matriX with observations (dependent variable)
# X: over-parameterized design matrix
# iFull: indices of design matrix for full model
# iReduced: indices of design matrix for reduced model
# iError: indices of design matrix for SSE
# id: name of effect being analyzed,
#     if id = "Residuals" then SSE, dfE, and MSE are returned.
# return: vector containing SS, df, MS, F, pVal (in order),
#         unless id = "Residuals" then vector containing
#         SSE, MSE, dfE (in order)
meanSquares = function(Y, X, iFull, iReduced, iError, id)
{
  # Design matrix used for error term.
  XE = X[, iError]
  
  # MS error
  ME = ppm(XE)
  I = diag(1, dim(ME))
  diffE = I - ME
  SSE = t(Y) %*% diffE %*% Y
  dfE = tr(diffE)
  MSE = SSE / dfE

  # Return MS effect or MS error depending on id argument.
  if(id != "Residuals")
  {
    # Design matrix for full model.
    XF = X[, iFull]
    
    # Design matrix for reduced model.
    XR = cbind(X[, iReduced])
    
    # MS effect
    MF = ppm(XF)
    MR = ppm(XR)
    diffA = MF - MR
    
    SSA = t(Y) %*% diffA %*% Y
    
    # Christensen, discussion of corollary B.49
    dfA = tr(diffA)
    
    MSA = SSA / dfA
    
    F = MSA / MSE
    pVal = pf(F, dfA, dfE, lower.tail = FALSE)
    
    out = c(round(dfA, 1), round(SSA, 3), round(MSA, 3), round(F, 4), round(pVal, 6))
  }
  else
  {
    out = c(round(dfE, 1), round(SSE, 3), round(MSE, 3), NA, NA)
  }
  
  rowsCols = list(id, c("df", "SS", "MS", "F", "pVal"))
  
  matrix(out, 1, 5, dimnames = rowsCols)
}


# Convert explanatory variable into over-parameterized design matrix.
# A: data frame of data
# effectIndex: index of explanatory variable in data frame
# return: over-parameterized design matrix for specified explanatory variable
buildDesignEffect = function(A, effectIndex)
{
  # Find number of levels of effect.
  lvls = sapply(A,  levels)[[effectIndex]]
  lvls
  
  # Create first column of design matrix for first level of effect.
  X = as.integer(A[effectIndex] == lvls[1])
  
  # Add columns to design matrix for subsequent levels of effect.
  for(i in 2:length(lvls))
  {  
    Xi = as.integer(A[effectIndex] == lvls[i])
    X = cbind(X, Xi)
  }
  
  # Make levels names column names.
  colnames(X) = lvls
  
  X
}


# Convert data frame into over-parameterized design matrix
# with main effects only.
# A: data frame of data
# yIndex: index of response (dependent) variable in data frame
# xIndices: indices of explanatory (independent) variables in data frame
# return: over-parameterized design matrix
buildDesignAllMain = function(A, yIndex, xIndices)
{
  # Build column containing response variable.
  Xy = A[, yIndex]
  
  # Build column for overall mean parameter.
  Xm = rep(1, dim(A)[1])
  
  X = cbind(Xy, Xm)
  
  # Build columns for categorical explanatory variables.
  numEffects = length(xIndices)
  
  for(k in 1:numEffects)
  {
    Xk = buildDesignEffect(A, xIndices[k])
    X = cbind(X, Xk)
  }
  
  X
}


# Build over-parameterized design matrix for interaction between two main effects.
# Z: over-parameterized design matrix with main effects.
# indicesA: indices of first main effect in specified design matrix
# indicesB: indices of second main effect in specified design matrix
# return: over-parameterized design matrix
buildInteractionEffect = function(Z, indicesA, indicesB)
{
  lengthA = length(indicesA)
  lengthB = length(indicesB)
  
  for(i in 1:lengthA)
  {
    for(j in 1:lengthB)
    {
      # Use first two letters of main effect names for name of interaction.
      nameA = substr(str_trim(colnames(Z)[indicesA[i]]), 1, 2)
      nameB = substr(str_trim(colnames(Z)[indicesB[j]]), 1, 2)
      cNames = paste(nameA, "|", nameB, sep = "")
      
      Xtemp = cbind(Z[, indicesA[i]] * Z[, indicesB[j]])
      colnames(Xtemp) = cNames
      
      if(i == 1 && j == 1)
      {
        X = Xtemp
      }
      else
      {
        X = cbind(X, Xtemp)
      }
    }
  }
  
  X
}


# Display pretty version of QQ plot.
# resids: residuals from a linear model
# hdr: plot title
prettyQQ = function(resids, hdr)
{
  # Use qqnorm() function to get x and y variables in QQ plot
  qn = qqnorm(resids, plot.it = FALSE)
  dfQn = data.frame(qn)
  
  # Create linear model of normal quantiles vs. residual observed quantiles
  qObs = quantile(x = resids, probs = c(0.25, 0.75))
  qModel = lm(qObs ~ qnorm(c(0.25, 0.75)))
  
  # QQ plot values
  qq = ggplot(dfQn, aes(x, y)) + geom_point(shape = 21, size = 4, color = "black", fill = "orange2")
  qq = qq + ggtitle(hdr) + labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
  # Add fit line.
  qq = qq + geom_abline(slope = qModel$coef[2], intercept = qModel$coef[1], size = 0.8, color = "darkolivegreen")
  qq = qq + theme_bw()
  print(qq)
}  


# Display pretty categorized box plot.
# dfData: data frame containing observations and categorical variable (in order)
# iObs: index of column with observations
# iVarX: index of column with x-axis variable
# hdr: plot title
prettyBoxPlot = function(dfData, iObs, iVarX, hdr)
{
  cNames = colnames(dfData)
  obsName = cNames[iObs]
  xName = cNames[iVarX]
  
  bp = ggplot(dfData, aes_q(as.name(xName), as.name(obsName)))
  bp = bp + stat_boxplot(geom = "errorbar", color = "darkolivegreen4", size = 0.8)
  bp = bp + geom_boxplot(outlier.size = 4, outlier.shape = 8, color = "darkolivegreen", size = 0.8)
  bp = bp + scale_fill_brewer(palette = "Set1") + theme_bw()
  bp = bp + ggtitle(hdr)
  
  print(bp)
}


# Display pretty categorized scatter plot of residuals vs. fitted values.
# dfData: data frame containing fitted values, residuals, and
#         categorical variable (in order)
# iFit: index of column with fitted values
# iRes: index of column with residuals
# bStudentized: TRUE => studentized residuals uncategorized
#               FALSE => residuals categorized
# hdr: plot title
# iCat: index of column with variable used to categorization
prettyResVsFit = function(dfData, iFit, iRes, bStudentized, hdr, iCat)
{
  cNames = colnames(dfData)
  fitName = cNames[iFit]
  resName = cNames[iRes]
 
  sp = ggplot(dfData, mapping=aes_q(as.name(fitName), as.name(resName)))
  sp = sp + ggtitle(hdr) + labs(x = fitName, y = resName)
  
  # Categorize only for non-studentized residuals.
  if(bStudentized)
  {
    sp = sp + geom_point(size = 4, shape = 21, color = "black", fill = "orange2")
  }
  else
  {
    catName = cNames[iCat]
    sp = sp + geom_point(size = 4, shape = 21, color = "black", aes_q(fill = as.name(catName)))
    sp = sp + scale_fill_brewer(palette = "Set1")
    }

  sp = sp + theme_bw()
  
  print(sp)
}


# Display pretty interaction plot.
# dfData: data frame containing observations and two categorical variables.
# iObs: index of column with observations, used for y-axis
# iVarX: index of column with first categorical variable, used for x-axis
# iGroup: index of column with second categorical variable, used to group
# hdr: title of plot
prettyInteractionPlot = function(dfData, iObs, iVarX, iGroup, hdr)
{
  iData = with(dfData, aggregate(dfData[, iObs], list(dfData[, iVarX], dfData[, iGroup]), mean))
  colnames(iData) = c(colnames(dfData)[c(2, 3)], "AvgWeightGain") 
  
  xName = colnames(iData)[1]
  yName = colnames(iData)[3]
  groupName = colnames(iData)[2]
  
  sp = ggplot(iData, aes_q(as.name(xName), as.name(yName), color = as.name(groupName)))

  sp = sp + geom_point(size = 4, shape = 21, color = "black", aes_q(fill = as.name(groupName)))
  sp = sp + geom_line(size = 1, aes_q(group = as.name(groupName)))
  sp = sp + scale_fill_brewer(palette = "Set1")
  sp = sp + scale_color_brewer(palette = "Set1") + theme_bw()
  sp = sp + ggtitle("Mother:Litter Interaction")
  print(sp) 
}


