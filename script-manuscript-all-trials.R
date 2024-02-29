# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

## Prepare environment ----------------------------------------------------------------------------------
rm(list=ls())

# Package loading
library(dplyr)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(pROC)

# Code sourcing
source("functions-process-results.R")

# Data sourcing
load("results-logistic-effect0-n15000.Rdata")
data.i15000 <- SimulationData

load("results-logistic-effect30-n15000.Rdata")
data.ii15000 <- SimulationData

load("results-hill-effect0-n15000.Rdata")
data.iii15000 <- SimulationData

load("results-hill-effect30-n15000.Rdata")
data.iv15000 <- SimulationData

## Analyze simulation results ---------------------------------------------------------------------------

i15000 <- analyzeSimulation(data.i15000)

ii15000 <- analyzeSimulation(data.ii15000)

iii15000 <- analyzeSimulation(data.iii15000)

iv15000 <- analyzeSimulation(data.iv15000)

### Numerical results, Tables: --------------------------------------------------------------------------

## Table 1 ----
# i
data.i15000$truth$vaccinated$N + data.i15000$truth$control$N # total number of subjects
data.i15000$truth$PoDModel$predictionFun # true PoD model
i15000$trueVE # true VE (group A = older; group B = younger)
data.i15000$truth$vaccinated # vaccinated titers
data.i15000$truth$control # control titers
data.i15000$truth$PoDModel$valueParams # true PoD parameters, need to be transformed to the parameters in the manuscript (different PoD parametrization)

# ii
data.ii15000$truth$vaccinated$N + data.ii15000$truth$control$N # total number of subjects
data.ii15000$truth$PoDModel$predictionFun # true PoD model
ii15000$trueVE # true VE (group A = older; group B = younger)
data.ii15000$truth$vaccinated # vaccinated titers
data.ii15000$truth$control # control titers
data.ii15000$truth$PoDModel$valueParams # true PoD parameters, need to be transformed to the parameters in the manuscript (different PoD parametrization)

# iii
data.iii15000$truth$vaccinated$N + data.iii15000$truth$control$N # total number of subjects
data.iii15000$truth$PoDModel$predictionFun # true PoD model
iii15000$trueVE # true VE (group A = older; group B = younger)
data.iii15000$truth$vaccinated # vaccinated titers
data.iii15000$truth$control # control titers
data.iii15000$truth$PoDModel$valueParams # true PoD parameters

# iv
data.iv15000$truth$vaccinated$N + data.iv15000$truth$control$N # total number of subjects
data.iv15000$truth$PoDModel$predictionFun # true PoD model
iv15000$trueVE # true VE (group A = older; group B = younger)
data.iv15000$truth$vaccinated # vaccinated titers
data.iv15000$truth$control # control titers
data.iv15000$truth$PoDModel$valueParams # true PoD parameters

## Table 2 ----
# False positive and false negative rates 
FP.CoP.logistic <- i15000$positive$CoPbased.linear
FP.typical.logistic <- i15000$positive$typical
FN.CoP.logistic <- 1-ii15000$positive$CoPbased.linear
FN.typical.logistic <- 1-ii15000$positive$typical

FP.CoP.hill <- iii15000$positive$CoPbased.linear
FP.typical.hill <- iii15000$positive$typical
FN.CoP.hill <- 1-iv15000$positive$CoPbased.linear
FN.typical.hill <- 1-iv15000$positive$typical

# Positive predictive value and negative predictive value 
PPV.CoP.logistic <- (1-FN.CoP.logistic)/((1-FN.CoP.logistic)+FP.CoP.logistic) # TP/(TP+FP)
PPV.typical.logistic <- (1-FN.typical.logistic)/((1-FN.typical.logistic)+FP.typical.logistic)
NPV.CoP.logistic <- (1-FP.CoP.logistic)/((1-FP.CoP.logistic)+FN.CoP.logistic) # TN/(TN+FN)
NPV.typical.logistic <- (1-FP.typical.logistic)/((1-FP.typical.logistic)+FN.typical.logistic)

PPV.CoP.hill <- (1-FN.CoP.hill)/((1-FN.CoP.hill)+FP.CoP.hill) # TP/(TP+FP)
PPV.typical.hill <- (1-FN.typical.hill)/((1-FN.typical.hill)+FP.typical.hill)
NPV.CoP.hill <- (1-FP.CoP.hill)/((1-FP.CoP.hill)+FN.CoP.hill) # TN/(TN+FN)
NPV.typical.hill <- (1-FP.typical.hill)/((1-FP.typical.hill)+FN.typical.hill)

# i
1-FP.typical.logistic # true negative, typical
1-FP.CoP.logistic # true negative, CoP-based

# ii
1-FN.typical.logistic # true positive, typical
1-FN.CoP.logistic # true positive, CoP-based

# i & ii
PPV.typical.logistic # positive predictive value, typical
PPV.CoP.logistic # positive predictive value, CoP-based
NPV.typical.logistic # negative predictive value, typical
NPV.CoP.logistic  # negative predictive value, CoP-based

# iii
1-FP.typical.hill # true negative, typical
1-FP.CoP.hill # true negative, CoP-based

# iv
1-FN.typical.hill # true positive, typical
1-FN.CoP.hill # true positive, CoP-based

# iii & iv
PPV.typical.hill # positive predictive value, typical
PPV.CoP.hill # positive predictive value, CoP-based
NPV.typical.hill # negative predictive value, typical
NPV.CoP.hill  # negative predictive value, CoP-based

# ROC curves 
outcome <- c(rep(1,5000), rep(0,5000))
typicalTrial.ROC.logistic.CoP <- roc(outcome, c(ii15000$predictors.roc$CoP, i15000$predictors.roc$CoP))
typicalTrial.ROC.logistic.typical <- roc(outcome, c(ii15000$predictors.roc$typical, i15000$predictors.roc$typical))
typicalTrial.ROC.hill.CoP <- roc(outcome, c(iv15000$predictors.roc$CoP, iii15000$predictors.roc$CoP))
typicalTrial.ROC.hill.typical <- roc(outcome, c(iv15000$predictors.roc$typical, iii15000$predictors.roc$typical))

# i & ii
roc.test(typicalTrial.ROC.logistic.CoP, typicalTrial.ROC.logistic.typical) # AUC ROC, CoP-based and typical

# iii & iv
roc.test(typicalTrial.ROC.hill.CoP, typicalTrial.ROC.hill.typical) # AUC ROC, CoP-based and typical

## Table 3 ----
# i
i15000$selected.typical.model # number of simulated trials (add = model without interaction; multi = model with interaction), typical approach
i15000$selected.typical.model/5000*100 # % of simulated trials (add = model without interaction; multi = model with interaction), typical approach
i15000$selected.CoP.model  # number of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach
i15000$selected.CoP.model/5000*100  # % of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach

# ii
ii15000$selected.typical.model # number of simulated trials (add = model without interaction; multi = model with interaction), typical approach
ii15000$selected.typical.model/5000*100 # % of simulated trials (add = model without interaction; multi = model with interaction), typical approach
ii15000$selected.CoP.model  # number of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach
ii15000$selected.CoP.model/5000*100  # % of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach

# iii
iii15000$selected.typical.model # number of simulated trials (add = model without interaction; multi = model with interaction), typical approach
iii15000$selected.typical.model/5000*100 # % of simulated trials (add = model without interaction; multi = model with interaction), typical approach
iii15000$selected.CoP.model  # number of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach
iii15000$selected.CoP.model/5000*100  # % of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach

# iv
iv15000$selected.typical.model # number of simulated trials (add = model without interaction; multi = model with interaction), typical approach
iv15000$selected.typical.model/5000*100 # % of simulated trials (add = model without interaction; multi = model with interaction), typical approach
iv15000$selected.CoP.model  # number of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach
iv15000$selected.CoP.model/5000*100  # % of simulated trials (add = model without interaction; multi = model with interaction), CoP-based approach

## Supplementary Table S1 ----
# i
# MSE, younger
mse.selected.CoP.groupB <- 1/length(i15000$Results$selected.CoP.VE.groupB)*sum((i15000$Results$selected.CoP.VE.groupB*100-i15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.typical.groupB <- 1/length(i15000$Results$selected.typical.VE.groupB)*sum((i15000$Results$selected.typical.VE.groupB*100-i15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.casecount.groupB <- 1/length(i15000$Results$casecountVE.groupB)*sum((i15000$Results$casecountVE.groupB*100-i15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.CoP.groupB
mse.selected.typical.groupB 
mse.selected.casecount.groupB

# MSE, older
mse.selected.CoP.groupA <- 1/length(i15000$Results$selected.CoP.VE.groupA)*sum((i15000$Results$selected.CoP.VE.groupA*100-i15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.typical.groupA <- 1/length(i15000$Results$selected.typical.VE.groupA)*sum((i15000$Results$selected.typical.VE.groupA*100-i15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.casecount.groupA <- 1/length(i15000$Results$casecountVE.groupA)*sum((i15000$Results$casecountVE.groupA*100-i15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.CoP.groupA
mse.selected.typical.groupA 
mse.selected.casecount.groupA

# MSE, overall
mse.selected.CoP.all <- 1/length(i15000$Results$selected.CoP.VE.all)*sum((i15000$Results$selected.CoP.VE.all*100-i15000$Results$trueVE.all[1]*100)^2)
mse.selected.typical.all <- 1/length(i15000$Results$selected.typical.VE.all)*sum((i15000$Results$selected.typical.VE.all*100-i15000$Results$trueVE.all[1]*100)^2)
mse.selected.casecount.all <- 1/length(i15000$Results$casecountVE.all)*sum((i15000$Results$casecountVE.all*100-i15000$Results$trueVE.all[1]*100)^2)
mse.selected.CoP.all
mse.selected.typical.all 
mse.selected.casecount.all

# ii
# MSE, younger
mse.selected.CoP.groupB <- 1/length(ii15000$Results$selected.CoP.VE.groupB)*sum((ii15000$Results$selected.CoP.VE.groupB*100-ii15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.typical.groupB <- 1/length(ii15000$Results$selected.typical.VE.groupB)*sum((ii15000$Results$selected.typical.VE.groupB*100-ii15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.casecount.groupB <- 1/length(ii15000$Results$casecountVE.groupB)*sum((ii15000$Results$casecountVE.groupB*100-ii15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.CoP.groupB
mse.selected.typical.groupB 
mse.selected.casecount.groupB

# MSE, older
mse.selected.CoP.groupA <- 1/length(ii15000$Results$selected.CoP.VE.groupA)*sum((ii15000$Results$selected.CoP.VE.groupA*100-ii15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.typical.groupA <- 1/length(ii15000$Results$selected.typical.VE.groupA)*sum((ii15000$Results$selected.typical.VE.groupA*100-ii15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.casecount.groupA <- 1/length(ii15000$Results$casecountVE.groupA)*sum((ii15000$Results$casecountVE.groupA*100-ii15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.CoP.groupA
mse.selected.typical.groupA 
mse.selected.casecount.groupA

# MSE, overall
mse.selected.CoP.all <- 1/length(ii15000$Results$selected.CoP.VE.all)*sum((ii15000$Results$selected.CoP.VE.all*100-ii15000$Results$trueVE.all[1]*100)^2)
mse.selected.typical.all <- 1/length(ii15000$Results$selected.typical.VE.all)*sum((ii15000$Results$selected.typical.VE.all*100-ii15000$Results$trueVE.all[1]*100)^2)
mse.selected.casecount.all <- 1/length(ii15000$Results$casecountVE.all)*sum((ii15000$Results$casecountVE.all*100-ii15000$Results$trueVE.all[1]*100)^2)
mse.selected.CoP.all
mse.selected.typical.all 
mse.selected.casecount.all

# iii
# MSE, younger
mse.selected.CoP.groupB <- 1/length(iii15000$Results$selected.CoP.VE.groupB)*sum((iii15000$Results$selected.CoP.VE.groupB*100-iii15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.typical.groupB <- 1/length(iii15000$Results$selected.typical.VE.groupB)*sum((iii15000$Results$selected.typical.VE.groupB*100-iii15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.casecount.groupB <- 1/length(iii15000$Results$casecountVE.groupB)*sum((iii15000$Results$casecountVE.groupB*100-iii15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.CoP.groupB
mse.selected.typical.groupB 
mse.selected.casecount.groupB

# MSE, older
mse.selected.CoP.groupA <- 1/length(iii15000$Results$selected.CoP.VE.groupA)*sum((iii15000$Results$selected.CoP.VE.groupA*100-iii15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.typical.groupA <- 1/length(iii15000$Results$selected.typical.VE.groupA)*sum((iii15000$Results$selected.typical.VE.groupA*100-iii15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.casecount.groupA <- 1/length(iii15000$Results$casecountVE.groupA)*sum((iii15000$Results$casecountVE.groupA*100-iii15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.CoP.groupA
mse.selected.typical.groupA 
mse.selected.casecount.groupA

# MSE, overall
mse.selected.CoP.all <- 1/length(iii15000$Results$selected.CoP.VE.all)*sum((iii15000$Results$selected.CoP.VE.all*100-iii15000$Results$trueVE.all[1]*100)^2)
mse.selected.typical.all <- 1/length(iii15000$Results$selected.typical.VE.all)*sum((iii15000$Results$selected.typical.VE.all*100-iii15000$Results$trueVE.all[1]*100)^2)
mse.selected.casecount.all <- 1/length(iii15000$Results$casecountVE.all)*sum((iii15000$Results$casecountVE.all*100-iii15000$Results$trueVE.all[1]*100)^2)
mse.selected.CoP.all
mse.selected.typical.all 
mse.selected.casecount.all

# iv
# MSE, younger
mse.selected.CoP.groupB <- 1/length(iv15000$Results$selected.CoP.VE.groupB)*sum((iv15000$Results$selected.CoP.VE.groupB*100-iv15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.typical.groupB <- 1/length(iv15000$Results$selected.typical.VE.groupB)*sum((iv15000$Results$selected.typical.VE.groupB*100-iv15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.casecount.groupB <- 1/length(iv15000$Results$casecountVE.groupB)*sum((iv15000$Results$casecountVE.groupB*100-iv15000$Results$trueVE.groupB[1]*100)^2)
mse.selected.CoP.groupB
mse.selected.typical.groupB 
mse.selected.casecount.groupB

# MSE, older
mse.selected.CoP.groupA <- 1/length(iv15000$Results$selected.CoP.VE.groupA)*sum((iv15000$Results$selected.CoP.VE.groupA*100-iv15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.typical.groupA <- 1/length(iv15000$Results$selected.typical.VE.groupA)*sum((iv15000$Results$selected.typical.VE.groupA*100-iv15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.casecount.groupA <- 1/length(iv15000$Results$casecountVE.groupA)*sum((iv15000$Results$casecountVE.groupA*100-iv15000$Results$trueVE.groupA[1]*100)^2)
mse.selected.CoP.groupA
mse.selected.typical.groupA 
mse.selected.casecount.groupA

# MSE, overall
mse.selected.CoP.all <- 1/length(iv15000$Results$selected.CoP.VE.all)*sum((iv15000$Results$selected.CoP.VE.all*100-iv15000$Results$trueVE.all[1]*100)^2)
mse.selected.typical.all <- 1/length(iv15000$Results$selected.typical.VE.all)*sum((iv15000$Results$selected.typical.VE.all*100-iv15000$Results$trueVE.all[1]*100)^2)
mse.selected.casecount.all <- 1/length(iv15000$Results$casecountVE.all)*sum((iv15000$Results$casecountVE.all*100-iv15000$Results$trueVE.all[1]*100)^2)
mse.selected.CoP.all
mse.selected.typical.all 
mse.selected.casecount.all

## Supplementary Table S2 ----
# i
i15000$coverage$typical.add$all*100 # overall coverage probability (%), typical, without interaction
i15000$coverage$typical.multi$all*100 # overall coverage probability (%), typical, with interaction
i15000$coverage$CoP.linear.add$all*100 # overall coverage probability (%), CoP-based, linear titer, without interaction
i15000$coverage$CoP.linear.multi$all*100 # overall coverage probability (%), CoP-based, linear titer, with interaction
i15000$coverage$CoP.quad.add$all*100 # overall coverage probability (%), CoP-based, quadratic titer, without interaction
i15000$coverage$CoP.quad.multi$all*100 # overall coverage probability (%), CoP-based, quadratic titer, with interaction

# ii
ii15000$coverage$typical.add$all*100 # overall coverage probability (%), typical, without interaction
ii15000$coverage$typical.multi$all*100 # overall coverage probability (%), typical, with interaction
ii15000$coverage$CoP.linear.add$all*100 # overall coverage probability (%), CoP-based, linear titer, without interaction
ii15000$coverage$CoP.linear.multi$all*100 # overall coverage probability (%), CoP-based, linear titer, with interaction
ii15000$coverage$CoP.quad.add$all*100 # overall coverage probability (%), CoP-based, quadratic titer, without interaction
ii15000$coverage$CoP.quad.multi$all*100 # overall coverage probability (%), CoP-based, quadratic titer, with interaction

# iii
iii15000$coverage$typical.add$all*100 # overall coverage probability (%), typical, without interaction
iii15000$coverage$typical.multi$all*100 # overall coverage probability (%), typical, with interaction
iii15000$coverage$CoP.linear.add$all*100 # overall coverage probability (%), CoP-based, linear titer, without interaction
iii15000$coverage$CoP.linear.multi$all*100 # overall coverage probability (%), CoP-based, linear titer, with interaction
iii15000$coverage$CoP.quad.add$all*100 # overall coverage probability (%), CoP-based, quadratic titer, without interaction
iii15000$coverage$CoP.quad.multi$all*100 # overall coverage probability (%), CoP-based, quadratic titer, with interaction

# iv
iv15000$coverage$typical.add$all*100 # overall coverage probability (%), typical, without interaction
iv15000$coverage$typical.multi$all*100 # overall coverage probability (%), typical, with interaction
iv15000$coverage$CoP.linear.add$all*100 # overall coverage probability (%), CoP-based, linear titer, without interaction
iv15000$coverage$CoP.linear.multi$all*100 # overall coverage probability (%), CoP-based, linear titer, with interaction
iv15000$coverage$CoP.quad.add$all*100 # overall coverage probability (%), CoP-based, quadratic titer, without interaction
iv15000$coverage$CoP.quad.multi$all*100 # overall coverage probability (%), CoP-based, quadratic titer, with interaction

## Supplementary Table S3 ----
# i
i15000$coverage$typical.add$groupB*100 # younger coverage probability (%), typical, without interaction
i15000$coverage$typical.multi$groupB*100 # younger coverage probability (%), typical, with interaction
i15000$coverage$CoP.linear.add$groupB*100 # younger coverage probability (%), CoP-based, linear titer, without interaction
i15000$coverage$CoP.linear.multi$groupB*100 # younger coverage probability (%), CoP-based, linear titer, with interaction
i15000$coverage$CoP.quad.add$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, without interaction
i15000$coverage$CoP.quad.multi$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, with interaction

# ii
ii15000$coverage$typical.add$groupB*100 # younger coverage probability (%), typical, without interaction
ii15000$coverage$typical.multi$groupB*100 # younger coverage probability (%), typical, with interaction
ii15000$coverage$CoP.linear.add$groupB*100 # younger coverage probability (%), CoP-based, linear titer, without interaction
ii15000$coverage$CoP.linear.multi$groupB*100 # younger coverage probability (%), CoP-based, linear titer, with interaction
ii15000$coverage$CoP.quad.add$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, without interaction
ii15000$coverage$CoP.quad.multi$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, with interaction

# iii
iii15000$coverage$typical.add$groupB*100 # younger coverage probability (%), typical, without interaction
iii15000$coverage$typical.multi$groupB*100 # younger coverage probability (%), typical, with interaction
iii15000$coverage$CoP.linear.add$groupB*100 # younger coverage probability (%), CoP-based, linear titer, without interaction
iii15000$coverage$CoP.linear.multi$groupB*100 # younger coverage probability (%), CoP-based, linear titer, with interaction
iii15000$coverage$CoP.quad.add$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, without interaction
iii15000$coverage$CoP.quad.multi$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, with interaction

# iv
iv15000$coverage$typical.add$groupB*100 # younger coverage probability (%), typical, without interaction
iv15000$coverage$typical.multi$groupB*100 # younger coverage probability (%), typical, with interaction
iv15000$coverage$CoP.linear.add$groupB*100 # younger coverage probability (%), CoP-based, linear titer, without interaction
iv15000$coverage$CoP.linear.multi$groupB*100 # younger coverage probability (%), CoP-based, linear titer, with interaction
iv15000$coverage$CoP.quad.add$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, without interaction
iv15000$coverage$CoP.quad.multi$groupB*100 # younger coverage probability (%), CoP-based, quadratic titer, with interaction

## Supplementary Table S4 ----
# i
i15000$coverage$typical.add$groupA*100 # older coverage probability (%), typical, without interaction
i15000$coverage$typical.multi$groupA*100 # older coverage probability (%), typical, with interaction
i15000$coverage$CoP.linear.add$groupA*100 # older coverage probability (%), CoP-based, linear titer, without interaction
i15000$coverage$CoP.linear.multi$groupA*100 # older coverage probability (%), CoP-based, linear titer, with interaction
i15000$coverage$CoP.quad.add$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, without interaction
i15000$coverage$CoP.quad.multi$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, with interaction

# ii
ii15000$coverage$typical.add$groupA*100 # older coverage probability (%), typical, without interaction
ii15000$coverage$typical.multi$groupA*100 # older coverage probability (%), typical, with interaction
ii15000$coverage$CoP.linear.add$groupA*100 # older coverage probability (%), CoP-based, linear titer, without interaction
ii15000$coverage$CoP.linear.multi$groupA*100 # older coverage probability (%), CoP-based, linear titer, with interaction
ii15000$coverage$CoP.quad.add$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, without interaction
ii15000$coverage$CoP.quad.multi$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, with interaction

# iii
iii15000$coverage$typical.add$groupA*100 # older coverage probability (%), typical, without interaction
iii15000$coverage$typical.multi$groupA*100 # older coverage probability (%), typical, with interaction
iii15000$coverage$CoP.linear.add$groupA*100 # older coverage probability (%), CoP-based, linear titer, without interaction
iii15000$coverage$CoP.linear.multi$groupA*100 # older coverage probability (%), CoP-based, linear titer, with interaction
iii15000$coverage$CoP.quad.add$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, without interaction
iii15000$coverage$CoP.quad.multi$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, with interaction

# iv
iv15000$coverage$typical.add$groupA*100 # older coverage probability (%), typical, without interaction
iv15000$coverage$typical.multi$groupA*100 # older coverage probability (%), typical, with interaction
iv15000$coverage$CoP.linear.add$groupA*100 # older coverage probability (%), CoP-based, linear titer, without interaction
iv15000$coverage$CoP.linear.multi$groupA*100 # older coverage probability (%), CoP-based, linear titer, with interaction
iv15000$coverage$CoP.quad.add$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, without interaction
iv15000$coverage$CoP.quad.multi$groupA*100 # older coverage probability (%), CoP-based, quadratic titer, with interaction

## Supplementary Table S5 ----

# i
# MSE, younger
mse.quad.multi.CoP.groupB <- 1/length(i15000$Results$logisticQuadCoPVE.groupB.multi)*sum((i15000$Results$logisticQuadCoPVE.groupB.multi*100-i15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.typical.groupB <- 1/length(i15000$Results$logisticVE.groupB.multi)*sum((i15000$Results$logisticVE.groupB.multi*100-i15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.casecount.groupB <- 1/length(i15000$Results$casecountVE.groupB)*sum((i15000$Results$casecountVE.groupB*100-i15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.CoP.groupB
mse.quad.multi.typical.groupB
mse.quad.multi.casecount.groupB

# MSE, older
mse.quad.multi.CoP.groupA <- 1/length(i15000$Results$logisticQuadCoPVE.groupA.multi)*sum((i15000$Results$logisticQuadCoPVE.groupA.multi*100-i15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.typical.groupA <- 1/length(i15000$Results$logisticVE.groupA.multi)*sum((i15000$Results$logisticVE.groupA.multi*100-i15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.casecount.groupA <- 1/length(i15000$Results$casecountVE.groupA)*sum((i15000$Results$casecountVE.groupA*100-i15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.CoP.groupA
mse.quad.multi.typical.groupA 
mse.quad.multi.casecount.groupA

# MSE, overall
mse.quad.multi.CoP.all <- 1/length(i15000$Results$logisticQuadCoPVE.all.multi)*sum((i15000$Results$logisticQuadCoPVE.all.multi*100-i15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.typical.all <- 1/length(i15000$Results$logisticVE.all.multi)*sum((i15000$Results$logisticVE.all.multi*100-i15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.casecount.all <- 1/length(i15000$Results$casecountVE.all)*sum((i15000$Results$casecountVE.all*100-i15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.CoP.all
mse.quad.multi.typical.all
mse.quad.multi.casecount.all

# ii
# MSE, younger
mse.quad.multi.CoP.groupB <- 1/length(ii15000$Results$logisticQuadCoPVE.groupB.multi)*sum((ii15000$Results$logisticQuadCoPVE.groupB.multi*100-ii15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.typical.groupB <- 1/length(ii15000$Results$logisticVE.groupB.multi)*sum((ii15000$Results$logisticVE.groupB.multi*100-ii15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.casecount.groupB <- 1/length(ii15000$Results$casecountVE.groupB)*sum((ii15000$Results$casecountVE.groupB*100-ii15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.CoP.groupB
mse.quad.multi.typical.groupB
mse.quad.multi.casecount.groupB

# MSE, older
mse.quad.multi.CoP.groupA <- 1/length(ii15000$Results$logisticQuadCoPVE.groupA.multi)*sum((ii15000$Results$logisticQuadCoPVE.groupA.multi*100-ii15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.typical.groupA <- 1/length(ii15000$Results$logisticVE.groupA.multi)*sum((ii15000$Results$logisticVE.groupA.multi*100-ii15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.casecount.groupA <- 1/length(ii15000$Results$casecountVE.groupA)*sum((ii15000$Results$casecountVE.groupA*100-ii15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.CoP.groupA
mse.quad.multi.typical.groupA 
mse.quad.multi.casecount.groupA

# MSE, overall
mse.quad.multi.CoP.all <- 1/length(ii15000$Results$logisticQuadCoPVE.all.multi)*sum((ii15000$Results$logisticQuadCoPVE.all.multi*100-ii15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.typical.all <- 1/length(ii15000$Results$logisticVE.all.multi)*sum((ii15000$Results$logisticVE.all.multi*100-ii15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.casecount.all <- 1/length(ii15000$Results$casecountVE.all)*sum((ii15000$Results$casecountVE.all*100-ii15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.CoP.all
mse.quad.multi.typical.all
mse.quad.multi.casecount.all
 
# iii
# MSE, younger
mse.quad.multi.CoP.groupB <- 1/length(iii15000$Results$logisticQuadCoPVE.groupB.multi)*sum((iii15000$Results$logisticQuadCoPVE.groupB.multi*100-iii15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.typical.groupB <- 1/length(iii15000$Results$logisticVE.groupB.multi)*sum((iii15000$Results$logisticVE.groupB.multi*100-iii15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.casecount.groupB <- 1/length(iii15000$Results$casecountVE.groupB)*sum((iii15000$Results$casecountVE.groupB*100-iii15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.CoP.groupB
mse.quad.multi.typical.groupB
mse.quad.multi.casecount.groupB

# MSE, older
mse.quad.multi.CoP.groupA <- 1/length(iii15000$Results$logisticQuadCoPVE.groupA.multi)*sum((iii15000$Results$logisticQuadCoPVE.groupA.multi*100-iii15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.typical.groupA <- 1/length(iii15000$Results$logisticVE.groupA.multi)*sum((iii15000$Results$logisticVE.groupA.multi*100-iii15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.casecount.groupA <- 1/length(iii15000$Results$casecountVE.groupA)*sum((iii15000$Results$casecountVE.groupA*100-iii15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.CoP.groupA
mse.quad.multi.typical.groupA 
mse.quad.multi.casecount.groupA

# MSE, overall
mse.quad.multi.CoP.all <- 1/length(iii15000$Results$logisticQuadCoPVE.all.multi)*sum((iii15000$Results$logisticQuadCoPVE.all.multi*100-iii15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.typical.all <- 1/length(iii15000$Results$logisticVE.all.multi)*sum((iii15000$Results$logisticVE.all.multi*100-iii15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.casecount.all <- 1/length(iii15000$Results$casecountVE.all)*sum((iii15000$Results$casecountVE.all*100-iii15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.CoP.all
mse.quad.multi.typical.all
mse.quad.multi.casecount.all

# iv
# MSE, younger
mse.quad.multi.CoP.groupB <- 1/length(iv15000$Results$logisticQuadCoPVE.groupB.multi)*sum((iv15000$Results$logisticQuadCoPVE.groupB.multi*100-iv15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.typical.groupB <- 1/length(iv15000$Results$logisticVE.groupB.multi)*sum((iv15000$Results$logisticVE.groupB.multi*100-iv15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.casecount.groupB <- 1/length(iv15000$Results$casecountVE.groupB)*sum((iv15000$Results$casecountVE.groupB*100-iv15000$Results$trueVE.groupB[1]*100)^2)
mse.quad.multi.CoP.groupB
mse.quad.multi.typical.groupB
mse.quad.multi.casecount.groupB

# MSE, older
mse.quad.multi.CoP.groupA <- 1/length(iv15000$Results$logisticQuadCoPVE.groupA.multi)*sum((iv15000$Results$logisticQuadCoPVE.groupA.multi*100-iv15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.typical.groupA <- 1/length(iv15000$Results$logisticVE.groupA.multi)*sum((iv15000$Results$logisticVE.groupA.multi*100-iv15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.casecount.groupA <- 1/length(iv15000$Results$casecountVE.groupA)*sum((iv15000$Results$casecountVE.groupA*100-iv15000$Results$trueVE.groupA[1]*100)^2)
mse.quad.multi.CoP.groupA
mse.quad.multi.typical.groupA 
mse.quad.multi.casecount.groupA

# MSE, overall
mse.quad.multi.CoP.all <- 1/length(iv15000$Results$logisticQuadCoPVE.all.multi)*sum((iv15000$Results$logisticQuadCoPVE.all.multi*100-iv15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.typical.all <- 1/length(iv15000$Results$logisticVE.all.multi)*sum((iv15000$Results$logisticVE.all.multi*100-iv15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.casecount.all <- 1/length(iv15000$Results$casecountVE.all)*sum((iv15000$Results$casecountVE.all*100-iv15000$Results$trueVE.all[1]*100)^2)
mse.quad.multi.CoP.all
mse.quad.multi.typical.all
mse.quad.multi.casecount.all

#### Graphical results, Figures: ------------------------------------------------------------------------

## Figure 1 ----
x <- seq(0,10, 0.01)
y <- seq(-10,10,0.02)

# PoDBAY
pmax <- 0.5
et50 <- 4
gama <- 8
PoD.podbay <- pmax/(1+(et50/x)^(-gama))

# Scaled logistic
lambda <- 0.5
beta <- 3
alpha <- 4
PoD.coude <- lambda/(1+exp(beta*(x-alpha)))

# Logistic
beta0 <- -2
beta1 <- -0.8
PoD.log <- 1/(1+exp(-(beta0+beta1*y)))

# Logistic quadratic
beta0 <- -4 # -4
beta1 <- 0.04 # 0.04
beta2 <- -0.04 # -0.04
PoD.log.quad <- 1/(1+exp(-(beta0+beta1*y + beta2*y^2)))

max <- -beta1/(2*beta2) # titer fot which PoD is maximal
p_max <- 1/(1+exp(-(beta0-beta1^2/4/beta2))) # maximal PoD

d <- data.frame(x = x,
                y = y,
                PoD.podbay = PoD.podbay,
                PoD.coude = PoD.coude,
                PoD.log = PoD.log,
                PoD.log.quad = PoD.log.quad
)

PoD.podbay <- ggplot(d, aes(x=x,y=PoD.podbay)) +
  geom_line() +
  labs(fill="") +
  xlim(-0.2,8) +
  xlab("Log Titer") +
  ylab("Probability of Disease") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())

PoD.coudeville <- ggplot(d, aes(x=x,y=PoD.coude)) +
  geom_line() +
  labs(fill="") +
  xlim(-0.2,8) +
  xlab("Log Titer") +
  ylab("Probability of Disease") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())

PoD.logistic <- ggplot(d) +
  geom_line(data = d %>% subset(y >= 0),
            aes(x=y,
                y=PoD.log)) +
  geom_line(data = d %>% subset(y < 0),
            aes(x=y,
                y=PoD.log),
            linetype = "dotted")+
  labs(fill="") +
  xlim(-8,8) +
  xlab("Log Titer") +
  ylab("Probability of Disease") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())

PoD.logistic.quadratic <- ggplot(d) + 
  geom_line(data = d %>% subset(y >= 0),
            aes(x=y,
                y=PoD.log.quad),
            size = 1.5) +
  geom_line(data = d %>% subset(y < 0),
            aes(x=y,
                y=PoD.log.quad),
            linetype = "dotted",
            size = 1.5)+
  labs(fill="") +
  xlim(-8,10) +
  xlab("Log Titer") +
  ylab("Probability of Disease") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())

PoD.logistic # A, logistic with linear term 
PoD.logistic.quadratic # B, logistic with quadratic term
PoD.coudeville # C, Dunning, 2006; Coudeville et al., 2010 
PoD.podbay # D, Hill function

## Figure 2 ----
titer <- seq(0,20,0.01)
vaccinated <- rnorm(10000,10,2)
control <- rnorm(5000,5,2)
val <- c(control,
         vaccinated)
name <- c(rep("control", length(control)),
          rep("vaccine", length(vaccinated)))
df <- data.frame(value = val,
                 variable = name) 
beta <- -0.33
beta.covariate <- 1
delta <- -2
y.PoD.logistic.A <- 1/(1+exp(-(beta*(beta.covariate^1)*titer + delta)))
ndata <- data.frame("titer" = titer,
                    "PoD" = y.PoD.logistic.A)
val <- c(control,
         vaccinated)
name <- c(rep("control", length(control)),
          rep("vaccine", length(vaccinated)))
df <- data.frame(value = val,
                 variable = name)
MSPoD.logistic.i <- ggplot(ndata, aes(x = titer)) +
  labs(x = 'Log2 SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(aes(y = PoD), color = "gray", size=2) + 
  ggtitle("True PoD curve: Logistic function, no age effect") +
  xlim(-0.65,15) +
  ylim(0,0.15)

beta.covariate <- 0.442
y.PoD.logistic.B.older <- 1/(1+exp(-(beta*(beta.covariate^1)*titer + delta)))
y.PoD.logistic.B.younger <- 1/(1+exp(-(beta*(beta.covariate^0)*titer + delta)))
age <- c(rep("older", length(y.PoD.logistic.B.older)),rep("younger", length(y.PoD.logistic.B.younger)))
ndata <- data.frame("titer" = titer,
                    "PoD" = c(y.PoD.logistic.B.older,y.PoD.logistic.B.younger),
                    "age" = age)
MSPoD.logistic.ii <- ggplot(ndata, aes(x = titer)) +
  labs(x = 'Log2 SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(aes(y = PoD, color = age), size=2) + 
  ggtitle("True PoD curve: Logistic function, effect size = 30%") +
  xlim(-0.65,15) +
  ylim(0,0.15) + 
  theme(legend.position="top")+
  scale_color_manual(values=c("darkslategray3", "gray20"))

pmax <- 0.033
et50 <- 7.12
gamma <- 7
et50.covariate <- 1
y.PoD.hill.A <- ifelse(titer > 0, pmax /( 1+(titer/(et50*(et50.covariate^1)))^(gamma)), pmax)
ndata <- data.frame("titer" = titer,
                    "PoD" = y.PoD.hill.A)
MSPoD.hill.iii <- ggplot(ndata, aes(x = titer)) +
  labs(x = 'Log2 SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(aes(y = PoD), color = "gray", size=2) + 
  ggtitle("True PoD curve: Hill function, no age effect") +
  xlim(-0.65,15) +
  ylim(0,0.15)

et50.covariate <- 1.361
y.PoD.hill.C.older <- ifelse(titer > 0, pmax /( 1+(titer/(et50*(et50.covariate^1)))^(gamma)), pmax)
y.PoD.hill.C.younger <- ifelse(titer > 0, pmax /( 1+(titer/(et50*(et50.covariate^0)))^(gamma)), pmax)
ndata <- data.frame("titer" = titer,
                    "PoD" = c(y.PoD.hill.C.older,y.PoD.hill.C.younger),
                    "age" = age)
MSPoD.hill.iv <- ggplot(ndata, aes(x = titer)) +
  labs(x = 'Log2 SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(aes(y = PoD, color = age), size=2) + 
  ggtitle("True PoD curve: Hill function, effect size = 30%") +
  xlim(-0.65,15) +
  ylim(0,0.15) + 
  theme(legend.position="top")+
  scale_color_manual(values=c("darkslategray3", "gray20"))

MSPoD.logistic.i # i 
MSPoD.logistic.ii # ii
MSPoD.hill.iii # iii
MSPoD.hill.iv # iv

## Figure 3 ----
MSridgeline <- ggplot(df, aes(x = value, y = variable)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "gray20",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(-0.1,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Immunogenicity') +
  ylab("Vaccination group") 

MSridgeline 

## Figure 4 ----
i15000$boxplots$selected # i 
ii15000$boxplots$selected # ii
iii15000$boxplots$selected # iii
iv15000$boxplots$selected # iv

# median CoP-based, older, scenario iv
median(iv15000$Results$selected.CoP.VE.groupA*100)
# true, older, scenario iv
iv15000$Results$trueVE.groupA[1]*100

## Figure 5 ----
i15000$widthplot$logistic.CoP.selected # i 
ii15000$widthplot$logistic.CoP.selected # ii
iii15000$widthplot$logistic.CoP.selected # iii
iv15000$widthplot$logistic.CoP.selected # iv

unlist(i15000$widthdatasets$logistic.CoP.selected$narrower)*100 # i, % of simulated trial the CoP-based method provides narrower CI (group A = older; group B = younger)
unlist(ii15000$widthdatasets$logistic.CoP.selected$narrower)*100 # ii, % of simulated trial the CoP-based method provides narrower CI (group A = older; group B = younger)
unlist(iii15000$widthdatasets$logistic.CoP.selected$narrower)*100 # iii, % of simulated trial the CoP-based method provides narrower CI (group A = older; group B = younger)
unlist(iv15000$widthdatasets$logistic.CoP.selected$narrower)*100 # iv, % of simulated trial the CoP-based method provides narrower CI (group A = older; group B = younger)

## Supplementary Figure S1 ----
idata <- i15000$Results
idata <- idata %>%
  mutate(wselected.CoP.VE.all = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                          selected.CoP.model == "linear multi" ~ logisticCoPVE.all.multi,
                                          selected.CoP.model == "quad add" ~ NA_real_,
                                          selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.all.multi,
                                          TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupA = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                             selected.CoP.model == "linear multi" ~ logisticCoPVE.groupA.multi,
                                             selected.CoP.model == "quad add" ~ NA_real_,
                                             selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.groupA.multi,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupB = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                             selected.CoP.model == "linear multi" ~ logisticCoPVE.groupB.multi,
                                             selected.CoP.model == "quad add" ~ NA_real_,
                                             selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.groupB.multi,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.typical.VE.all = ifelse(selected.typical.model == "add", NA_real_, logisticVE.all.multi)) %>%
  mutate(wselected.typical.VE.groupA = ifelse(selected.typical.model == "add", NA_real_, logisticVE.groupA.multi)) %>%
  mutate(wselected.typical.VE.groupB = ifelse(selected.typical.model == "add", NA_real_, logisticVE.groupB.multi)) %>%
  mutate(wselected.casecount.VE.groupA = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                                   selected.CoP.model == "linear multi" ~ casecountVE.groupA,
                                                   selected.CoP.model == "quad add" ~ NA_real_,
                                                   selected.CoP.model == "quad multi" ~ casecountVE.groupA,
                                                   TRUE ~ NA_real_))
i15000.wscatterplot <- ggplot(idata, aes(x=wselected.CoP.VE.groupA, y=wselected.casecount.VE.groupA-wselected.CoP.VE.groupA)) + 
  geom_point(color = "darkslategray3") +
  coord_cartesian(ylim = c(-0.5,0.5), xlim = c(0,1)) +
  geom_hline(aes(yintercept=0)) +
  themePlot 
efficacyDF<- wprepareDataVE(CoP.all = idata$wselected.CoP.VE.all,
                            typical.all = idata$wselected.typical.VE.all,
                            CoP.groupA = idata$wselected.CoP.VE.groupA,
                            typical.groupA = idata$wselected.typical.VE.groupA,
                            CoP.groupB = idata$wselected.CoP.VE.groupB,
                            typical.groupB = idata$wselected.typical.VE.groupB)
i15000.wviolinplot <- wggplotEfficacy(efficacyDF,
                                      idata$trueVE.all[1],
                                      idata$trueVE.groupA[1],
                                      idata$trueVE.groupB[1],
                                      "",
                                      "VE")

iidata <- ii15000$Results
iidata <- iidata %>%
  mutate(wselected.CoP.VE.all = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.all.add, 
                                          selected.CoP.model == "linear multi" ~ NA_real_,
                                          selected.CoP.model == "quad add" ~ logisticQuadCoPVE.all.add,
                                          selected.CoP.model == "quad multi" ~ NA_real_,
                                          TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupA = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.groupA.add, 
                                             selected.CoP.model == "linear multi" ~ NA_real_,
                                             selected.CoP.model == "quad add" ~ logisticQuadCoPVE.groupA.add,
                                             selected.CoP.model == "quad multi" ~ NA_real_,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupB = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.groupB.add, 
                                             selected.CoP.model == "linear multi" ~ NA_real_,
                                             selected.CoP.model == "quad add" ~ logisticQuadCoPVE.groupB.add,
                                             selected.CoP.model == "quad multi" ~ NA_real_,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.typical.VE.all = ifelse(selected.typical.model == "add", logisticVE.all.add, NA_real_)) %>%
  mutate(wselected.typical.VE.groupA = ifelse(selected.typical.model == "add", logisticVE.groupA.add, NA_real_)) %>%
  mutate(wselected.typical.VE.groupB = ifelse(selected.typical.model == "add", logisticVE.groupB.add, NA_real_)) %>%
  mutate(wselected.casecount.VE.groupA = case_when(selected.CoP.model == "linear add" ~ casecountVE.groupA, 
                                                   selected.CoP.model == "linear multi" ~ NA_real_,
                                                   selected.CoP.model == "quad add" ~ casecountVE.groupA,
                                                   selected.CoP.model == "quad multi" ~ NA_real_,
                                                   TRUE ~ NA_real_))
ii15000.wscatterplot <- ggplot(iidata, aes(x=wselected.CoP.VE.groupA, y=wselected.casecount.VE.groupA-wselected.CoP.VE.groupA)) + 
  geom_point(color = "darkslategray3") +
  coord_cartesian(ylim = c(-0.5,0.5), xlim = c(0,1)) +
  geom_hline(aes(yintercept=0)) +
  themePlot 
efficacyDF<- wprepareDataVE(CoP.all = iidata$wselected.CoP.VE.all,
                            typical.all = iidata$wselected.typical.VE.all,
                            CoP.groupA = iidata$wselected.CoP.VE.groupA,
                            typical.groupA = iidata$wselected.typical.VE.groupA,
                            CoP.groupB = iidata$wselected.CoP.VE.groupB,
                            typical.groupB = iidata$wselected.typical.VE.groupB)
ii15000.wviolinplot <- wggplotEfficacy(efficacyDF,
                                       iidata$trueVE.all[1],
                                       iidata$trueVE.groupA[1],
                                       iidata$trueVE.groupB[1],
                                       "",
                                       "VE")

iiidata <- iii15000$Results
iiidata <- iiidata %>%
  mutate(wselected.CoP.VE.all = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                          selected.CoP.model == "linear multi" ~ logisticCoPVE.all.multi,
                                          selected.CoP.model == "quad add" ~ NA_real_,
                                          selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.all.multi,
                                          TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupA = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                             selected.CoP.model == "linear multi" ~ logisticCoPVE.groupA.multi,
                                             selected.CoP.model == "quad add" ~ NA_real_,
                                             selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.groupA.multi,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupB = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                             selected.CoP.model == "linear multi" ~ logisticCoPVE.groupB.multi,
                                             selected.CoP.model == "quad add" ~ NA_real_,
                                             selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.groupB.multi,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.typical.VE.all = ifelse(selected.typical.model == "add", NA_real_, logisticVE.all.multi)) %>%
  mutate(wselected.typical.VE.groupA = ifelse(selected.typical.model == "add", NA_real_, logisticVE.groupA.multi)) %>%
  mutate(wselected.typical.VE.groupB = ifelse(selected.typical.model == "add", NA_real_, logisticVE.groupB.multi)) %>%
  mutate(wselected.casecount.VE.groupA = case_when(selected.CoP.model == "linear add" ~ NA_real_, 
                                                   selected.CoP.model == "linear multi" ~ casecountVE.groupA,
                                                   selected.CoP.model == "quad add" ~ NA_real_,
                                                   selected.CoP.model == "quad multi" ~ casecountVE.groupA,
                                                   TRUE ~ NA_real_))
iii15000.wscatterplot <- ggplot(iiidata, aes(x=wselected.CoP.VE.groupA, y=wselected.casecount.VE.groupA-wselected.CoP.VE.groupA)) + 
  geom_point(color = "darkslategray3") +
  coord_cartesian(ylim = c(-0.5,0.5), xlim = c(0,1)) +
  geom_hline(aes(yintercept=0)) +
  themePlot 
efficacyDF<- wprepareDataVE(CoP.all = iiidata$wselected.CoP.VE.all,
                            typical.all = iiidata$wselected.typical.VE.all,
                            CoP.groupA = iiidata$wselected.CoP.VE.groupA,
                            typical.groupA = iiidata$wselected.typical.VE.groupA,
                            CoP.groupB = iiidata$wselected.CoP.VE.groupB,
                            typical.groupB = iiidata$wselected.typical.VE.groupB)
iii15000.wviolinplot <- wggplotEfficacy(efficacyDF,
                                        iiidata$trueVE.all[1],
                                        iiidata$trueVE.groupA[1],
                                        iiidata$trueVE.groupB[1],
                                        "",
                                        "VE")

ivdata <- iv15000$Results
ivdata <- ivdata %>%
  mutate(wselected.CoP.VE.all = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.all.add, 
                                          selected.CoP.model == "linear multi" ~ NA_real_,
                                          selected.CoP.model == "quad add" ~ logisticQuadCoPVE.all.add,
                                          selected.CoP.model == "quad multi" ~ NA_real_,
                                          TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupA = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.groupA.add, 
                                             selected.CoP.model == "linear multi" ~ NA_real_,
                                             selected.CoP.model == "quad add" ~ logisticQuadCoPVE.groupA.add,
                                             selected.CoP.model == "quad multi" ~ NA_real_,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.CoP.VE.groupB = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.groupB.add, 
                                             selected.CoP.model == "linear multi" ~ NA_real_,
                                             selected.CoP.model == "quad add" ~ logisticQuadCoPVE.groupB.add,
                                             selected.CoP.model == "quad multi" ~ NA_real_,
                                             TRUE ~ NA_real_)) %>%
  mutate(wselected.typical.VE.all = ifelse(selected.typical.model == "add", logisticVE.all.add, NA_real_)) %>%
  mutate(wselected.typical.VE.groupA = ifelse(selected.typical.model == "add", logisticVE.groupA.add, NA_real_)) %>%
  mutate(wselected.typical.VE.groupB = ifelse(selected.typical.model == "add", logisticVE.groupB.add, NA_real_)) %>%
  mutate(wselected.casecount.VE.groupA = case_when(selected.CoP.model == "linear add" ~ casecountVE.groupA, 
                                                   selected.CoP.model == "linear multi" ~ NA_real_,
                                                   selected.CoP.model == "quad add" ~ casecountVE.groupA,
                                                   selected.CoP.model == "quad multi" ~ NA_real_,
                                                   TRUE ~ NA_real_))
iv15000.wscatterplot <-ggplot(ivdata, aes(x=wselected.CoP.VE.groupA, y=wselected.casecount.VE.groupA-wselected.CoP.VE.groupA)) + 
  geom_point(color = "darkslategray3") +
  coord_cartesian(ylim = c(-0.5,0.5), xlim = c(0,1)) +
  geom_hline(aes(yintercept=0)) +
  themePlot 
efficacyDF<- wprepareDataVE(CoP.all = ivdata$wselected.CoP.VE.all,
                            typical.all = ivdata$wselected.typical.VE.all,
                            CoP.groupA = ivdata$wselected.CoP.VE.groupA,
                            typical.groupA = ivdata$wselected.typical.VE.groupA,
                            CoP.groupB = ivdata$wselected.CoP.VE.groupB,
                            typical.groupB = ivdata$wselected.typical.VE.groupB)
iv15000.wviolinplot <- wggplotEfficacy(efficacyDF,
                                       ivdata$trueVE.all[1],
                                       ivdata$trueVE.groupA[1],
                                       ivdata$trueVE.groupB[1],
                                       "",
                                       "VE")

# A: Violin plots
i15000.wviolinplot # i 
ii15000.wviolinplot # ii
iii15000.wviolinplot # iii
iv15000.wviolinplot # iv

# B: Scatter plots
i15000.wscatterplot # i 
ii15000.wscatterplot # ii
iii15000.wscatterplot # iii
iv15000.wscatterplot # iv

## Supplementary Figure S2 ----
i15000$widthplot$casecount.logistic.selected # i 
ii15000$widthplot$casecount.logistic.selected # ii
iii15000$widthplot$casecount.logistic.selected # iii
iv15000$widthplot$casecount.logistic.selected # iv

unlist(i15000$narrowerCI$casecount.logistic.selected)*100 # i, % of simulated trial the typical logistic regression method provides narrower CI than case-counting (group A = older; group B = younger)
unlist(ii15000$narrowerCI$casecount.logistic.selected)*100 # ii, % of simulated trial the typical logistic regression method provides narrower CI than case-counting (group A = older; group B = younger)
unlist(iii15000$narrowerCI$casecount.logistic.selected)*100 # iii, % of simulated trial the typical logistic regression method provides narrower CI than case-counting (group A = older; group B = younger)
unlist(iv15000$narrowerCI$casecount.logistic.selected)*100 # iv, % of simulated trial the typical logistic regression method provides narrower CI than case-counting (group A = older; group B = younger)

## Supplementary Figure S3 ----
i15000$boxplots$linear.add # i 
ii15000$boxplots$linear.add # ii
iii15000$boxplots$linear.add # iii
iv15000$boxplots$linear.add # iv

## Supplementary Figure S4 --
i15000$boxplots$linear.multi # i 
ii15000$boxplots$linear.multi # ii
iii15000$boxplots$linear.multi # iii
iv15000$boxplots$linear.multi # iv

## Supplementary Figure S5 --
i15000$boxplots$quad.add # i 
ii15000$boxplots$quad.add # ii
iii15000$boxplots$quad.add # iii
iv15000$boxplots$quad.add # iv

## Supplementary Figure S6 --
i15000$boxplots$quad.multi # i 
ii15000$boxplots$quad.multi # ii
iii15000$boxplots$quad.multi # iii
iv15000$boxplots$quad.multi # iv

# median CoP-based, older, scenario ii
median(ii15000$Results$logisticQuadCoPVE.groupA.multi*100)
# true, older, scenario ii
ii15000$Results$trueVE.groupA[1]*100

# median CoP-based, older, scenario iv
median(iv15000$Results$logisticQuadCoPVE.groupA.multi*100)
# true, older, scenario iv
iv15000$Results$trueVE.groupA[1]*100

