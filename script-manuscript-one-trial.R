# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

## Prepare environment ----------------------------------------------------------------------------------
rm(list=ls())

# Package loading
library(Rcpp)
library(dplyr)
library(tidyr)
library(pracma)
library(modelr)
library(rlang)
library(tidyverse)
library(MASS)
library(vaxpmx) # version 0.0.3
library(ggplot2)
library(ggridges)
library(gridExtra)
library(ggdist)

# Code sourcing
source("functions-simulation-supplementary.R")
source("functions-process-results.R")
source("functions-define-truth-hill.R")

# Data generation, simulation ---------------------------------------------------------------------------
seed <- 10
set.seed(seed)
truth <- defineTruthExample()
vaccinated <- generatePopulation(truth$vaccinated$N,
                                 truth$vaccinated$mean,
                                 truth$vaccinated$stdDev)
control <- generatePopulation(truth$control$N,
                              truth$control$mean,
                              truth$control$stdDev)

# Generate covariate data of vaccinated and control populations based on the simulation inputs
vaccinated$addSubjectAttributes(
  covariateData(covariateInfo = truth$vaccinated$covariateInfo,
                Nrow = truth$vaccinated$N,
                vaccStatus = "vaccinated")
)
control$addSubjectAttributes(
  covariateData(covariateInfo = truth$control$covariateInfo,
                Nrow = truth$control$N,
                vaccStatus = "control")
)

# Assign probability of disease(PoD) to each subject based on their titer (and covariate data) - assignPoD()
vaccinated$assignPoD( 
  PoD(titers = vaccinated$titers,
      attributes = vaccinated$subjectAttributes,
      PoDModel = truth$PoDModel,
      valueParams = truth$PoDModel$valueParams,
      adjustTiters = truth$adjustTiters,
      adjustFrom = truth$adjustFrom,
      adjustTo = truth$adjustTo)
)
control$assignPoD(
  PoD(titer = control$titers,
      attributes = control$subjectAttributes,
      PoDModel = truth$PoDModel,
      valueParams = truth$PoDModel$valueParams,
      adjustTiters = truth$adjustTiters,
      adjustFrom = truth$adjustFrom,
      adjustTo = truth$adjustTo)
)

# Calculate the true efficacy - using true PoD curve and true distributions
trueEfficacyAll <- EfficacyComputation_MV(vaccPOP = vaccinated, # this has to be a population
                                          contPOP = control, # this has to be a population
                                          covInfo = truth$vaccinated$covariateInfo,
                                          paramIn = truth$PoDModel$valueParams,
                                          podModelObject = truth$PoDModel)
covInfo <- truth$vaccinated$covariateInfo
covInfoGroupA <- covInfo 
covInfoGroupA$var.covariateBinary$prob <- c(1,0)
covInfoGroupB <- covInfo 
covInfoGroupB$var.covariateBinary$prob <- c(0,1)
trueEfficacyGroupA <- EfficacyComputation_MV(vaccPOP = vaccinated, 
                                             contPOP = control, 
                                             covInfo = covInfoGroupA,
                                             paramIn = truth$PoDModel$valueParams,
                                             podModelObject = truth$PoDModel)
trueEfficacyGroupB <- EfficacyComputation_MV(vaccPOP = vaccinated, 
                                             contPOP = control, 
                                             covInfo = covInfoGroupB,
                                             paramIn = truth$PoDModel$valueParams,
                                             podModelObject = truth$PoDModel)

## Simulate clinical trial ------------------------------------------------------------------------------
# All
simulatedTest <- ClinicalTrial(vaccinated, control)

# Group A, Older
vaccinatedGroupA <- getPopulationGroup(population1 = vaccinated,
                                       attribute = "var.covariateBinary",
                                       group = "valueA")
controlGroupA <- getPopulationGroup(population1 = control,
                                    attribute = "var.covariateBinary",
                                    group = "valueA")
casecountEfficacyGroupA <- 1 - (vaccinatedGroupA$getDiseasedCount() / vaccinatedGroupA$N) /
  (controlGroupA$getDiseasedCount() / controlGroupA$N)
confidenceIntervalGroupA <- waldCI(vaccinatedGroupA, controlGroupA, 0.95)

# Group B, Younger
vaccinatedGroupB <- getPopulationGroup(population1 = vaccinated,
                                       attribute = "var.covariateBinary",
                                       group = "valueB")
controlGroupB <- getPopulationGroup(population1 = control,
                                    attribute = "var.covariateBinary",
                                    group = "valueB")
casecountEfficacyGroupB <- 1 - (vaccinatedGroupB$getDiseasedCount() / vaccinatedGroupB$N) /
  (controlGroupB$getDiseasedCount() / controlGroupB$N)
confidenceIntervalGroupB <- waldCI(vaccinatedGroupB, controlGroupB, 0.95)

## Create immunogenicity sample 
diseasedAll <- ExtractDiseased(vaccinated = simulatedTest$vaccinated, 
                               control = simulatedTest$control)
diseasedGroupA <- getPopulationGroup(population1 = diseasedAll,
                                     attribute = "var.covariateBinary",
                                     group = "valueA")
diseasedGroupB <- getPopulationGroup(population1 = diseasedAll,
                                     attribute = "var.covariateBinary",
                                     group = "valueB")
nondiseasedAll <- ExtractNondiseased(vaccinated = simulatedTest$vaccinated, 
                                     control = simulatedTest$control)
nondiseasedGroupA <- getPopulationGroup(population1 = nondiseasedAll,
                                        attribute = "var.covariateBinary",
                                        group = "valueA")
nondiseasedGroupB <- getPopulationGroup(population1 = nondiseasedAll,
                                        attribute = "var.covariateBinary",
                                        group = "valueB")
immunogenicityAll <- BlindSampling(diseased = diseasedAll,
                                   nondiseased = nondiseasedAll,
                                   method = truth$method)
immunogenicityGroupA <- list()
immunogenicityGroupA$immunogenicityNondiseased <- getPopulationGroup(population1 = immunogenicityAll$immunogenicityNondiseased,
                                                                     attribute = "var.covariateBinary",
                                                                     group = "valueA")
immunogenicityGroupA$immunogenicityVaccinated <- getPopulationGroup(population1 = immunogenicityAll$immunogenicityVaccinated,
                                                                    attribute = "var.covariateBinary",
                                                                    group = "valueA")
immunogenicityGroupA$immunogenicityControl <- getPopulationGroup(population1 = immunogenicityAll$immunogenicityControl,
                                                                 attribute = "var.covariateBinary",
                                                                 group = "valueA")
immunogenicityGroupB <- list()
immunogenicityGroupB$immunogenicityNondiseased <- getPopulationGroup(population1 = immunogenicityAll$immunogenicityNondiseased,
                                                                     attribute = "var.covariateBinary",
                                                                     group = "valueB")
immunogenicityGroupB$immunogenicityVaccinated <- getPopulationGroup(population1 = immunogenicityAll$immunogenicityVaccinated,
                                                                    attribute = "var.covariateBinary",
                                                                    group = "valueB")
immunogenicityGroupB$immunogenicityControl <- getPopulationGroup(population1 = immunogenicityAll$immunogenicityControl,
                                                                 attribute = "var.covariateBinary",
                                                                 group = "valueB")

## Clinical trial results -------------------------------------------------------------------------------
NDiseased <- list(
  all = diseasedAll$N,
  groupA = diseasedGroupA$N,
  groupB = diseasedGroupB$N
)

NNondiseased <- list(
  all = nondiseasedAll$N,
  groupA = nondiseasedGroupA$N,
  groupB = nondiseasedGroupB$N
)

NImmunogenicityNondiseased <- list(
  all = immunogenicityAll$immunogenicityNondiseased$N,
  groupA = immunogenicityGroupA$immunogenicityNondiseased$N,
  groupB = immunogenicityGroupB$immunogenicityNondiseased$N
)

trueEfficacy <- list(
  all = trueEfficacyAll,
  groupA = trueEfficacyGroupA,
  groupB = trueEfficacyGroupB
)

resultTrial <- list(
  seed = seed,
  NDiseased = NDiseased,
  NNondiseased = NNondiseased,
  NImmunogenicityNondiseased = NImmunogenicityNondiseased,
  trueEfficacy = trueEfficacy
)

VE <- list(
  all = simulatedTest$efficacy,
  groupA = casecountEfficacyGroupA,
  groupB = casecountEfficacyGroupB
)

VECI <- list(
  all = simulatedTest$confidenceInterval,
  groupA = confidenceIntervalGroupA,
  groupB = confidenceIntervalGroupB
)

resultCaseCount <- list(
  VE = VE,
  VECI = VECI
)

## Logistic regression ----------------------------------------------------------------------------------
data <- populationToDF(vaccinated, control)
data <- data %>%
  mutate(var.quadtiter = var.titer*var.titer) %>%
  mutate(vaccine = as.numeric(var.vaccStatus=="vaccinated"))

data.control <- data %>% filter(var.vaccStatus == "control")
data.vaccinated <- data %>% filter(var.vaccStatus == "vaccinated")
data.groupA <- data %>% filter(var.covariateBinary == "valueA")
data.groupB <- data %>% filter(var.covariateBinary == "valueB")
data.groupA.vaccinated <- data.groupA %>% filter(var.vaccStatus == "vaccinated")
data.groupA.control <- data.groupA %>% filter(var.vaccStatus == "control")
data.groupB.vaccinated <- data.groupB %>% filter(var.vaccStatus == "vaccinated")
data.groupB.control <- data.groupB %>% filter(var.vaccStatus == "control")

## Titer not used in the model (Typical approach)
# Without interaction
logisticFit.add <- glm(y.diseaseStatus ~ var.vaccStatus + var.covariateBinary, data = data, family = binomial()) 
summary(logisticFit.add)
p.logisticFit.add.vaccStatuscontrol <- coef(summary(logisticFit.add))[[2,4]]
p.logisticFit.add.covariateBinaryvalueB <- coef(summary(logisticFit.add))[[3,4]]
AIC.logistic.add <- logisticFit.add$aic
# Overall VE
VE.logisticFit.add <- ve(logisticFit.add, data, nboot = 2000)
# Group A, Older
VE.logisticFit.add.groupA <- ve(logisticFit.add, data.groupA, nboot = 2000)
# Group B, Younger
VE.logisticFit.add.groupB <- ve(logisticFit.add, data.groupB, nboot = 2000)

# With interaction
logisticFit.multi <- glm(y.diseaseStatus ~ var.vaccStatus * var.covariateBinary, data = data, family = binomial()) 
summary(logisticFit.multi)
p.logisticFit.multi.vaccStatuscontrol <- coef(summary(logisticFit.multi))[[2,4]]
p.logisticFit.multi.covariateBinaryvalueB <- coef(summary(logisticFit.multi))[[3,4]]
p.logisticFit.multi.interaction <- coef(summary(logisticFit.multi))[[4,4]]
AIC.logistic.multi <- logisticFit.multi$aic
# Overall VE
VE.logisticFit.multi <- ve(logisticFit.multi, data, nboot = 2000)
# Group A, Older
VE.logisticFit.multi.groupA <- ve(logisticFit.multi, data.groupA, nboot = 2000)
# Group B, Younger
VE.logisticFit.multi.groupB <- ve(logisticFit.multi, data.groupB, nboot = 2000)

## Titer used in the model (CoP-based approach) 
# Linear term for titer, without interaction
logisticFit.CoP.add <- glm(y.diseaseStatus ~ var.titer + var.covariateBinary, data = data, family = binomial()) 
summary(logisticFit.CoP.add)
p.logisticFit.CoP.add.titer <- coef(summary(logisticFit.CoP.add))[[2,4]]
p.logisticFit.CoP.add.covariateBinaryvalueB <- coef(summary(logisticFit.CoP.add))[[3,4]]
AIC.logistic.CoP.add <- logisticFit.CoP.add$aic
# Overall VE
VE.logisticFit.CoP.add <- ve(logisticFit.CoP.add, data, nboot = 2000)
# Group A, Older
VE.logisticFit.CoP.add.groupA <- ve(logisticFit.CoP.add, data.groupA, nboot = 2000)
# Group B, Younger
VE.logisticFit.CoP.add.groupB <- ve(logisticFit.CoP.add, data.groupB, nboot = 2000)

# Linear term for titer, with interaction
logisticFit.CoP.multi <- glm(y.diseaseStatus ~ var.titer*var.covariateBinary, data = data, family = binomial()) 
summary(logisticFit.CoP.multi)
p.logisticFit.CoP.multi.titer <- coef(summary(logisticFit.CoP.multi))[[2,4]]
p.logisticFit.CoP.multi.covariateBinaryvalueB <- coef(summary(logisticFit.CoP.multi))[[3,4]]
p.logisticFit.CoP.multi.interaction <- coef(summary(logisticFit.CoP.multi))[[4,4]]
AIC.logistic.CoP.multi <- logisticFit.CoP.multi$aic
# Overall VE
VE.logisticFit.CoP.multi <- ve(logisticFit.CoP.multi, data, nboot = 2000)
# Group A, Older
VE.logisticFit.CoP.multi.groupA <- ve(logisticFit.CoP.multi, data.groupA, nboot = 2000)
# Group B, Younger
VE.logisticFit.CoP.multi.groupB <- ve(logisticFit.CoP.multi, data.groupB, nboot = 2000)

# Quadratic term for titer, without interaction
logisticFit.quadCoP.add <- glm(y.diseaseStatus ~ var.titer + var.quadtiter + var.covariateBinary, data = data, family = binomial()) 
summary(logisticFit.quadCoP.add)
p.logisticFit.quadCoP.add.titer <- coef(summary(logisticFit.quadCoP.add))[[2,4]]
p.logisticFit.quadCoP.add.quadtiter <- coef(summary(logisticFit.quadCoP.add))[[3,4]]
p.logisticFit.quadCoP.add.covariateBinaryvalueB <- coef(summary(logisticFit.quadCoP.add))[[4,4]]
AIC.logistic.quadCoP.add <- logisticFit.quadCoP.add$aic
# Overall VE
VE.logisticFit.quadCoP.add <- ve(logisticFit.quadCoP.add, data, nboot = 2000)
# Group A, Older
VE.logisticFit.quadCoP.add.groupA <- ve(logisticFit.quadCoP.add, data.groupA, nboot = 2000)
# Group B, Younger
VE.logisticFit.quadCoP.add.groupB <- ve(logisticFit.quadCoP.add, data.groupB, nboot = 2000)

# Quadratic term for titer, with interaction
logisticFit.quadCoP.multi <- glm(y.diseaseStatus ~ var.titer + var.quadtiter*var.covariateBinary, data = data, family = binomial()) 
summary(logisticFit.quadCoP.multi)
p.logisticFit.quadCoP.multi.titer <- coef(summary(logisticFit.quadCoP.multi))[[2,4]]
p.logisticFit.quadCoP.multi.quadtiter <- coef(summary(logisticFit.quadCoP.multi))[[3,4]]
p.logisticFit.quadCoP.multi.covariateBinaryvalueB <- coef(summary(logisticFit.quadCoP.multi))[[4,4]]
p.logisticFit.quadCoP.multi.interaction <- coef(summary(logisticFit.quadCoP.multi))[[5,4]]
AIC.logistic.quadCoP.multi <- logisticFit.quadCoP.multi$aic
# Overall VE
VE.logisticFit.quadCoP.multi <- ve(logisticFit.quadCoP.multi, data, nboot = 2000)
# Group A, Older
VE.logisticFit.quadCoP.multi.groupA <- ve(logisticFit.quadCoP.multi, data.groupA, nboot = 2000)
# Group B, Younger
VE.logisticFit.quadCoP.multi.groupB <- ve(logisticFit.quadCoP.multi, data.groupB, nboot = 2000)

## Logistic regression results --------------------------------------------------------------------------

# Logistic regression, typical approach
pval <- list(
  vaccStatus.add = p.logisticFit.add.vaccStatuscontrol,
  covariateBinary.add = p.logisticFit.add.covariateBinaryvalueB,
  vaccStatus.multi = p.logisticFit.multi.vaccStatuscontrol,
  covariateBinary.multi = p.logisticFit.multi.covariateBinaryvalueB,
  interaction.multi = p.logisticFit.multi.interaction
)
aic <- list(
  aic.add = AIC.logistic.add,
  aic.multi = AIC.logistic.multi
)
VE <- list(all.add = VE.logisticFit.add[[1]],
           groupA.add = VE.logisticFit.add.groupA[[1]],
           groupB.add =  VE.logisticFit.add.groupB[[1]],
           all.multi = VE.logisticFit.multi[[1]],
           groupA.multi = VE.logisticFit.multi.groupA[[1]],
           groupB.multi =  VE.logisticFit.multi.groupB[[1]])
VECI <- list(all.add = VE.logisticFit.add$CI,
             groupA.add = VE.logisticFit.add.groupA$CI,
             groupB.add =  VE.logisticFit.add.groupB$CI,
             all.multi = VE.logisticFit.multi$CI,
             groupA.multi = VE.logisticFit.multi.groupA$CI,
             groupB.multi =  VE.logisticFit.multi.groupB$CI)
resultLogistic <- list(
  pval = pval,
  aic = aic,
  VE = VE,
  VECI = VECI
)

# Logistic regression, CoP-based approach, linear titer
pval <- list(
  titer.add = p.logisticFit.CoP.add.titer, 
  covariateBinary.add = p.logisticFit.CoP.add.covariateBinaryvalueB,
  titer.multi = p.logisticFit.CoP.multi.titer, 
  covariateBinary.multi = p.logisticFit.CoP.multi.covariateBinaryvalueB,
  interaction.multi = p.logisticFit.CoP.multi.interaction 
)
aic <- list(
  aic.add = AIC.logistic.CoP.add,
  aic.multi = AIC.logistic.CoP.multi
)
VE <- list(
  all.add = VE.logisticFit.CoP.add[[1]],
  groupA.add = VE.logisticFit.CoP.add.groupA[[1]],
  groupB.add = VE.logisticFit.CoP.add.groupB[[1]],
  all.multi = VE.logisticFit.CoP.multi[[1]],
  groupA.multi = VE.logisticFit.CoP.multi.groupA[[1]],
  groupB.multi = VE.logisticFit.CoP.multi.groupB[[1]]
)
VECI <- list(all.add = VE.logisticFit.CoP.add$CI,
             groupA.add = VE.logisticFit.CoP.add.groupA$CI,
             groupB.add =  VE.logisticFit.CoP.add.groupB$CI,
             all.multi = VE.logisticFit.CoP.multi$CI,
             groupA.multi = VE.logisticFit.CoP.multi.groupA$CI,
             groupB.multi =  VE.logisticFit.CoP.multi.groupB$CI
)
resultLogisticCoP <- list(
  VE = VE,
  VECI = VECI,
  pval = pval,
  aic = aic
)

# Logistic regression, CoP-based approach, quadratic titer
pval <- list(
  titer.add = p.logisticFit.quadCoP.add.titer, 
  quadtiter.add = p.logisticFit.quadCoP.add.quadtiter,
  covariateBinary.add = p.logisticFit.quadCoP.add.covariateBinaryvalueB,
  titer.multi = p.logisticFit.quadCoP.multi.titer, 
  quadtiter.multi = p.logisticFit.quadCoP.multi.quadtiter, 
  covariateBinary.multi = p.logisticFit.quadCoP.multi.covariateBinaryvalueB,
  interaction.multi = p.logisticFit.quadCoP.multi.interaction 
)
aic <- list(
  aic.add = AIC.logistic.quadCoP.add,
  aic.multi = AIC.logistic.quadCoP.multi
)
VE <- list(
  all.add = VE.logisticFit.quadCoP.add[[1]],
  groupA.add = VE.logisticFit.quadCoP.add.groupA[[1]],
  groupB.add = VE.logisticFit.quadCoP.add.groupB[[1]],
  all.multi = VE.logisticFit.quadCoP.multi[[1]],
  groupA.multi = VE.logisticFit.quadCoP.multi.groupA[[1]],
  groupB.multi = VE.logisticFit.quadCoP.multi.groupB[[1]]
)
VECI <- list(all.add = VE.logisticFit.quadCoP.add$CI,
             groupA.add = VE.logisticFit.quadCoP.add.groupA$CI,
             groupB.add =  VE.logisticFit.quadCoP.add.groupB$CI,
             all.multi = VE.logisticFit.quadCoP.multi$CI,
             groupA.multi = VE.logisticFit.quadCoP.multi.groupA$CI,
             groupB.multi =  VE.logisticFit.quadCoP.multi.groupB$CI
)
resultLogisticQuadCoP <- list(
  VE = VE,
  VECI = VECI,
  pval = pval,
  aic = aic
)

## Compile all results ----------------------------------------------------------------------------------
result <- list(
  resultTrial = resultTrial, 
  resultCaseCount = resultCaseCount, 
  resultLogistic = resultLogistic,
  resultLogisticCoP = resultLogisticCoP,
  resultLogisticQuadCoP = resultLogisticQuadCoP
)

## Numerical results, Tables: ---------------------------------------------------------------------------

## Table 4 ----
summary(logisticFit.add) # model 15
summary(logisticFit.multi) # model 16
summary(logisticFit.CoP.add) # model 17
summary(logisticFit.CoP.multi) # model 18
summary(logisticFit.quadCoP.add) # model 19
summary(logisticFit.quadCoP.multi) # model 20

## Table 5 ----

## Younger

# Control
nrow(data.groupB.control) # number of subjects in control group, also number of subjects in immunogenicity subset of control group
data.groupB.control.cases <- data.groupB.control %>% filter(y.diseaseStatus == 1)
nrow(data.groupB.control.cases) # number of cases in control group
GMT.groupB.control <- 2^(mean(controlGroupB$titers))
ConvertNormalToGeometric(mean(controlGroupB$titers),
                         sd(controlGroupB$titers),
                         controlGroupB$N) # immunogenicity in control group

# Vaccinated
nrow(data.groupB.vaccinated) # number of subjects in vaccinated group, also number of subjects in immunogenicity subset of vaccinated group
data.groupB.vaccinated.cases <- data.groupB.vaccinated %>% filter(y.diseaseStatus == 1)
nrow(data.groupB.vaccinated.cases) # number of cases in vaccinated group
GMT.groupB.vaccinated <- 2^(mean(vaccinatedGroupB$titers))
ConvertNormalToGeometric(mean(vaccinatedGroupB$titers),
                         sd(vaccinatedGroupB$titers),
                         vaccinatedGroupB$N) # immunogenicity in vaccinated group

# Vaccine efficacy, %
resultCaseCount$VE$groupB*100
unlist(resultCaseCount$VECI$groupB)*100 # case-count

resultLogistic$VE$groupB.multi
resultLogistic$VECI$groupB.multi # logistic, typical

resultLogisticQuadCoP$VE$groupB.multi
resultLogisticQuadCoP$VECI$groupB.multi # logistic, immunogenicity-based

## Older

# Control
nrow(data.groupA.control) # number of subjects in control group, also number of subjects in immunogenicity subset of control group
data.groupA.control.cases <- data.groupA.control %>% filter(y.diseaseStatus == 1)
nrow(data.groupA.control.cases) # number of cases in control group
GMT.groupA.control <- 2^(mean(controlGroupA$titers))
ConvertNormalToGeometric(mean(controlGroupA$titers),
                         sd(controlGroupA$titers),
                         controlGroupA$N) # immunogenicity in control group

# Vaccinated
nrow(data.groupA.vaccinated) # number of subjects in vaccinated group, also number of subjects in immunogenicity subset of vaccinated group
data.groupA.vaccinated.cases <- data.groupA.vaccinated %>% filter(y.diseaseStatus == 1)
nrow(data.groupA.vaccinated.cases) # number of cases in vaccinated group
GMT.groupB.vaccinated <- 2^(mean(vaccinatedGroupA$titers))
ConvertNormalToGeometric(mean(vaccinatedGroupA$titers),
                         sd(vaccinatedGroupA$titers),
                         vaccinatedGroupA$N) # immunogenicity in vaccinated group

# Vaccine efficacy, %
resultCaseCount$VE$groupA*100
unlist(resultCaseCount$VECI$groupA)*100 # case-count

resultLogistic$VE$groupA.multi
resultLogistic$VECI$groupA.multi # logistic, typical

resultLogisticQuadCoP$VE$groupA.multi
resultLogisticQuadCoP$VECI$groupA.multi # logistic, immunogenicity-based

## Graphical results, Figures: ------------------------------------------------------------------------

## Figure 6 ----
titers <- seq(0,15,0.01)
quadtiters <- titers^2
age.older <- rep("valueA",length(titers))
age.younger <- rep("valueB",length(titers))
df.older <- data.frame(var.titer = titers,
                       var.quadtiter = quadtiters,
                       var.covariateBinary = age.older)
df.younger <- data.frame(var.titer = titers,
                         var.quadtiter = quadtiters,
                         var.covariateBinary = age.younger)
df.older$var.covariateBinary <- factor(df.older$var.covariateBinary ,
                                     levels = c("valueA","valueB"),ordered = TRUE)
df.younger$var.covariateBinary <- factor(df.younger$var.covariateBinary ,
                                     levels = c("valueA","valueB"),ordered = TRUE)

ilink <- family(logisticFit.quadCoP.multi)$linkinv
ndata <- bind_cols(df.older, setNames(as_tibble(predict(logisticFit.quadCoP.multi, df.older, se.fit = TRUE)[1:2]),
                                     c('fit.link.older','se.link.older')))
ndata <- bind_cols(ndata, setNames(as_tibble(predict(logisticFit.quadCoP.multi, df.younger, se.fit = TRUE)[1:2]),
                                   c('fit.link.younger','se.link.younger')))
ndata <- mutate(ndata,
                PoD.logistic.older  = ilink(fit.link.older),
                CI.upr.logistic.older = ilink(fit.link.older + (1.96 * se.link.older)),
                CI.lwr.logistic.older = ilink(fit.link.older - (1.96 * se.link.older)),
                PoD.logistic.younger  = ilink(fit.link.younger),
                CI.upr.logistic.younger = ilink(fit.link.younger + (1.96 * se.link.younger)),
                CI.lwr.logistic.younger = ilink(fit.link.younger - (1.96 * se.link.younger)))
PoD.logistic.older <- ggplot(ndata, aes(x = var.titer)) +
  labs(x = '', y = '') +
  theme_bw() +
  geom_line(aes(y = PoD.logistic.older), color = "darkslategray3", size=2) + 
  geom_ribbon(data = ndata,
              aes(ymin = CI.lwr.logistic.older, ymax = CI.upr.logistic.older),
              alpha = 0.3,
              fill = "darkslategray3")+
  ggtitle("Older") +
  xlim(-0.1,15) +
  ylim(0,0.15) +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank()  #remove y axis labels
  )
MSridgeline.older <- ggplot(data.groupA, aes(x = var.titer, y = factor(var.vaccStatus))) +
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
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank()  #remove y axis labels
  )
PoD.logistic.younger <- ggplot(ndata, aes(x = var.titer)) +
  #  geom_line(aes(y = PoD.logistic)) +
  labs(x = '', y = '') +
  theme_bw() +
  geom_line(aes(y = PoD.logistic.younger), color = "gray20", size = 2) + 
  geom_ribbon(data = ndata,
              aes(ymin = CI.lwr.logistic.younger, ymax = CI.upr.logistic.younger),
              alpha = 0.3,
              fill = "gray20")+
  ggtitle("Younger") +
  xlim(-0.1,15) +
  ylim(0,0.15) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank()  #remove y axis labels
  )
MSridgeline.younger <-ggplot(data.groupB, aes(x = var.titer, y = factor(var.vaccStatus))) +
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
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.text.y=element_blank()  #remove y axis labels
  )

grid.arrange(PoD.logistic.older, MSridgeline.older, nrow = 2) # older 
grid.arrange(PoD.logistic.younger, MSridgeline.younger, nrow = 2) # younger

## Revision 7th February, 2024 -------------------------------------------------------------

## In the typical approach, fit a probit model to estimate VE
## Titer not used in the model 

set.seed(seed)
# Without interaction
probitFit.add <- glm(y.diseaseStatus ~ var.vaccStatus + var.covariateBinary, data = data, family = binomial(link=probit)) 
summary(probitFit.add)
p.probitFit.add.vaccStatuscontrol <- coef(summary(probitFit.add))[[2,4]]
p.probitFit.add.covariateBinaryvalueB <- coef(summary(probitFit.add))[[3,4]]
AIC.probit.add <- probitFit.add$aic
# Overall VE
VE.probitFit.add <- ve(probitFit.add, data, nboot = 2000)
# Group A, Older
VE.probitFit.add.groupA <- ve(probitFit.add, data.groupA, nboot = 2000)
# Group B, Younger
VE.probitFit.add.groupB <- ve(probitFit.add, data.groupB, nboot = 2000)

# With interaction
probitFit.multi <- glm(y.diseaseStatus ~ var.vaccStatus * var.covariateBinary, data = data, family = binomial(link=probit)) 
summary(probitFit.multi)
p.probitFit.multi.vaccStatuscontrol <- coef(summary(probitFit.multi))[[2,4]]
p.probitFit.multi.covariateBinaryvalueB <- coef(summary(probitFit.multi))[[3,4]]
p.probitFit.multi.interaction <- coef(summary(probitFit.multi))[[4,4]]
AIC.probit.multi <- probitFit.multi$aic
# Overall VE
VE.probitFit.multi <- ve(probitFit.multi, data, nboot = 2000)
# Group A, Older
VE.probitFit.multi.groupA <- ve(probitFit.multi, data.groupA, nboot = 2000)
# Group B, Younger
VE.probitFit.multi.groupB <- ve(probitFit.multi, data.groupB, nboot = 2000)
