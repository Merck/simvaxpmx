# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

Simulation <- function(seed = 1) {
  
  ## Loading trial inputs -----------------------------------------------------
  set.seed(seed)
  
  # Load the truth for the simulation (from the functions-define-truth.R file)
  truth <- defineTruth()

  ## Generate Populations -------------------------------------------------------
  # Generate titer data of vaccinated and control populations based on the simulation inputs
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

  ## Simulate clinical trial ---------------------------------------------------
  
  # Assign disease status (DS) based on the PoD of each subject 
  # Estimate case-count efficacy and confidence intervals 
  
  # All
  simulatedTest <- ClinicalTrial(vaccinated, control)
  
  # Group A
  vaccinatedGroupA <- getPopulationGroup(population1 = vaccinated,
                                         attribute = "var.covariateBinary",
                                         group = "valueA")
  controlGroupA <- getPopulationGroup(population1 = control,
                                      attribute = "var.covariateBinary",
                                      group = "valueA")
  casecountEfficacyGroupA <- 1 - (vaccinatedGroupA$getDiseasedCount() / vaccinatedGroupA$N) /
    (controlGroupA$getDiseasedCount() / controlGroupA$N)
  confidenceIntervalGroupA <- waldCI(vaccinatedGroupA, controlGroupA, 0.95)
  
  # Group B  
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
  
  # Titers and attributes of the diseased population are pooled together 
  # all diseased subjects (disease status = T) are extracted
  diseasedAll <- ExtractDiseased(vaccinated = simulatedTest$vaccinated, 
                                 control = simulatedTest$control)
  
  diseasedGroupA <- getPopulationGroup(population1 = diseasedAll,
                                       attribute = "var.covariateBinary",
                                       group = "valueA")
  
  diseasedGroupB <- getPopulationGroup(population1 = diseasedAll,
                                       attribute = "var.covariateBinary",
                                       group = "valueB")
  
  # Titers and atributes of the non-diseased population are pooled together 
  # all non-diseased subjects (disease status = F) are extracted
  nondiseasedAll <- ExtractNondiseased(vaccinated = simulatedTest$vaccinated, 
                                       control = simulatedTest$control)
  
  nondiseasedGroupA <- getPopulationGroup(population1 = nondiseasedAll,
                                          attribute = "var.covariateBinary",
                                          group = "valueA")
  
  nondiseasedGroupB <- getPopulationGroup(population1 = nondiseasedAll,
                                          attribute = "var.covariateBinary",
                                          group = "valueB")
  
  # Immunogenicity subset is created by sampling from the whole non-diseased population 
  # Either "Full", or "Fixed" blindsampling methods, as defined in the functions-define-truth.R file
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
  
  ## Clinical trial results -----------------------------------------------------
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
  
  ## Logistic regression ------------------------------------------------------

  #data <- populationToDF(immunogenicityAll$immunogenicityNondiseased, diseasedAll) # for immunogenicity subset simulations
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

  
  ## Titer not taken into account (Typical approach)
  logisticFit.add <- glm(y.diseaseStatus ~ var.vaccStatus + var.covariateBinary, data = data, family = binomial()) 
  summary(logisticFit.add)
  p.logisticFit.add.vaccStatuscontrol <- coef(summary(logisticFit.add))[[2,4]]
  p.logisticFit.add.covariateBinaryvalueB <- coef(summary(logisticFit.add))[[3,4]]
  AIC.logistic.add <- logisticFit.add$aic
  # Overall VE
  VE.logisticFit.add <- ve(logisticFit.add, data, nboot = 2000)
  # Group A
  VE.logisticFit.add.groupA <- ve(logisticFit.add, data.groupA, nboot = 2000)
  # Group B
  VE.logisticFit.add.groupB <- ve(logisticFit.add, data.groupB, nboot = 2000)
  
  logisticFit.multi <- glm(y.diseaseStatus ~ var.vaccStatus * var.covariateBinary, data = data, family = binomial()) 
  summary(logisticFit.multi)
  p.logisticFit.multi.vaccStatuscontrol <- coef(summary(logisticFit.multi))[[2,4]]
  p.logisticFit.multi.covariateBinaryvalueB <- coef(summary(logisticFit.multi))[[3,4]]
  p.logisticFit.multi.interaction <- coef(summary(logisticFit.multi))[[4,4]]
  AIC.logistic.multi <- logisticFit.multi$aic
  # Overall VE
  VE.logisticFit.multi <- ve(logisticFit.multi, data, nboot = 2000)
  # Group A
  VE.logisticFit.multi.groupA <- ve(logisticFit.multi, data.groupA, nboot = 2000)
  # Group B
  VE.logisticFit.multi.groupB <- ve(logisticFit.multi, data.groupB, nboot = 2000)
  
  ## Titer taken into account (CoP-based approach) 
  # Linear term for titer
  logisticFit.CoP.add <- glm(y.diseaseStatus ~ var.titer + var.covariateBinary, data = data, family = binomial()) 
  summary(logisticFit.CoP.add)
  p.logisticFit.CoP.add.titer <- coef(summary(logisticFit.CoP.add))[[2,4]]
  p.logisticFit.CoP.add.covariateBinaryvalueB <- coef(summary(logisticFit.CoP.add))[[3,4]]
  AIC.logistic.CoP.add <- logisticFit.CoP.add$aic
  # Overall VE
  VE.logisticFit.CoP.add <- ve(logisticFit.CoP.add, data, nboot = 2000)
  # Group A
  VE.logisticFit.CoP.add.groupA <- ve(logisticFit.CoP.add, data.groupA, nboot = 2000)
  # Group B
  VE.logisticFit.CoP.add.groupB <- ve(logisticFit.CoP.add, data.groupB, nboot = 2000)
  
  logisticFit.CoP.multi <- glm(y.diseaseStatus ~ var.titer*var.covariateBinary, data = data, family = binomial()) 
  summary(logisticFit.CoP.multi)
  p.logisticFit.CoP.multi.titer <- coef(summary(logisticFit.CoP.multi))[[2,4]]
  p.logisticFit.CoP.multi.covariateBinaryvalueB <- coef(summary(logisticFit.CoP.multi))[[3,4]]
  p.logisticFit.CoP.multi.interaction <- coef(summary(logisticFit.CoP.multi))[[4,4]]
  AIC.logistic.CoP.multi <- logisticFit.CoP.multi$aic
  # Overall VE
  VE.logisticFit.CoP.multi <- ve(logisticFit.CoP.multi, data, nboot = 2000)
  # Group A
  VE.logisticFit.CoP.multi.groupA <- ve(logisticFit.CoP.multi, data.groupA, nboot = 2000)
  # Group B
  VE.logisticFit.CoP.multi.groupB <- ve(logisticFit.CoP.multi, data.groupB, nboot = 2000)
  
  # Quadratic term for titer
  logisticFit.quadCoP.add <- glm(y.diseaseStatus ~ var.titer + var.quadtiter + var.covariateBinary, data = data, family = binomial()) 
  summary(logisticFit.quadCoP.add)
  p.logisticFit.quadCoP.add.titer <- coef(summary(logisticFit.quadCoP.add))[[2,4]]
  p.logisticFit.quadCoP.add.quadtiter <- coef(summary(logisticFit.quadCoP.add))[[3,4]]
  p.logisticFit.quadCoP.add.covariateBinaryvalueB <- coef(summary(logisticFit.quadCoP.add))[[4,4]]
  AIC.logistic.quadCoP.add <- logisticFit.quadCoP.add$aic
  # Overall VE
  VE.logisticFit.quadCoP.add <- ve(logisticFit.quadCoP.add, data, nboot = 2000)
  # Group A
  VE.logisticFit.quadCoP.add.groupA <- ve(logisticFit.quadCoP.add, data.groupA, nboot = 2000)
  # Group B
  VE.logisticFit.quadCoP.add.groupB <- ve(logisticFit.quadCoP.add, data.groupB, nboot = 2000)
  
  logisticFit.quadCoP.multi <- glm(y.diseaseStatus ~ var.titer + var.quadtiter*var.covariateBinary, data = data, family = binomial()) 
  summary(logisticFit.quadCoP.multi)
  p.logisticFit.quadCoP.multi.titer <- coef(summary(logisticFit.quadCoP.multi))[[2,4]]
  p.logisticFit.quadCoP.multi.quadtiter <- coef(summary(logisticFit.quadCoP.multi))[[3,4]]
  p.logisticFit.quadCoP.multi.covariateBinaryvalueB <- coef(summary(logisticFit.quadCoP.multi))[[4,4]]
  p.logisticFit.quadCoP.multi.interaction <- coef(summary(logisticFit.quadCoP.multi))[[5,4]]
  AIC.logistic.quadCoP.multi <- logisticFit.quadCoP.multi$aic
  # Overall VE
  VE.logisticFit.quadCoP.multi <- ve(logisticFit.quadCoP.multi, data, nboot = 2000)
  # Group A
  VE.logisticFit.quadCoP.multi.groupA <- ve(logisticFit.quadCoP.multi, data.groupA, nboot = 2000)
  # Group B
  VE.logisticFit.quadCoP.multi.groupB <- ve(logisticFit.quadCoP.multi, data.groupB, nboot = 2000)
 
  ## Logistic regression results ----------------------------------------------
  
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

  ## Compile all results ------------------------------------------------------
  result <- list(
    resultTrial = resultTrial, 
    resultCaseCount = resultCaseCount, 
    resultLogistic = resultLogistic,
    resultLogisticCoP = resultLogisticCoP,
    resultLogisticQuadCoP = resultLogisticQuadCoP
    )
}