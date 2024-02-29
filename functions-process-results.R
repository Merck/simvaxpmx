# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

### Supplementary functions

### Summarizing simulation results ----

analyzeSimulation <- function(SimulationData){
  
  Results <- data.frame(seed = SimulationData$results$resultTrial$seed,
                        NDiseased.all = SimulationData$results$resultTrial$NDiseased.all,
                        NDiseased.groupA = SimulationData$results$resultTrial$NDiseased.groupA,
                        NDiseased.groupB = SimulationData$results$resultTrial$NDiseased.groupB,
                        NNondiseased.all = SimulationData$results$resultTrial$NNondiseased.all,
                        NNondiseased.groupA = SimulationData$results$resultTrial$NNondiseased.groupA,
                        NNondiseased.groupB = SimulationData$results$resultTrial$NNondiseased.groupB,
                        trueVE.all = SimulationData$results$resultTrial$trueEfficacy.all,
                        trueVE.groupA = SimulationData$results$resultTrial$trueEfficacy.groupA,
                        trueVE.groupB = SimulationData$results$resultTrial$trueEfficacy.groupB,
                        #
                        casecountVE.all = SimulationData$results$resultCaseCount$VE.all,
                        casecountCI.LB.all = SimulationData$results$resultCaseCount$VECI.all.lowerBound,
                        casecountCI.UB.all = SimulationData$results$resultCaseCount$VECI.all.upperBound,
                        casecountVE.groupA = SimulationData$results$resultCaseCount$VE.groupA,
                        casecountCI.LB.groupA = SimulationData$results$resultCaseCount$VECI.groupA.lowerBound,
                        casecountCI.UB.groupA = SimulationData$results$resultCaseCount$VECI.groupA.upperBound,
                        casecountVE.groupB = SimulationData$results$resultCaseCount$VE.groupB,
                        casecountCI.LB.groupB = SimulationData$results$resultCaseCount$VECI.groupB.lowerBound,
                        casecountCI.UB.groupB = SimulationData$results$resultCaseCount$VECI.groupB.upperBound,
                        #
                        logisticVE.all.add = SimulationData$results$resultLogistic$VE.all.add/100,
                        logisticCI.LB.all.add = SimulationData$results$resultLogistic$VECI.all.add.LB/100,
                        logisticCI.UB.all.add = SimulationData$results$resultLogistic$VECI.all.add.UB/100,
                        logisticVE.all.multi = SimulationData$results$resultLogistic$VE.all.multi/100,
                        logisticCI.LB.all.multi = SimulationData$results$resultLogistic$VECI.all.multi.LB/100,
                        logisticCI.UB.all.multi = SimulationData$results$resultLogistic$VECI.all.multi.UB/100,
                        logisticVE.groupA.add = SimulationData$results$resultLogistic$VE.groupA.add/100,
                        logisticCI.LB.groupA.add = SimulationData$results$resultLogistic$VECI.groupA.add.LB/100,
                        logisticCI.UB.groupA.add = SimulationData$results$resultLogistic$VECI.groupA.add.UB/100,
                        logisticVE.groupA.multi = SimulationData$results$resultLogistic$VE.groupA.multi/100,
                        logisticCI.LB.groupA.multi = SimulationData$results$resultLogistic$VECI.groupA.multi.LB/100,
                        logisticCI.UB.groupA.multi = SimulationData$results$resultLogistic$VECI.groupA.multi.UB/100,
                        logisticVE.groupB.add = SimulationData$results$resultLogistic$VE.groupB.add/100,
                        logisticCI.LB.groupB.add = SimulationData$results$resultLogistic$VECI.groupB.add.LB/100,
                        logisticCI.UB.groupB.add = SimulationData$results$resultLogistic$VECI.groupB.add.UB/100,
                        logisticVE.groupB.multi = SimulationData$results$resultLogistic$VE.groupB.multi/100,
                        logisticCI.LB.groupB.multi = SimulationData$results$resultLogistic$VECI.groupB.multi.LB/100,
                        logisticCI.UB.groupB.multi = SimulationData$results$resultLogistic$VECI.groupB.multi.UB/100,
                        logisticP.VS.add = SimulationData$results$resultLogistic$pval.vaccStatus.add,
                        logisticP.A.add = SimulationData$results$resultLogistic$pval.covariateBinary.add,
                        logisticP.VS.multi = SimulationData$results$resultLogistic$pval.vaccStatus.multi,
                        logisticP.A.multi = SimulationData$results$resultLogistic$pval.covariateBinary.multi,
                        logisticP.inter.multi = SimulationData$results$resultLogistic$pval.interaction.multi,
                        logisticAIC.add = SimulationData$results$resultLogistic$aic.aic.add,
                        logisticAIC.multi = SimulationData$results$resultLogistic$aic.aic.multi,
                        #
                        logisticCoPVE.all.add = SimulationData$results$resultLogisticCoP$VE.all.add/100,
                        logisticCoPCI.LB.all.add = SimulationData$results$resultLogisticCoP$VECI.all.add.LB/100,
                        logisticCoPCI.UB.all.add = SimulationData$results$resultLogisticCoP$VECI.all.add.UB/100,
                        logisticCoPVE.all.multi = SimulationData$results$resultLogisticCoP$VE.all.multi/100,
                        logisticCoPCI.LB.all.multi = SimulationData$results$resultLogisticCoP$VECI.all.multi.LB/100,
                        logisticCoPCI.UB.all.multi = SimulationData$results$resultLogisticCoP$VECI.all.multi.UB/100,
                        logisticCoPVE.groupA.add = SimulationData$results$resultLogisticCoP$VE.groupA.add/100,
                        logisticCoPCI.LB.groupA.add = SimulationData$results$resultLogisticCoP$VECI.groupA.add.LB/100,
                        logisticCoPCI.UB.groupA.add = SimulationData$results$resultLogisticCoP$VECI.groupA.add.UB/100,
                        logisticCoPVE.groupA.multi = SimulationData$results$resultLogisticCoP$VE.groupA.multi/100,
                        logisticCoPCI.LB.groupA.multi = SimulationData$results$resultLogisticCoP$VECI.groupA.multi.LB/100,
                        logisticCoPCI.UB.groupA.multi = SimulationData$results$resultLogisticCoP$VECI.groupA.multi.UB/100,
                        logisticCoPVE.groupB.add = SimulationData$results$resultLogisticCoP$VE.groupB.add/100,
                        logisticCoPCI.LB.groupB.add = SimulationData$results$resultLogisticCoP$VECI.groupB.add.LB/100,
                        logisticCoPCI.UB.groupB.add = SimulationData$results$resultLogisticCoP$VECI.groupB.add.UB/100,
                        logisticCoPVE.groupB.multi = SimulationData$results$resultLogisticCoP$VE.groupB.multi/100,
                        logisticCoPCI.LB.groupB.multi = SimulationData$results$resultLogisticCoP$VECI.groupB.multi.LB/100,
                        logisticCoPCI.UB.groupB.multi = SimulationData$results$resultLogisticCoP$VECI.groupB.multi.UB/100,
                        logisticCoPP.T.add = SimulationData$results$resultLogisticCoP$pval.titer.add,
                        logisticCoPP.A.add = SimulationData$results$resultLogisticCoP$pval.covariateBinary.add,
                        logisticCoPP.T.multi = SimulationData$results$resultLogisticCoP$pval.titer.multi,
                        logisticCoPP.A.multi = SimulationData$results$resultLogisticCoP$pval.covariateBinary.multi,
                        logisticCoPP.inter.multi = SimulationData$results$resultLogisticCoP$pval.interaction.multi,
                        logisticCoPAIC.add = SimulationData$results$resultLogisticCoP$aic.aic.add,
                        logisticCoPAIC.multi = SimulationData$results$resultLogisticCoP$aic.aic.multi,
                        #
                        logisticQuadCoPVE.all.add = SimulationData$results$resultLogisticQuadCoP$VE.all.add/100,
                        logisticQuadCoPCI.LB.all.add = SimulationData$results$resultLogisticQuadCoP$VECI.all.add.LB/100,
                        logisticQuadCoPCI.UB.all.add = SimulationData$results$resultLogisticQuadCoP$VECI.all.add.UB/100,
                        logisticQuadCoPVE.all.multi = SimulationData$results$resultLogisticQuadCoP$VE.all.multi/100,
                        logisticQuadCoPCI.LB.all.multi = SimulationData$results$resultLogisticQuadCoP$VECI.all.multi.LB/100,
                        logisticQuadCoPCI.UB.all.multi = SimulationData$results$resultLogisticQuadCoP$VECI.all.multi.UB/100,
                        logisticQuadCoPVE.groupA.add = SimulationData$results$resultLogisticQuadCoP$VE.groupA.add/100,
                        logisticQuadCoPCI.LB.groupA.add = SimulationData$results$resultLogisticQuadCoP$VECI.groupA.add.LB/100,
                        logisticQuadCoPCI.UB.groupA.add = SimulationData$results$resultLogisticQuadCoP$VECI.groupA.add.UB/100,
                        logisticQuadCoPVE.groupA.multi = SimulationData$results$resultLogisticQuadCoP$VE.groupA.multi/100,
                        logisticQuadCoPCI.LB.groupA.multi = SimulationData$results$resultLogisticQuadCoP$VECI.groupA.multi.LB/100,
                        logisticQuadCoPCI.UB.groupA.multi = SimulationData$results$resultLogisticQuadCoP$VECI.groupA.multi.UB/100,
                        logisticQuadCoPVE.groupB.add = SimulationData$results$resultLogisticQuadCoP$VE.groupB.add/100,
                        logisticQuadCoPCI.LB.groupB.add = SimulationData$results$resultLogisticQuadCoP$VECI.groupB.add.LB/100,
                        logisticQuadCoPCI.UB.groupB.add = SimulationData$results$resultLogisticQuadCoP$VECI.groupB.add.UB/100,
                        logisticQuadCoPVE.groupB.multi = SimulationData$results$resultLogisticQuadCoP$VE.groupB.multi/100,
                        logisticQuadCoPCI.LB.groupB.multi = SimulationData$results$resultLogisticQuadCoP$VECI.groupB.multi.LB/100,
                        logisticQuadCoPCI.UB.groupB.multi = SimulationData$results$resultLogisticQuadCoP$VECI.groupB.multi.UB/100,
                        logisticQuadCoPP.T.add = SimulationData$results$resultLogisticQuadCoP$pval.titer.add,
                        logisticQuadCoPP.T2.add = SimulationData$results$resultLogisticQuadCoP$pval.quadtiter.add,
                        logisticQuadCoPP.A.add = SimulationData$results$resultLogisticQuadCoP$pval.covariateBinary.add,
                        logisticQuadCoPP.T.multi = SimulationData$results$resultLogisticQuadCoP$pval.titer.multi,
                        logisticQuadCoPP.T2.multi = SimulationData$results$resultLogisticQuadCoP$pval.quadtiter.multi,
                        logisticQuadCoPP.A.multi = SimulationData$results$resultLogisticQuadCoP$pval.covariateBinary.multi,
                        logisticQuadCoPP.inter.multi = SimulationData$results$resultLogisticQuadCoP$pval.interaction.multi,
                        logisticQuadCoPAIC.add = SimulationData$results$resultLogisticQuadCoP$aic.aic.add,
                        logisticQuadCoPAIC.multi = SimulationData$results$resultLogisticQuadCoP$aic.aic.multi
                        )

  ## Positive test for a covariate effect, typical approach ----
  # test is positive either if the interaction term in (in the multiplicative model) is significant
  # or id the covariate term (in the additive model) is significant
  Results <- Results %>% mutate(positive.typical = (logisticP.A.add < 0.05 | logisticP.inter.multi < 0.05))
  
  ## Positive test for a covariate effect, CoP-based approach, linear term for titer
  Results <- Results %>% 
    mutate(positive.CoP = (logisticCoPP.A.add < 0.05 | logisticCoPP.inter.multi < 0.05)) %>%
    mutate(positive.quadCoP = (logisticQuadCoPP.A.add < 0.05 | logisticQuadCoPP.inter.multi < 0.05))
  
  positive <- list("typical" = sum(Results$positive.typical)/nrow(Results),
                   "CoPbased.linear" = sum(Results$positive.CoP)/nrow(Results),
                   "CoPbased.quadratic" = sum(Results$positive.quadCoP)/nrow(Results))
  
  positive2500 <- list("typical" = sum(Results$positive.typical[1:2500])/2500,
                       "CoPbased.linear" = sum(Results$positive.CoP[1:2500])/2500,
                       "CoPbased.quadratic" = sum(Results$positive.quadCoP[1:2500])/1:2500)
  
  #### 1a Results based on typical approach, model without an interaction term ----
  
  ## Width of CI plot
  # Compare logistic additive model and case-count
  WidthDataset.1a <- prepareDataCIWidth(CILB1.groupA = Results$logisticCI.LB.groupA.add*100,
                                     CILB1.groupB = Results$logisticCI.LB.groupB.add*100,
                                     CIUB1.groupA = Results$logisticCI.UB.groupA.add*100,
                                     CIUB1.groupB = Results$logisticCI.UB.groupB.add*100,
                                     CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                     CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.logistic.add <- ggplotCIBA(WidthDataset.1a$data)
  narrower.casecount.logistic.add <- WidthDataset.1a$narrower
  
  #### 1b Results based on typical approach, model with an interaction term ----
  
  ## Width of CI plot
  # Compare logistic multiplicative model and case-count
  WidthDataset.1b <- prepareDataCIWidth(CILB1.groupA = Results$logisticCI.LB.groupA.multi*100,
                                     CILB1.groupB = Results$logisticCI.LB.groupB.multi*100,
                                     CIUB1.groupA = Results$logisticCI.UB.groupA.multi*100,
                                     CIUB1.groupB = Results$logisticCI.UB.groupB.multi*100,
                                     CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                     CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.logistic.multi <- ggplotCIBA(WidthDataset.1b$data)
  narrower.casecount.logistic.multi <- WidthDataset.1b$narrower
  
  #### 1final Results based on typical approach, selected model ----
  ## Model selection, typical approach
  Results <- Results %>% mutate(selected.typical.model = "multi")
  Results$selected.typical.model[Results$logisticAIC.add < Results$logisticAIC.multi] <- "add"
  
  selected.typical.model <- table(Results$selected.typical.model)
  
  ## VE
  Results <- Results %>% 
    mutate(selected.typical.VE.all = ifelse(selected.typical.model == "add", logisticVE.all.add, logisticVE.all.multi)) %>%
    mutate(selected.typical.CI.LB.all = ifelse(selected.typical.model == "add", logisticCI.LB.all.add, logisticCI.LB.all.multi)) %>%
    mutate(selected.typical.CI.UB.all = ifelse(selected.typical.model == "add", logisticCI.UB.all.add, logisticCI.UB.all.multi)) %>%
    mutate(selected.typical.VE.groupA = ifelse(selected.typical.model == "add", logisticVE.groupA.add, logisticVE.groupA.multi)) %>%
    mutate(selected.typical.CI.LB.groupA = ifelse(selected.typical.model == "add", logisticCI.LB.groupA.add, logisticCI.LB.groupA.multi)) %>%
    mutate(selected.typical.CI.UB.groupA = ifelse(selected.typical.model == "add", logisticCI.UB.groupA.add, logisticCI.UB.groupA.multi)) %>%
    mutate(selected.typical.VE.groupB = ifelse(selected.typical.model == "add", logisticVE.groupB.add, logisticVE.groupB.multi)) %>%
    mutate(selected.typical.VE.groupB = ifelse(selected.typical.model == "add", logisticVE.groupB.add, logisticVE.groupB.multi)) %>%
    mutate(selected.typical.CI.LB.groupB = ifelse(selected.typical.model == "add", logisticCI.LB.groupB.add, logisticCI.LB.groupB.multi)) %>%
    mutate(selected.typical.CI.UB.groupB = ifelse(selected.typical.model == "add", logisticCI.UB.groupB.add, logisticCI.UB.groupB.multi))
  
  # Width of CI plot
  # Compare case-count and logistic selected model
  WidthDataset.1f1 <- prepareDataCIWidth(CILB1.groupA = Results$selected.typical.CI.LB.groupA*100,
                                         CILB1.groupB = Results$selected.typical.CI.LB.groupB*100,
                                         CIUB1.groupA = Results$selected.typical.CI.UB.groupA*100,
                                         CIUB1.groupB = Results$selected.typical.CI.UB.groupB*100,
                                         CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                         CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                         CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                         CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.logistic.selected <- ggplotCIBA(WidthDataset.1f1$data)
  narrower.casecount.logistic.selected <- WidthDataset.1f1$narrower
  
  #### 2a Results based on linear CoP-based approach, model without an interaction term ----
  ## VE boxplots
  efficacyDF<- prepareDataVE(CoP.all = Results$logisticCoPVE.all.add,
                             typical.all = Results$logisticVE.all.add,
                             CoP.groupA = Results$logisticCoPVE.groupA.add,
                             typical.groupA = Results$logisticVE.groupA.add,
                             CoP.groupB = Results$logisticCoPVE.groupB.add,
                             typical.groupB = Results$logisticVE.groupB.add,
                             casecount.all = Results$casecountVE.all,
                             casecount.groupA =  Results$casecountVE.groupA,
                             casecount.groupB =  Results$casecountVE.groupB)
  boxplots.linear.add <- ggplotEfficacy(efficacyDF,
                                      Results$trueVE.all[1],
                                      Results$trueVE.groupA[1],
                                      Results$trueVE.groupB[1],
                                      "",
                                      "VE")
  
  ## Width of CI plot
  # Compare logistic additive model and logistic with linear titer additive model
  WidthDataset.2a1 <- prepareDataCIWidth(CILB1.groupA = Results$logisticCoPCI.LB.groupA.add*100,
                                     CILB1.groupB = Results$logisticCoPCI.LB.groupB.add*100,
                                     CIUB1.groupA = Results$logisticCoPCI.UB.groupA.add*100,
                                     CIUB1.groupB = Results$logisticCoPCI.UB.groupB.add*100,
                                     CILB2.groupA = Results$logisticCI.LB.groupA.add*100,
                                     CILB2.groupB = Results$logisticCI.LB.groupB.add*100,
                                     CIUB2.groupA = Results$logisticCI.UB.groupA.add*100,
                                     CIUB2.groupB = Results$logisticCI.UB.groupB.add*100)
  widthplot.logistic.CoP.linear.add <- ggplotCIBA(WidthDataset.2a1$data)
  narrower.logistic.CoP.linear.add <- WidthDataset.2a1$narrower
  
  # Compare case-count and logistic with linear titer additive model
  WidthDataset.2a2 <- prepareDataCIWidth(CILB1.groupA = Results$logisticCoPCI.LB.groupA.add*100,
                                     CILB1.groupB = Results$logisticCoPCI.LB.groupB.add*100,
                                     CIUB1.groupA = Results$logisticCoPCI.UB.groupA.add*100,
                                     CIUB1.groupB = Results$logisticCoPCI.UB.groupB.add*100,
                                     CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                     CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.CoP.linear.add <- ggplotCIBA(WidthDataset.2a2$data)
  narrower.casecount.CoP.linear.add <- WidthDataset.2a2$narrower
  
  #### 2b Results based on linear CoP-based approach, model with an interaction term ----
  ## VE boxplots
  efficacyDF<- prepareDataVE(CoP.all = Results$logisticCoPVE.all.multi,
                             typical.all = Results$logisticVE.all.multi,
                             CoP.groupA = Results$logisticCoPVE.groupA.multi,
                             typical.groupA = Results$logisticVE.groupA.multi,
                             CoP.groupB = Results$logisticCoPVE.groupB.multi,
                             typical.groupB = Results$logisticVE.groupB.multi,
                             casecount.all = Results$casecountVE.all,
                             casecount.groupA =  Results$casecountVE.groupA,
                             casecount.groupB =  Results$casecountVE.groupB)
  boxplots.linear.multi <- ggplotEfficacy(efficacyDF,
                                        Results$trueVE.all[1],
                                        Results$trueVE.groupA[1],
                                        Results$trueVE.groupB[1],
                                        "",
                                        "VE")
  ## Width of CI plot
  # Compare logistic multiplicative model and logistic with linear titer multiplicative model
  WidthDataset.2b1 <- prepareDataCIWidth(CILB1.groupA = Results$logisticCoPCI.LB.groupA.multi*100,
                                     CILB1.groupB = Results$logisticCoPCI.LB.groupB.multi*100,
                                    CIUB1.groupA = Results$logisticCoPCI.UB.groupA.multi*100,
                                     CIUB1.groupB = Results$logisticCoPCI.UB.groupB.multi*100,
                                     CILB2.groupA = Results$logisticCI.LB.groupA.multi*100,
                                     CILB2.groupB = Results$logisticCI.LB.groupB.multi*100,
                                     CIUB2.groupA = Results$logisticCI.UB.groupA.multi*100,
                                     CIUB2.groupB = Results$logisticCI.UB.groupB.multi*100)
  widthplot.logistic.CoP.linear.multi <- ggplotCIBA(WidthDataset.2b1$data)
  narrower.logistic.CoP.linear.multi <- WidthDataset.2b1$narrower
  
  # Compare case-count and logistic with linear titer multiplicative model
  WidthDataset.2b2 <- prepareDataCIWidth(CILB1.groupA = Results$logisticCoPCI.LB.groupA.multi*100,
                                     CILB1.groupB = Results$logisticCoPCI.LB.groupB.multi*100,
                                    CIUB1.groupA = Results$logisticCoPCI.UB.groupA.multi*100,
                                     CIUB1.groupB = Results$logisticCoPCI.UB.groupB.multi*100,
                                     CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                      CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.CoP.linear.multi <- ggplotCIBA(WidthDataset.2b2$data)
  narrower.casecount.CoP.linear.multi <- WidthDataset.2b2$narrower
  
  #### 2c Results based on quadratic CoP-based approach, model without an interaction term ----
  ## VE boxplots
  efficacyDF<- prepareDataVE(CoP.all = Results$logisticQuadCoPVE.all.add,
                             typical.all = Results$logisticVE.all.add,
                             CoP.groupA = Results$logisticQuadCoPVE.groupA.add,
                             typical.groupA = Results$logisticVE.groupA.add,
                             CoP.groupB = Results$logisticQuadCoPVE.groupB.add,
                             typical.groupB = Results$logisticVE.groupB.add,
                             casecount.all = Results$casecountVE.all,
                             casecount.groupA =  Results$casecountVE.groupA,
                             casecount.groupB =  Results$casecountVE.groupB)
  boxplots.quad.add <- ggplotEfficacy(efficacyDF,
                                        Results$trueVE.all[1],
                                        Results$trueVE.groupA[1],
                                        Results$trueVE.groupB[1],
                                        "",
                                        "VE")
  
  ## Width of CI plot
  # Compare logistic additive model and logistic with quadratic titer additive model
  WidthDataset.2c1 <- prepareDataCIWidth(CILB1.groupA = Results$logisticQuadCoPCI.LB.groupA.add*100,
                                     CILB1.groupB = Results$logisticQuadCoPCI.LB.groupB.add*100,
                                      CIUB1.groupA = Results$logisticQuadCoPCI.UB.groupA.add*100,
                                     CIUB1.groupB = Results$logisticQuadCoPCI.UB.groupB.add*100,
                                     CILB2.groupA = Results$logisticCI.LB.groupA.add*100,
                                     CILB2.groupB = Results$logisticCI.LB.groupB.add*100,
                                     CIUB2.groupA = Results$logisticCI.UB.groupA.add*100,
                                     CIUB2.groupB = Results$logisticCI.UB.groupB.add*100)
  widthplot.logistic.CoP.quad.add <- ggplotCIBA(WidthDataset.2c1$data)
  narrower.logistic.CoP.quad.add <- WidthDataset.2c1$narrower
  
  # Compare case-count and logistic with quadratic titer additive model
  WidthDataset.2c2 <- prepareDataCIWidth(CILB1.groupA = Results$logisticQuadCoPCI.LB.groupA.add*100,
                                     CILB1.groupB = Results$logisticQuadCoPCI.LB.groupB.add*100,
                                     CIUB1.groupA = Results$logisticQuadCoPCI.UB.groupA.add*100,
                                     CIUB1.groupB = Results$logisticQuadCoPCI.UB.groupB.add*100,
                                     CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                     CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.CoP.quad.add <- ggplotCIBA(WidthDataset.2c2$data)
  narrower.casecount.CoP.quad.add <- WidthDataset.2c2$narrower

  #### 2d Results based on quadratic CoP-based approach, model with an interaction term ----
  ## VE boxplots
  efficacyDF<- prepareDataVE(CoP.all = Results$logisticQuadCoPVE.all.multi,
                             typical.all = Results$logisticVE.all.multi,
                             CoP.groupA = Results$logisticQuadCoPVE.groupA.multi,
                             typical.groupA = Results$logisticVE.groupA.multi,
                             CoP.groupB = Results$logisticQuadCoPVE.groupB.multi,
                             typical.groupB = Results$logisticVE.groupB.multi,
                             casecount.all = Results$casecountVE.all,
                             casecount.groupA =  Results$casecountVE.groupA,
                             casecount.groupB =  Results$casecountVE.groupB)
  boxplots.quad.multi <- ggplotEfficacy(efficacyDF,
                                      Results$trueVE.all[1],
                                      Results$trueVE.groupA[1],
                                      Results$trueVE.groupB[1],
                                      "",
                                      "VE")
  
  ## Width of CI plot
  # Compare logistic multiplicative model and logistic with quadratic titer multiplicative model
  WidthDataset.2d1 <- prepareDataCIWidth(CILB1.groupA = Results$logisticQuadCoPCI.LB.groupA.multi*100,
                                     CILB1.groupB = Results$logisticQuadCoPCI.LB.groupB.multi*100,
                                     CIUB1.groupA = Results$logisticQuadCoPCI.UB.groupA.multi*100,
                                     CIUB1.groupB = Results$logisticQuadCoPCI.UB.groupB.multi*100,
                                     CILB2.groupA = Results$logisticCI.LB.groupA.multi*100,
                                     CILB2.groupB = Results$logisticCI.LB.groupB.multi*100,
                                     CIUB2.groupA = Results$logisticCI.UB.groupA.multi*100,
                                     CIUB2.groupB = Results$logisticCI.UB.groupB.multi*100)
  widthplot.logistic.CoP.quad.multi <- ggplotCIBA(WidthDataset.2d1$data)
  narrower.logistic.CoP.quad.multi <- WidthDataset.2d1$narrower
  
  # Compare logistic additive model and logistic with quadratic titer additive model
  WidthDataset.2d2 <- prepareDataCIWidth(CILB1.groupA = Results$logisticQuadCoPCI.LB.groupA.multi*100,
                                     CILB1.groupB = Results$logisticQuadCoPCI.LB.groupB.multi*100,
                                    CIUB1.groupA = Results$logisticQuadCoPCI.UB.groupA.multi*100,
                                     CIUB1.groupB = Results$logisticQuadCoPCI.UB.groupB.multi*100,
                                     CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                     CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.CoP.quad.multi <- ggplotCIBA(WidthDataset.2d2$data)
  narrower.casecount.CoP.quad.multi <- WidthDataset.2d2$narrower
  
  #### 2final Results based on CoP-based approach, selected model ----
  ## Model selection, CoP-based approach
  Results <- Results %>% mutate(minAIC = pmin(logisticCoPAIC.add, 
                                              logisticCoPAIC.multi,
                                              logisticQuadCoPAIC.add,
                                              logisticQuadCoPAIC.multi))
  
  
  Results$selected.CoP.model[Results$logisticCoPAIC.add == Results$minAIC] <- "linear add"
  Results$selected.CoP.model[Results$logisticCoPAIC.multi == Results$minAIC] <- "linear multi"
  Results$selected.CoP.model[Results$logisticQuadCoPAIC.add == Results$minAIC] <- "quad add"
  Results$selected.CoP.model[Results$logisticQuadCoPAIC.multi == Results$minAIC] <- "quad multi"
  
  selected.CoP.model <- table(Results$selected.CoP.model)
  
  ## VE
  Results <- Results %>% 
    mutate(selected.CoP.VE.all = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.all.add, 
                                           selected.CoP.model == "linear multi" ~ logisticCoPVE.all.multi,
                                           selected.CoP.model == "quad add" ~ logisticQuadCoPVE.all.add,
                                           selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.all.multi,
                                           TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.CI.LB.all = case_when(selected.CoP.model == "linear add" ~ logisticCoPCI.LB.all.add, 
                                              selected.CoP.model == "linear multi" ~ logisticCoPCI.LB.all.multi,
                                              selected.CoP.model == "quad add" ~ logisticQuadCoPCI.LB.all.add,
                                              selected.CoP.model == "quad multi" ~ logisticQuadCoPCI.LB.all.multi,
                                              TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.CI.UB.all = case_when(selected.CoP.model == "linear add" ~ logisticCoPCI.UB.all.add, 
                                              selected.CoP.model == "linear multi" ~ logisticCoPCI.UB.all.multi,
                                              selected.CoP.model == "quad add" ~ logisticQuadCoPCI.UB.all.add,
                                              selected.CoP.model == "quad multi" ~ logisticQuadCoPCI.UB.all.multi,
                                              TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.VE.groupA = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.groupA.add, 
                                              selected.CoP.model == "linear multi" ~ logisticCoPVE.groupA.multi,
                                              selected.CoP.model == "quad add" ~ logisticQuadCoPVE.groupA.add,
                                              selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.groupA.multi,
                                              TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.CI.LB.groupA = case_when(selected.CoP.model == "linear add" ~ logisticCoPCI.LB.groupA.add, 
                                                 selected.CoP.model == "linear multi" ~ logisticCoPCI.LB.groupA.multi,
                                                 selected.CoP.model == "quad add" ~ logisticQuadCoPCI.LB.groupA.add,
                                                 selected.CoP.model == "quad multi" ~ logisticQuadCoPCI.LB.groupA.multi,
                                                 TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.CI.UB.groupA = case_when(selected.CoP.model == "linear add" ~ logisticCoPCI.UB.groupA.add, 
                                                 selected.CoP.model == "linear multi" ~ logisticCoPCI.UB.groupA.multi,
                                                 selected.CoP.model == "quad add" ~ logisticQuadCoPCI.UB.groupA.add,
                                                 selected.CoP.model == "quad multi" ~ logisticQuadCoPCI.UB.groupA.multi,
                                                 TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.VE.groupB = case_when(selected.CoP.model == "linear add" ~ logisticCoPVE.groupB.add,
                                              selected.CoP.model == "linear multi" ~ logisticCoPVE.groupB.multi,
                                              selected.CoP.model == "quad add" ~ logisticQuadCoPVE.groupB.add,
                                              selected.CoP.model == "quad multi" ~ logisticQuadCoPVE.groupB.multi,
                                              TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.CI.LB.groupB = case_when(selected.CoP.model == "linear add" ~ logisticCoPCI.LB.groupB.add,
                                                 selected.CoP.model == "linear multi" ~ logisticCoPCI.LB.groupB.multi,
                                                 selected.CoP.model == "quad add" ~ logisticQuadCoPCI.LB.groupB.add,
                                                 selected.CoP.model == "quad multi" ~ logisticQuadCoPCI.LB.groupB.multi,
                                                 TRUE ~ NA_real_)) %>%
    mutate(selected.CoP.CI.UB.groupB = case_when(selected.CoP.model == "linear add" ~ logisticCoPCI.UB.groupB.add,
                                                 selected.CoP.model == "linear multi" ~ logisticCoPCI.UB.groupB.multi,
                                                 selected.CoP.model == "quad add" ~ logisticQuadCoPCI.UB.groupB.add,
                                                 selected.CoP.model == "quad multi" ~ logisticQuadCoPCI.UB.groupB.multi,
                                                 TRUE ~ NA_real_))
  
  ## VE boxplots
  efficacyDF<- prepareDataVE(CoP.all = Results$selected.CoP.VE.all,
                             typical.all = Results$selected.typical.VE.all,
                             CoP.groupA = Results$selected.CoP.VE.groupA,
                             typical.groupA = Results$selected.typical.VE.groupA,
                             CoP.groupB = Results$selected.CoP.VE.groupB,
                             typical.groupB = Results$selected.typical.VE.groupB,
                             casecount.all = Results$casecountVE.all,
                             casecount.groupA =  Results$casecountVE.groupA,
                             casecount.groupB =  Results$casecountVE.groupB)
  boxplots.selected <- ggplotEfficacy(efficacyDF,
                                      Results$trueVE.all[1],
                                      Results$trueVE.groupA[1],
                                      Results$trueVE.groupB[1],
                                      "",
                                      "VE")

  # Width of CI plot
  # Compare logistic selected model and logistic with titer selected model
  WidthDataset.2f1 <- prepareDataCIWidth(CILB1.groupA = Results$selected.CoP.CI.LB.groupA*100,
                                     CILB1.groupB = Results$selected.CoP.CI.LB.groupB*100,
                                     CIUB1.groupA = Results$selected.CoP.CI.UB.groupA*100,
                                     CIUB1.groupB = Results$selected.CoP.CI.UB.groupB*100,
                                     CILB2.groupA = Results$selected.typical.CI.LB.groupA*100,
                                     CILB2.groupB = Results$selected.typical.CI.LB.groupB*100,
                                    CIUB2.groupA = Results$selected.typical.CI.UB.groupA*100,
                                     CIUB2.groupB = Results$selected.typical.CI.UB.groupB*100)
  widthplot.logistic.CoP.selected <- ggplotCIBA(WidthDataset.2f1$data)
  narrower.logistic.CoP.selected <- WidthDataset.2f1$narrower
  
  # Compare case-count and logistic with titer selected model
  WidthDataset.2f2 <- prepareDataCIWidth(CILB1.groupA = Results$selected.CoP.CI.LB.groupA*100,
                                     CILB1.groupB = Results$selected.CoP.CI.LB.groupB*100,
                                    CIUB1.groupA = Results$selected.CoP.CI.UB.groupA*100,
                                     CIUB1.groupB = Results$selected.CoP.CI.UB.groupB*100,
                                    CILB2.groupA = Results$casecountCI.LB.groupA*100,
                                     CILB2.groupB = Results$casecountCI.LB.groupB*100,
                                      CIUB2.groupA = Results$casecountCI.UB.groupA*100,
                                     CIUB2.groupB = Results$casecountCI.UB.groupB*100)
  widthplot.casecount.CoP.selected <- ggplotCIBA(WidthDataset.2f2$data)
  narrower.casecount.CoP.selected <- WidthDataset.2f2$narrower

  #### Summary of results ----
  ## Boxplots of VE
  boxplots <- list(linear.add = boxplots.linear.add,
                   linear.multi = boxplots.linear.multi,
                   quad.add = boxplots.quad.add,
                   quad.multi = boxplots.quad.multi,
                   selected = boxplots.selected)
  
  ## Widths of VE CI
  widthdatasets <- list(casecount.logistic.selected = WidthDataset.1f1,
                     casecount.CoP.selected = WidthDataset.2f2,
                     logistic.CoP.selected = WidthDataset.2f1,
                     casecount.logistic.add = WidthDataset.1a,
                     casecount.logistic.multi = WidthDataset.1b,
                     logistic.CoP.linear.add = WidthDataset.2a1,
                     casecount.CoP.linear.add = WidthDataset.2a2,
                     logistic.CoP.linear.multi = WidthDataset.2b1,
                     casecount.CoP.linear.multi = WidthDataset.2b2,
                     logistic.CoP.quad.add = WidthDataset.2c1,
                     casecount.CoP.quad.add = WidthDataset.2c2,
                     logistic.CoP.quad.multi = WidthDataset.2d1,
                     casecount.CoP.quad.multi = WidthDataset.2d2
  )
  
  
  widthplots <- list(casecount.logistic.selected = widthplot.casecount.logistic.selected,
                     casecount.CoP.selected = widthplot.casecount.CoP.selected,
                     logistic.CoP.selected = widthplot.logistic.CoP.selected,
                     casecount.logistic.add = widthplot.casecount.logistic.add,
                     casecount.logistic.multi = widthplot.casecount.logistic.multi,
                     logistic.CoP.linear.add = widthplot.logistic.CoP.linear.add,
                     casecount.CoP.linear.add = widthplot.casecount.CoP.linear.add,
                     logistic.CoP.linear.multi = widthplot.logistic.CoP.linear.multi,
                     casecount.CoP.linear.multi = widthplot.casecount.CoP.linear.multi,
                     logistic.CoP.quad.add = widthplot.logistic.CoP.quad.add,
                     casecount.CoP.quad.add = widthplot.casecount.CoP.quad.add,
                     logistic.CoP.quad.multi = widthplot.logistic.CoP.quad.multi,
                     casecount.CoP.quad.multi = widthplot.casecount.CoP.quad.multi
                     )
  
  narrower <- list(casecount.logistic.selected = narrower.casecount.logistic.selected,
                   casecount.CoP.selected = narrower.casecount.CoP.selected,
                   logistic.CoP.selected = narrower.logistic.CoP.selected,
                   casecount.logistic.add = narrower.casecount.logistic.add,
                    casecount.logistic.multi = narrower.casecount.logistic.multi,
                     logistic.CoP.linear.add = narrower.logistic.CoP.linear.add,
                    casecount.CoP.linear.add = narrower.casecount.CoP.linear.add,
                     logistic.CoP.linear.multi = narrower.logistic.CoP.linear.multi,
                   casecount.CoP.linear.multi = narrower.casecount.CoP.linear.multi,
                    logistic.CoP.quad.add = narrower.logistic.CoP.quad.add,
                     casecount.CoP.quad.add = narrower.casecount.CoP.quad.add,
                    logistic.CoP.quad.multi = narrower.logistic.CoP.quad.multi,
                     casecount.CoP.quad.multi = narrower.casecount.CoP.quad.multi
  )
  
  ## Median VE
  medianVE <- list("casecount" = list(all = median(Results$casecountVE.all),
                                      groupA =  median(Results$casecountVE.groupA),
                                      groupB =  median(Results$casecountVE.groupB)),
                   "typical.add" = list(all = median(Results$logisticVE.all.add),
                                         groupA = median(Results$logisticVE.groupA.add),
                                         groupB = median(Results$logisticVE.groupB.add)),
                   "typical.multi" = list(all = median(Results$logisticVE.all.multi),
                                             groupA = median(Results$logisticVE.groupA.multi),
                                             groupB = median(Results$logisticVE.groupB.multi)),
                   "typical.selected" = list(all = median(Results$selected.typical.VE.all),
                                             groupA = median(Results$selected.typical.VE.groupA),
                                             groupB = median(Results$selected.typical.VE.groupB)),
                   "CoP.linear.add" = list(all = median(Results$logisticCoPVE.all.add),
                                           groupA = median(Results$logisticCoPVE.groupA.add),
                                           groupB = median(Results$logisticCoPVE.groupB.add)),
                   "CoP.linear.multi" = list(all = median(Results$logisticCoPVE.all.multi),
                                         groupA = median(Results$logisticCoPVE.groupA.multi),
                                         groupB = median(Results$logisticCoPVE.groupB.multi)),
                   "CoP.quad.add" = list(all = median(Results$logisticQuadCoPVE.all.add),
                                         groupA = median(Results$logisticQuadCoPVE.groupA.add),
                                         groupB = median(Results$logisticQuadCoPVE.groupB.add)),
                   "CoP.quad.multi" = list(all = median(Results$logisticQuadCoPVE.all.multi),
                                         groupA = median(Results$logisticQuadCoPVE.groupA.multi),
                                         groupB = median(Results$logisticQuadCoPVE.groupB.multi)),
                   "CoP.selected" = list(all = median(Results$selected.CoP.VE.all),
                                         groupA = median(Results$selected.CoP.VE.groupA),
                                         groupB = median(Results$selected.CoP.VE.groupB)))

  ## Coverage of VE CI
  coverage <- list("casecount" = list(all = sum((Results$casecountCI.LB.all < Results$trueVE.all[1])&(Results$casecountCI.UB.all > Results$trueVE.all[1]))/nrow(Results),
                                      groupA = sum((Results$casecountCI.LB.groupA < Results$trueVE.groupA[1])&(Results$casecountCI.UB.groupA > Results$trueVE.groupA[1]))/nrow(Results),
                                      groupB = sum((Results$casecountCI.LB.groupB < Results$trueVE.groupB[1])&(Results$casecountCI.UB.groupB > Results$trueVE.groupB[1]))/nrow(Results)),
                   "typical.add" = list(all = sum((Results$logisticCI.LB.all.add < Results$trueVE.all[1])&(Results$logisticCI.UB.all.add > Results$trueVE.all[1]))/nrow(Results),
                                             groupA = sum((Results$logisticCI.LB.groupA.add < Results$trueVE.groupA[1])&(Results$logisticCI.UB.groupA.add > Results$trueVE.groupA[1]))/nrow(Results),
                                             groupB = sum((Results$logisticCI.LB.groupB.add < Results$trueVE.groupB[1])&(Results$logisticCI.UB.groupB.add > Results$trueVE.groupB[1]))/nrow(Results)),
                   "typical.multi" = list(all = sum((Results$logisticCI.LB.all.multi < Results$trueVE.all[1])&(Results$logisticCI.UB.all.multi > Results$trueVE.all[1]))/nrow(Results),
                                             groupA = sum((Results$logisticCI.LB.groupA.multi < Results$trueVE.groupA[1])&(Results$logisticCI.UB.groupA.multi > Results$trueVE.groupA[1]))/nrow(Results),
                                             groupB = sum((Results$logisticCI.LB.groupB.multi < Results$trueVE.groupB[1])&(Results$logisticCI.UB.groupB.multi > Results$trueVE.groupB[1]))/nrow(Results)),
                   "typical.selected" = list(all = sum((Results$selected.typical.CI.LB.all < Results$trueVE.all[1])&(Results$selected.typical.CI.UB.all > Results$trueVE.all[1]))/nrow(Results),
                                             groupA = sum((Results$selected.typical.CI.LB.groupA < Results$trueVE.groupA[1])&(Results$selected.typical.CI.UB.groupA > Results$trueVE.groupA[1]))/nrow(Results),
                                             groupB = sum((Results$selected.typical.CI.LB.groupB < Results$trueVE.groupB[1])&(Results$selected.typical.CI.UB.groupB > Results$trueVE.groupB[1]))/nrow(Results)),
                   "CoP.linear.add" = list(all = sum((Results$logisticCoPCI.LB.all.add < Results$trueVE.all[1])&(Results$logisticCoPCI.UB.all.add > Results$trueVE.all[1]))/nrow(Results),
                                         groupA = sum((Results$logisticCoPCI.LB.groupA.add < Results$trueVE.groupA[1])&(Results$logisticCoPCI.UB.groupA.add > Results$trueVE.groupA[1]))/nrow(Results),
                                         groupB = sum((Results$logisticCoPCI.LB.groupB.add < Results$trueVE.groupB[1])&(Results$logisticCoPCI.UB.groupB.add > Results$trueVE.groupB[1]))/nrow(Results)),
                   "CoP.linear.multi" = list(all = sum((Results$logisticCoPCI.LB.all.multi < Results$trueVE.all[1])&(Results$logisticCoPCI.UB.all.multi > Results$trueVE.all[1]))/nrow(Results),
                                         groupA = sum((Results$logisticCoPCI.LB.groupA.multi < Results$trueVE.groupA[1])&(Results$logisticCoPCI.UB.groupA.multi > Results$trueVE.groupA[1]))/nrow(Results),
                                         groupB = sum((Results$logisticCoPCI.LB.groupB.multi < Results$trueVE.groupB[1])&(Results$logisticCoPCI.UB.groupB.multi > Results$trueVE.groupB[1]))/nrow(Results)),
                   "CoP.quad.add" = list(all = sum((Results$logisticQuadCoPCI.LB.all.add < Results$trueVE.all[1])&(Results$logisticQuadCoPCI.UB.all.add > Results$trueVE.all[1]))/nrow(Results),
                                           groupA = sum((Results$logisticCoPCI.LB.groupA.add < Results$trueVE.groupA[1])&(Results$logisticCoPCI.UB.groupA.add > Results$trueVE.groupA[1]))/nrow(Results),
                                           groupB = sum((Results$logisticQuadCoPCI.LB.groupB.add < Results$trueVE.groupB[1])&(Results$logisticQuadCoPCI.UB.groupB.add > Results$trueVE.groupB[1]))/nrow(Results)),
                   "CoP.quad.multi" = list(all = sum((Results$logisticQuadCoPCI.LB.all.multi < Results$trueVE.all[1])&(Results$logisticQuadCoPCI.UB.all.multi > Results$trueVE.all[1]))/nrow(Results),
                                             groupA = sum((Results$logisticQuadCoPCI.LB.groupA.multi < Results$trueVE.groupA[1])&(Results$logisticQuadCoPCI.UB.groupA.multi > Results$trueVE.groupA[1]))/nrow(Results),
                                             groupB = sum((Results$logisticQuadCoPCI.LB.groupB.multi < Results$trueVE.groupB[1])&(Results$logisticQuadCoPCI.UB.groupB.multi > Results$trueVE.groupB[1]))/nrow(Results)),
                   "CoP.selected" = list(all = sum((Results$selected.CoP.CI.LB.all < Results$trueVE.all[1])&(Results$selected.CoP.CI.UB.all > Results$trueVE.all[1]))/nrow(Results),
                                         groupA = sum((Results$selected.CoP.CI.LB.groupA < Results$trueVE.groupA[1])&(Results$selected.CoP.CI.UB.groupA > Results$trueVE.groupA[1]))/nrow(Results),
                                         groupB = sum((Results$selected.CoP.CI.LB.groupB < Results$trueVE.groupB[1])&(Results$selected.CoP.CI.UB.groupB > Results$trueVE.groupB[1]))/nrow(Results)))
                   
  
  ## Truth used for simulating trial data
  truth <- SimulationData$truth
  
  ## True VE
  trueVE <- list("all" = Results$trueVE.all[1],
                 "groupA" = Results$trueVE.groupA[1],
                 "groupB" = Results$trueVE.groupB[1])

  ## Summarizing data to be used by the roc function (from pROC package)
  df <- data.frame(Results$logisticCoPP.inter.multi,Results$logisticCoPP.A.add)
  predictor.roc.CoP <- apply(df, 1, FUN = min)
  dff <- data.frame(Results$logisticP.inter.multi,Results$logisticP.A.add)
  predictor.roc.typical <- apply(dff, 1, FUN = min)

  predictors.roc <- list("CoP" = predictor.roc.CoP,
                         "typical" = predictor.roc.typical)
  
  return(list("trueVE" = trueVE,
              "medianVE" = medianVE,
              "coverage" = coverage,
              "positive" = positive,
              "positive2500" = positive2500,
              "boxplots" = boxplots,
              "truth" = truth,
              "predictors.roc" = predictors.roc,
              "widthplots" = widthplots,
              "widthdatasets" = widthdatasets,
              "narrowerCI" = narrower,
              "selected.typical.model" = selected.typical.model,
              "selected.CoP.model" = selected.CoP.model,
              "Results" = Results))
}

### Numerical results ----

## Coverage
# Calculates coverage from lower & upper CI, in comparison with true efficacy - based on PoDBA package data structure
CoverageFun <- function(LowerCI,
                        UpperCI,
                        TrueEfficacy,
                        NumberOfSimulations){

  coverage <- sum(LowerCI < TrueEfficacy & UpperCI > TrueEfficacy)/NumberOfSimulations
  return(coverage)
}

# Coverage of both case-count Efficacay and PoDBA efficacy - based on PoDBA package data structure
CoverageDataset <- function(results1, results2, ClinicalTrial, includeNAs = T, DecimalRound = 2, group = "all", mod = "add", method1, method2){

  CILB1 <- as.numeric(unlist(results1[paste0("VECI.",group,".",mod,".LB")]))/100
  CILB2 <- unlist(results2[paste0("VECI.",group,".lowerBound")])
  CIUB1 <- as.numeric(unlist(results1[paste0("VECI.",group,".",mod,".UB")]))/100
  CIUB2 <- unlist(results2[paste0("VECI.",group,".upperBound")])
  TrueEfficacy <- unlist(ClinicalTrial[paste0("trueEfficacy.",group)])
    
  if(includeNAs){
  noSim <- length(TrueEfficacy)
  }else{
    noSim <- length(na.omit(TrueEfficacy))
  }

  results1$VECI.CILow <- CILB1[!is.na(TrueEfficacy)]
  results1$VECI.CIHigh <- CIUB1[!is.na(TrueEfficacy)]
  results2$VECI.lowerBound <- CILB2[!is.na(TrueEfficacy)]
  results2$VECI.upperBound <- CIUB2[!is.na(TrueEfficacy)]
  ClinicalTrial$trueEfficacy <- ClinicalTrial$trueEfficacy.all[!is.na(TrueEfficacy)]

  data <- c(results1[c("VECI.CILow",
                         "VECI.CIHigh")],
            results2[c("VECI.lowerBound",
                            "VECI.upperBound")],
            ClinicalTrial[paste0("trueEfficacy.",group)])
  
  data <- lapply(data, function(x) ifelse(is.na(x), 1, x) )

  Coverage1Eff95 <- with(data,
                             CoverageFun(VECI.CILow,
                                         VECI.CIHigh,
                                         unique(na.omit(TrueEfficacy)),
                                         noSim
                             ))

  Coverage2Eff95   <- with(data,
                             CoverageFun(VECI.lowerBound,
                                         VECI.upperBound,
                                         unique(na.omit(TrueEfficacy)),
                                         noSim
                             ))

  out <- list()
  name1 <- paste0("coverage.", method1,".",group)
  name2 <- paste0("coverage.", method2,".",group)
  out[[name1]] <- round(Coverage1Eff95*100,DecimalRound)
  out[[name2]] <- round(Coverage2Eff95*100,DecimalRound)
  return(out)
}

### Visualization ----
size <- 12
themePlot <-
  theme(
    plot.title = element_text(size=size,hjust = 0.5),
    axis.text.y   = element_text(size=size),
    axis.text.x   = element_text(size=size),
    axis.title.y  = element_text(size=size),
    axis.title.x  = element_text(size=size),
    legend.text = element_text(size = size),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5)
  )

## Efficacy boxplot ----

# function for data preparation
prepareDataVE <- function(CoP.all,
                          typical.all,
                          CoP.groupA,
                          typical.groupA,
                          CoP.groupB,
                          typical.groupB,
                          casecount.all = Results$casecountVE.all,
                          casecount.groupA =  Results$casecountVE.groupA,
                          casecount.groupB =  Results$casecountVE.groupB){

  df <- data.frame(value = c(CoP.all,
                             typical.all,
                             CoP.groupA,
                             typical.groupA,
                             CoP.groupB,
                             typical.groupB,
                             casecount.all,
                             casecount.groupA,
                             casecount.groupB),
                   type = c(rep("CoP-based", length(casecount.all)),
                            rep("typical",  length(casecount.all)),
                            rep("CoP-based", length(casecount.all)),
                            rep("typical",  length(casecount.all)),
                            rep("CoP-based", length(casecount.all)),
                            rep("typical",  length(casecount.all)),
                            rep("case-count",  3*length(casecount.all))),
                   group = c(rep("overall", 2*length(casecount.all)),
                             rep("groupA", length(casecount.all)),
                             rep("groupA",length(casecount.all)),
                             rep("groupB", length(casecount.all)),
                             rep("groupB",length(casecount.all)),
                             rep("overall", length(casecount.all)),
                             rep("groupA",length(casecount.all)),
                             rep("groupB", length(casecount.all)))
  )
 df$type <- factor(df$type,levels = c("CoP-based","typical","case-count"),ordered = TRUE)
 df$group <- factor(df$group,levels = c("groupA","groupB","overall"),ordered = TRUE)
 return(df)
}

wprepareDataVE <- function(CoP.all,
                          typical.all,
                          CoP.groupA,
                          typical.groupA,
                          CoP.groupB,
                          typical.groupB){
  
  df <- data.frame(value = c(CoP.all,
                             typical.all,
                             CoP.groupA,
                             typical.groupA,
                             CoP.groupB,
                             typical.groupB),
                   type = c(rep("CoP-based", length(CoP.all)),
                            rep("typical",  length(typical.all)),
                            rep("CoP-based", length(CoP.groupA)),
                            rep("typical",  length(typical.groupA)),
                            rep("CoP-based", length(CoP.groupB)),
                            rep("typical",  length(typical.groupB))),
                   group = c(rep("overall", length(CoP.all)),
                             rep("overall", length(typical.all)),
                             rep("groupA", length(CoP.groupA)),
                             rep("groupA",length(typical.groupA)),
                             rep("groupB", length(CoP.groupB)),
                             rep("groupB",length(typical.groupB)))
  )
  df$type <- factor(df$type,levels = c("CoP-based","typical"),ordered = TRUE)
  df$group <- factor(df$group,levels = c("groupA","groupB","overall"),ordered = TRUE)
  return(df)
}

# ggplot function for efficacy point estimates
ggplotEfficacy <- function(data, trueVE, trueVE.groupA, trueVE.groupB, title = " ", ylabel = "Efficacy") {
  ggplot(data, aes(x = type, y = value, col=group)) +

    geom_boxplot(outlier.shape= 1, outlier.size=1, outlier.fill = NULL) +
    # geom_jitter(aes(x = type, y = value, fill = NULL), width = 0.15, shape = 1, size = 2) +
    geom_hline(aes(yintercept = trueVE.groupA),
               color = "darkslategray3",
               size = 0.5,
               linetype = "dashed") +
    geom_hline(aes(yintercept = trueVE.groupB),
               color = "gray20",
               size = 0.5,
               linetype = "dashed") +
    geom_hline(aes(yintercept = trueVE),
               color = "gray",
               size = 0.5,
               linetype = "dashed") +
    labs(x = NULL, y = ylabel) +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                       limits = c(-5,1)) +
    ggtitle(title) +
    themePlot +
    scale_color_manual(values=c("darkslategray3", "gray20", "gray"))
}

# ggplot function for efficacy point estimates using "wrongly" selected model
wggplotEfficacy <- function(data, trueVE, trueVE.groupA, trueVE.groupB, title = " ", ylabel = "Efficacy") {
  ggplot(data, aes(x = type, y = value, col=group)) +
    
    geom_violin(outlier.shape= 1, outlier.size=1, outlier.fill = NULL) +
    # geom_jitter(aes(x = type, y = value, fill = NULL), width = 0.15, shape = 1, size = 2) +
    geom_hline(aes(yintercept = trueVE.groupA),
               color = "darkslategray3",
               size = 0.5,
               linetype = "dashed") +
    geom_hline(aes(yintercept = trueVE.groupB),
               color = "gray20",
               size = 0.5,
               linetype = "dashed") +
    geom_hline(aes(yintercept = trueVE),
               color = "gray",
               size = 0.5,
               linetype = "dashed") +
    labs(x = NULL, y = ylabel) +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                       limits = c(-5,1)) +
    ggtitle(title) +
    themePlot +
    scale_color_manual(values=c("darkslategray3", "gray20", "gray")) 
}

## Widths of VE CI ----
prepareDataCIWidth <- function(CILB1.groupA,
                               CILB1.groupB,
                               CIUB1.groupA,
                               CIUB1.groupB,
                               CILB2.groupA,
                               CILB2.groupB,
                               CIUB2.groupA,
                               CIUB2.groupB){
  

  width1.groupA <- CIUB1.groupA-CILB1.groupA
  width1.groupB <- CIUB1.groupB-CILB1.groupB
  width2.groupA <- CIUB2.groupA-CILB2.groupA
  width2.groupB <- CIUB2.groupB-CILB2.groupB
  
  
  df <- data.frame(width1 = unlist(c(width1.groupA,
                                     width1.groupB)),
                   width2 = unlist(c(width2.groupA,
                                     width2.groupB)),
                   group = c(rep("groupA", length(width1.groupA)),
                             rep("groupB", length(width1.groupB)))
  )
  df$group <- factor(df$group,levels = c("groupA","groupB"))
  
  narrower.groupA <- sum(width1.groupA<width2.groupA, na.rm = TRUE)/length(width1.groupA)
  narrower.groupB <- sum(width1.groupB<width2.groupB, na.rm = TRUE)/length(width1.groupB)
  narrower <- list("groupA" = narrower.groupA,
                   "groupB" = narrower.groupB)
  
  return(list("data" = df,
              "narrower" = narrower,
              "median.width1.groupA" = median(width1.groupA),
              "median.width1.groupB" = median(width1.groupB),
              "median.width2.groupA" = median(width2.groupA),
              "median.width2.groupB" = median(width2.groupB),
              "median.diff.groupA" = median(width2.groupA-width1.groupA),
              "median.diff.groupB" = median(width2.groupB-width1.groupB))
  )
}

ggplotCIWidth <- function(df){
  ggplot(df, aes(x = width1, y = width2, color=group)) +
    geom_point(size = 0.5) +
    geom_abline(slope=1, intercept=0)+
    coord_cartesian(ylim = c(0,100), xlim=c(0,100)) +
    # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
    #                    limits = c(-5,1)) +
    themePlot +
    scale_color_manual(values=c("darkslategray3", "gray20", "gray"))
}

ggplotCIBA <- function(df){
  ggplot(df, aes(x = width1, y = width2-width1, color=group)) +
    geom_point(size = 0.1) +
    geom_hline(aes(yintercept=0)) +
    #geom_abline(slope=1, intercept=0)+
    coord_cartesian(ylim = c(-20,80), xlim=c(0,80)) +
    # scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
    #                    limits = c(-5,1)) +
    themePlot +
    scale_color_manual(values=c("darkslategray3", "gray20", "gray"))
}

## Ridgeline plots ----
ridgelineSummary <- function(results, truth, trueVE, title){
  # Point estimates of VE, pmax, gamma, et50 normalized using the true values
  # Point estimates of pmax.covariate, gamma.covariate, et50.covariate (without normalization)
  val <- c(results$VE.all/trueVE,
           results$PoD.param.pmax.covariateBinary,
           results$PoD.param.pmax/truth$PoDModel$valueParams$param.pmax,
           results$PoD.param.gamma.covariateBinary,
           results$PoD.param.gamma/truth$PoDModel$valueParams$param.gamma,
           results$PoD.param.et50.covariateBinary,
           results$PoD.param.et50/truth$PoDModel$valueParams$param.et50)
  name <- c(rep("VE", length(results$VE.all)),
            rep("pmax_ageGroup", length(results$PoD.param.pmax.covariateBinary)),
            rep("pmax", length(results$PoD.param.pmax)),
            rep("gamma_ageGroup", length(results$PoD.param.gamma.covariateBinary)),
            rep("gamma", length(results$PoD.param.gamma)),
            rep("et50_ageGroup", length(results$PoD.param.et50.covariateBinary)),
            rep("et50", length(results$PoD.param.et50)))
  
  df <- data.frame(value = val,
                   variable = name)
  
  truth.lines <- data.frame(variable = name,
                            x0 = c(1,
                                   truth$PoDModel$valueParams$param.pmax.covariateBinary,
                                   1,
                                   truth$PoDModel$valueParams$param.gamma.covariateBinary,
                                   1,
                                   truth$PoDModel$valueParams$param.et50.covariateBinary,
                                   1))
    
  ggplot(df,
         aes(x = value, y = variable, fill = stat(abs(1-x)))) +
    geom_density_ridges_gradient(scale = 1.5, alpha=0.5, stat="binline", bins=200) +
    scale_fill_viridis(name = "value") +
    theme_ridges(grid = FALSE) + 
    theme(legend.position = "none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 8)) +
    geom_vline(aes(xintercept=1)) +
    labs(title = paste(title,': All (',length(results$VE.all),') simulated trials', sep = "")) +
    geom_vline(aes(xintercept=truth$PoDModel$valueParams$param.pmax.covariateBinary),
               linetype = "dashed", 
               colour = "gray") +
    geom_vline(aes(xintercept=truth$PoDModel$valueParams$param.gamma.covariateBinary),
               linetype = "dashed", 
               colour = "gray") +
    geom_vline(aes(xintercept=truth$PoDModel$valueParams$param.et50.covariateBinary),
               linetype = "dashed", 
               colour = "gray") +
    xlim(0,2.5)+
    annotate("text", x = 2.3, y = 1.5, label = paste("# trials = ",sum(!is.na(results$PoD.param.et50)),sep = ""))+
    annotate("text", x = 2.3, y = 2.5, label = paste("# trials = ",sum(!is.na(results$PoD.param.et50.covariateBinary)),sep = ""))+
    annotate("text", x = 2.3, y = 3.5, label = paste("# trials = ",sum(!is.na(results$PoD.param.gamma)),sep = ""))+
    annotate("text", x = 2.3, y = 4.5, label = paste("# trials = ",sum(!is.na(results$PoD.param.gamma.covariateBinary)),sep = ""))+
    annotate("text", x = 2.3, y = 5.5, label = paste("# trials = ",sum(!is.na(results$PoD.param.pmax)),sep = ""))+
    annotate("text", x = 2.3, y = 6.5, label = paste("# trials = ",sum(!is.na(results$PoD.param.pmax.covariateBinary)),sep = ""))+
    annotate("text", x = 2.3, y = 7.5, label = paste("# trials = ",sum(!is.na(results$VE.all)),sep = ""))
  
  
}

ridgelineOneSimulatedTrial <- function(id=1, results, truth, trueVE, seeds, title){
  idDetail <- which(seeds == id)
  VE <- results$VECI_raw[[idDetail]]$VECI_raw$all/trueVE
  pmax.covariate <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.pmax.covariateBinary
  pmax <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.pmax/truth$PoDModel$valueParams$param.pmax
  gamma.covariate <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.gamma.covariateBinary
  gamma <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.gamma/truth$PoDModel$valueParams$param.gamma
  et50.covariate <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.et50.covariateBinary
  et50 <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.et50/truth$PoDModel$valueParams$param.et50
  
  val <- c(VE,
           pmax.covariate,
           pmax,
           gamma.covariate,
           gamma,
           et50.covariate,
           et50)
  name <- c(rep("VE", length(VE)),
            rep("pmax_ageGroup", length(pmax.covariate)),
            rep("pmax", length(pmax)),
            rep("gamma_ageGroup", length(gamma.covariate)),
            rep("gamma", length(gamma)),
            rep("et50_ageGroup", length(et50.covariate)),
            rep("et50", length(et50)))
  
  df <- data.frame(value = val,
                   variable = name)
  
  ggplot(df,
         aes(x = value, y = variable, fill = stat(abs(1-x)))) +
    geom_density_ridges_gradient(scale = 1.5, alpha=0.5) +
    scale_fill_viridis(name = "value", option = "B") +
    theme_ridges(grid = FALSE) + 
    theme(legend.position = "none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 8)) +
    geom_vline(aes(xintercept=1)) +
    labs(title = paste(title,': One simulated trial', sep = "")) +
    geom_vline(aes(xintercept=truth$PoDModel$valueParams$param.pmax.covariateBinary),
               linetype = "dashed", 
               color = "gray") +
    geom_vline(aes(xintercept=truth$PoDModel$valueParams$param.gamma.covariateBinary),
               linetype = "dashed", 
               colour = "gray") +
    geom_vline(aes(xintercept=truth$PoDModel$valueParams$param.et50.covariateBinary),
               linetype = "dashed",
               colour = "gray") +
    xlim(0,2.5)
}

ridgelineQuartiles <- function (id=1, seeds, results){
  idDetail <- which(seeds == id)
  pmax.covariate <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.pmax.covariateBinary

  gamma.covariate <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.gamma.covariateBinary

  et50.covariate <- results$PoDCI_raw[[idDetail]]$PoDCI_raw$results$param.et50.covariateBinary


  val <- c(pmax.covariate,
           gamma.covariate,
           et50.covariate)
  name <- c(rep("pmax_ageGroup", length(pmax.covariate)),
            rep("gamma_ageGroup", length(gamma.covariate)),
            rep("et50_ageGroup", length(et50.covariate)))

  df <- data.frame(value = val,
                   variable = name)

  # Highlight the tails of the distributions
  ggplot(df, aes(x = value, y = variable)) +
    geom_density_ridges(
      jittered_points = TRUE, quantile_lines = TRUE, scale = 0.7, alpha = 0.7,
      vline_size = 1, vline_color = "gray20",
      point_size = 0.4, point_alpha = 1,
      position = position_raincloud(adjust_vlines = TRUE)
    )+
    theme_ridges(grid = FALSE) + 
    theme(legend.position = "none",
          panel.spacing = unit(0.1, "lines"),
          strip.text.x = element_text(size = 8)) +
    geom_vline(aes(xintercept=1)) +
    xlim(0,2.5)

}


## Truth definition for the illustrative example in the manuscript
defineTruthExample <- function() {
  name <- "Illustrative_example"
  
  variable <- list()
  parameter <- list()
  vaccinated <- list()
  control <- list()
  
    ## PoD model for data generation ----
  PoD <- "I(ifelse(var.titer > 0, param.pmax*(param.pmax.covariateBinary^var.covariateBinary) /(
  1+(var.titer/(param.et50*(param.et50.covariateBinary^var.covariateBinary))
  )^(param.gamma*(param.gamma.covariateBinary^var.covariateBinary) )), param.pmax))"
  
  # Independent variables
  variable$name <- c("var.titer", "var.covariateBinary")
  
  # Parameters
  parameter$name <- c("param.pmax",
                      "param.pmax.covariateBinary",
                      "param.et50",
                      "param.et50.covariateBinary",
                      "param.gamma",
                      "param.gamma.covariateBinary")
  parameter$value <- list(param.pmax = 0.033,
                          param.pmax.covariateBinary = 1,
                          param.et50 = 7.12,
                          param.et50.covariateBinary = 1.361, #1.131; 1.361 
                          param.gamma = 7,
                          param.gamma.covariateBinary = 1)
  
  ## Populations ----
  
  # Vaccinated 
  
  # titer distribution
  vaccinated$N = 10000
  vaccinated$mean = 10
  vaccinated$stdDev = 2
  
  # covariate information
  # the first value in the interval will be adjusted in the code to 1, the other into 0 
  vaccinated$covariateInfo$var.covariateBinary <- list("name" = "var.covariateBinary",
                                                       "values" = c("valueA","valueB"),
                                                       "prob" = c(0.25, 0.75))
  
  # Control
  # titer distribution
  control$N = 5000
  control$mean = 5
  control$stdDev = 2
  
  # covariate information
  control$covariateInfo$var.covariateBinary <- list("name" = "var.covariateBinary",
                                                    "values" = c("valueA","valueB"),
                                                    "prob" = c(0.25, 0.75))
  
  
  ## Immunogenicity subsetting methods: "Full", "Fixed" ----
  method <- list(name = "Full",
                 value = NA)

  # method <- list(name = "Fixed",
  #                value = 1600)
  
  PoDModel <- MakeModelObj(PoD,
                           variable$name,
                           variable$range,
                           parameter$name,
                           parameter$range,
                           parameter$value)
  
  return(list(
    name = name,
    vaccinated = vaccinated,
    control = control,
    PoDModel = PoDModel,
    method = method,
    adjustTiters = FALSE, # this titer adjustment improves convergence:
    adjustFrom = 0, # titers shouldn't be negative, otherwise the negLL fails sometimes
    adjustTo = 0
  ))
}

ConvertNormalToGeometric <- function(mean,
                                     sd,
                                     n,
                                     ci = 0.95){
  se <- sd/sqrt(n) 
  GMT.upperBound <- 2^(mean+(se * qt(ci + (1 - ci)/2, df = n - 1)))
  GMT.lowerBound <- 2^(mean-(se * qt(ci + (1 - ci)/2, df = n - 1)))
  GMT <- 2^mean
  
  return(list(GMT = GMT,
              GMT.upperBound = GMT.upperBound,
              GMT.lowerBound = GMT.lowerBound))
}
