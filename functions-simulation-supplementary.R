# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

getPopulationGroup <- function(population1, 
                               attribute,
                               group,
                               fromDataDist = FALSE) {
  
  new <- generatePopulation(0)
  
  if (length(population1$PoDs)>0 & length(population1$diseaseStatus)>0){
    df <- cbind(population1$titers, population1$diseaseStatus, population1$PoDs, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.diseaseStatus"
    colnames(df)[3] <- "y.PoD"
    newDf <- filter(df, UQ(sym(attribute)) == group)
    newTiters <- newDf$var.titer
    newDf$var.titer <- NULL  
    newDS <- newDf$y.diseaseStatus
    newDf$y.diseaseStatus <- NULL
    newPoDs <- newDf$y.PoD
    newDf$y.PoD <- NULL
    newAttributes <- newDf
    
    new$PoDs <- newPoDs
    new$diseaseStatus <- newDS
    
  } else if (length(population1$diseaseStatus)>0) {
    df <- cbind(population1$titers, population1$diseaseStatus, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.diseaseStatus"
    newDf <- filter(df, UQ(sym(attribute)) == group)
    newTiters <- newDf$var.titer
    newDf$var.titer <- NULL  
    newDSs <- newDf$y.diseaseStatus
    newDf$y.diseaseStatus <- NULL
    newAttributes <- newDf
    
    new$diseaseStatus <- newDS
    
  } else if (length(population1$PoDs)>0) {
    df <- cbind(population1$titers, population1$PoDs, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.PoD"
    newDf <- filter(df, UQ(sym(attribute)) == group)
    newTiters <- newDf$var.titer
    newDf$var.titer <- NULL  
    newPoDs <- newDf$y.PoD
    newDf$y.PoD <- NULL
    newAttributes <- newDf
    
    new$PoDs <- newPoDs
    
  } else {
    df <- cbind(population1$titers, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    newDf <- filter(df, UQ(sym(attribute)) == group)
    newTiters <- newDf$var.titer
    newDf$var.titer <- NULL  
    newAttributes <- newDf
  }
  
  new$titers <- newTiters
  new$N <- length(new$titers)
  new$mean <- population1$mean
  new$stdDev <- population1$stdDev
  new$subjectAttributes <-newAttributes
  
  if (fromDataDist) {
    new$mean <- mean(newTiters)
    new$stdDev <- sd(newTiters)
  }
  
  return(new)
}

populationToDF <- function(vaccPOP, contPOP) {
  population1 <- mergePopulations(vaccPOP, contPOP)
  if (length(population1$PoDs)>0 & length(population1$diseaseStatus)>0){
    df <- cbind(population1$titers, as.numeric(population1$diseaseStatus), population1$PoDs, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.diseaseStatus"
    colnames(df)[3] <- "y.PoD"
    
  } else if (length(population1$diseaseStatus)>0) {
    df <- cbind(population1$titers, as.numeric(population1$diseaseStatus), population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.diseaseStatus"
    
  } else if (length(population1$PoDs)>0) {
    df <- cbind(population1$titers, population1$PoDs, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.PoD"
    
  } else {
    df <- cbind(population1$titers, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
  }
  
  return(df)
}

#' @title Clinical trial: estimation of case-count efficacy
#'
#' @description Function assigns disease status (DS) to vaccinated and control groups and based on that calculates the case-count efficacy. Vaccinated and control groups are provided in the form of population class objects (see the \code{Population-class} function for more details).
#'
#' Input populations need to contain information about Probability of disease (PoD) for each subject - calculated using \code{population$assignPoD(PoD(x))}. See \code{PoD} function for further details.
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects with assigned PoD
#' @param control \code{Population-class} object: control subjects with assigned PoD
#' @param CI numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return
#' \itemize{
#'   \item vaccinated: vaccinated subjects with assigned DS, \code{Population-class} object
#'
#'   \item control: control subjects with assigned DS, \code{Population-class} object
#'
#'   \item efficacy: case-count efficacy
#'
#'   \item confidenceInterval: case-count efficacy confidence interval calculated with \code{waldCI()} function
#' }
#'
#' @usage
#' ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' @examples
#' # Loading vaccinated, control population data with PoD information
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with 95\% confidence interval
#' CT <- ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' CT$efficacy
#' CT$confidenceInterval
#'
#' CT$vaccinated
#'
#' @export
ClinicalTrial <- function(vaccinated,
                          control,
                          CI = 0.95){
  
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  dsVacc <- mapply(rbinom, n = 1, size = 1, prob = vaccinated$PoDs)
  vaccinated$diseaseStatus <- sapply(dsVacc, numToBool)
  dsCont <- mapply(rbinom, n = 1, size = 1, prob = control$PoDs)
  control$diseaseStatus <- sapply(dsCont, numToBool)
  
  casecountEfficacy <- 1 - (vaccinated$getDiseasedCount() / vaccinated$N) /
    (control$getDiseasedCount() / control$N)
  
  confidenceInterval <- waldCI(vaccinated, control, CI)
  
  return(list(
    vaccinated = vaccinated,
    control = control,
    efficacy = casecountEfficacy,
    confidenceInterval = confidenceInterval
  ))
}

#' @title Clinical trial function expanded for usage in simulations when the calculation of coverage probability is needed for three confidence intervals: 80\%, 90\%, and user-defined  
#' 
#' 
#' @description Function works the same way as \code{ClinicalTrial} function but it also calculates 80\% and 90\% confidence intervals.
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects with assigned PoD
#' @param control \code{Population-class} object: control subjects with assigned PoD
#' @param CI numeric: value from (0, 1) interval, confidence level of interest 
#'
#' @return
#' \itemize{
#'   \item vaccinated: vaccinated subjects with assigned DS, \code{Population-class} object
#'
#'   \item control: control subjects with assigned DS, \code{Population-class} object
#'
#'   \item efficacy: case-count efficacy
#'
#'   \item confidenceInterval: confidence interval calculated with \code{waldCI} function
#'
#'   \item confidenceInterval90: 90\% confidence interval calculated with \code{waldCI} function
#'
#'   \item confidenceInterval80: 80\% confidence interval calculated with \code{waldCI} function
#' }
#'
#' @usage
#' ClinicalTrialCoverage(vaccinated, control, CI = 0.95)
#'
#' @export
ClinicalTrialCoverage <- function(vaccinated,
                                  control,
                                  CI = 0.95){
  
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  dsVacc <- mapply(rbinom, n = 1, size = 1, prob = vaccinated$PoDs)
  vaccinated$diseaseStatus <- sapply(dsVacc, numToBool)
  dsCont <- mapply(rbinom, n = 1, size = 1, prob = control$PoDs)
  control$diseaseStatus <- sapply(dsCont, numToBool)
  
  casecountEfficacy <- 1 - (vaccinated$getDiseasedCount() / vaccinated$N) /
    (control$getDiseasedCount() / control$N)
  
  confidenceInterval95 <- waldCI(vaccinated, control, CI)
  confidenceInterval90 <- waldCI(vaccinated, control, 0.90)
  confidenceInterval80 <- waldCI(vaccinated, control, 0.80)
  
  return(list(
    vaccinated = vaccinated,
    control = control,
    efficacy = casecountEfficacy,
    confidenceInterval95 = confidenceInterval95,
    confidenceInterval90 = confidenceInterval90,
    confidenceInterval80 = confidenceInterval80
  ))
}

#' @title Wald confidence interval estimation
#'
#' @description Function calculates and returns case-count efficacy confidence intervals estimated using Wald's method.
#'
#' Input data need to contain information about disease status on individual level.
#'
#' @param vaccinated \code{Population-class} object: vaccinated subjects, containing information about disease status
#' @param control \code{Population-class} object: control subjects, containing information about disease status
#' @param confLevel numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return Named list of lower and upper confidence interval bound
#'
#' @details  Confidence interval of the relative risk is calculated using the Wald method. (Wald, A. Tests of statistical hypotheses concerning several parameters when the number of observations is large. Transactions of the American Mathematical Society 54, 426-482 (1943)).
#'
#' @examples
#' # Loading vaccinated and control populations data with PoD information
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with 95\% confidence interval
#' set.seed(1)
#' CT <- ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' waldCI(vaccinated, control)
#'
#' @export
waldCI <- function(vaccinated, 
                   control, 
                   confLevel = 0.95){
  
  if (!is(vaccinated, "Population")) {
    incorrectPopulationInput("vaccinated")
  }
  if (!is(control, "Population")) {
    incorrectPopulationInput("control")
  }
  
  # significance level z value
  N. <- 1 - ((1 - confLevel)/2)
  z <- stats::qnorm(N., mean = 0, sd = 1)
  
  # Get matrix of vaccine/control + diseased/non-diseased 
  A <- vaccinated$getDiseasedCount()
  B <- vaccinated$getNondiseasedCount()
  C <- control$getDiseasedCount()
  D <- control$getNondiseasedCount()
  
  # Get size of vaccinated and control
  N_VACC <- A + B # N VACC
  N_CONTROL <- C + D # N CONT
  
  # calculate RR point estimate
  RR_point_t <- (A/N_VACC) / (C/N_CONTROL)
  ln_RR_point <- log(RR_point_t)
  
  # Wald RR se  
  ln_RR_point_se <- sqrt((1/A) - (1/N_VACC) + (1/C) - (1/N_CONTROL))
  RR_point_se <- exp(ln_RR_point_se)
  
  # Wald RR confidence intervals
  RR_lower <- 1 - exp(ln_RR_point + (z * ln_RR_point_se))
  RR_upper <- 1 - exp(ln_RR_point - (z * ln_RR_point_se))
  RR_point <- 1 - RR_point_t
  
  # summary
  wald_ci_point <- RR_point
  wald_ci_lower <- RR_lower
  wald_ci_upper <- RR_upper
  
  return(list(
    lowerBound = wald_ci_lower,
    upperBound = wald_ci_upper
  )
  )
}

#' @title Immunogenicity subset: vaccinated, control, non-diseased
#'
#' @description
#' Function creates non-diseased immunogenicity subset, and vaccinated and control immunogenicity subsets based on chosen method. The immunogenicity subsets are provided in the form of population class objects (see the \code{Population-class} function for more details).
#' Works with covariate data about the population (stored in subjectAttributes field of the population).
#'
#' @param diseased \code{Population-class} object: diseased subjects, created using \code{ExtractDiseased} function; should have data in the field subjectAttributes
#' @param nondiseased \code{Population-class} object: non-diseased subjects, created using \code{ExtractNondiseased} function; should have data in the field subjectAttributes
#' @param method named list: "name" possible inputs "Full", "Fixed";
#'
#' "value" = numeric value
#'
#' @return
#' \itemize{
#'   \item ImmunogenicityVaccinated: vaccinated subjects in the immunogenicity subset, \code{Population-class} object (N, mean, stdDev, titers)
#'
#'   \item ImmunogenicityControl: control subjects in the immunogenicity subset, \code{Population-class} object (N, mean, stdDev, titers)
#'
#'   \item ImmunogenicityNondiseased: non-diseased subjects in the immunogenicity subset, \code{Population-class} object (N, mean, stdDev, titers)
#' }
#'
#' @usage
#' BlindSampling(diseased, 
#'               nondiseased,  
#'               method = list(name = "Full", value = NA))
#'
#' @details
#'
#' For details about the method parameter see \code{ImmunogenicitySubset} function.
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' ## Example 1
#' # Creating immunogenicity subset, method = "Full"
#' ImmunogenicitySubsetFull <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Full", 
#'                                 value = NA))
#'
#' ## Example 2
#' # Creating of immunogenicity subset, method = "Fixed"
#' ImmunogenicitySubsetFixed <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Fixed", 
#'                                 value = 100))
#'
#' @importFrom dplyr filter
#' @export
BlindSampling <- function(diseased, nondiseased, method = list(name = "Full",
                                                               value = NA)) {
  
  if (!is(diseased, "Population")) {
    incorrectPopulationInput("diseased")
  }
  if (!is(nondiseased, "Population")) {
    incorrectPopulationInput("nondiseased")
  }
  
  # Create immunogenicity population
  immunogenicitySample <- ImmunogenicitySubset(diseased = diseased,
                                               nondiseased = nondiseased,
                                               method = method)
  
  # Create a dataframe with titers, DS and attributes
  df <- cbind(immunogenicitySample$titers, immunogenicitySample$diseaseStatus, immunogenicitySample$PoDs, immunogenicitySample$subjectAttributes)
  colnames(df)[1] <- "var.titer"
  colnames(df)[2] <- "y.diseaseStatus"
  colnames(df)[3] <- "y.PoDs"
  
  # Vaccinated df and immunogenicity population
  dfVaccinated <- filter(df,var.vaccStatus=="vaccinated")
  immunogenicityVaccinatedTiters <- dfVaccinated$var.titer
  immunogenicityVaccinatedDS <- dfVaccinated$y.diseaseStatus
  immunogenicityVaccinatedPoDs <- dfVaccinated$y.PoDs
  dfVaccinated$var.titer <- NULL
  dfVaccinated$y.diseaseStatus <- NULL
  dfVaccinated$y.PoDs <- NULL
  immunogenicityVaccinatedAttributes <- dfVaccinated
  
  immunogenicityVaccinated                   <- generatePopulation(0)
  immunogenicityVaccinated$N                 <- length(immunogenicityVaccinatedTiters)
  immunogenicityVaccinated$stdDev            <- sd(immunogenicityVaccinatedTiters)
  immunogenicityVaccinated$mean              <- mean(immunogenicityVaccinatedTiters)
  immunogenicityVaccinated$titers            <- immunogenicityVaccinatedTiters
  immunogenicityVaccinated$diseaseStatus     <- immunogenicityVaccinatedDS
  immunogenicityVaccinated$PoDs              <- immunogenicityVaccinatedPoDs
  immunogenicityVaccinated$subjectAttributes <- immunogenicityVaccinatedAttributes
  
  # Control df and immunogenicity population
  dfControl <- filter(df,var.vaccStatus=="control")
  immunogenicityControlTiters <- dfControl$var.titer
  immunogenicityControlDS <- dfControl$y.diseaseStatus
  immunogenicityControlPoDs <- dfControl$y.PoDs
  dfControl$var.titer <- NULL
  dfControl$y.diseaseStatus <- NULL
  dfControl$y.PoDs <- NULL
  immunogenicityControlAttributes <- dfControl
  
  immunogenicityControl                   <- generatePopulation(0)
  immunogenicityControl$N                 <- length(immunogenicityControlTiters)
  immunogenicityControl$stdDev            <- sd(immunogenicityControlTiters)
  immunogenicityControl$mean              <- mean(immunogenicityControlTiters)
  immunogenicityControl$titers            <- immunogenicityControlTiters
  immunogenicityControl$diseaseStatus     <- immunogenicityControlDS
  immunogenicityControl$PoDs              <- immunogenicityControlPoDs
  immunogenicityControl$subjectAttributes <- immunogenicityControlAttributes
  
  # Non-diseased immunogenicity population
  immunogenicityNondiseased <- ExtractNondiseased(immunogenicityVaccinated, immunogenicityControl)
  
  return(list(
    immunogenicityNondiseased   = immunogenicityNondiseased,
    immunogenicityVaccinated    = immunogenicityVaccinated,
    immunogenicityControl       = immunogenicityControl
  ))
}

#' @title Immunogenicity subset
#'
#' @description
#' Function creates the immunogenicity subset based on the chosen method.
#' Works with covariate data about the population (stored in subjectAttributes field of the population).
#'
#' @param diseased \code{Population-class} object: diseased subjects with assigned vaccination status; should have data in the field subjectAttributes
#' @param nondiseased \code{Population-class} object: non-diseased subjects with assigned vacination status; should have data in the field subjectAttributes
#' @param method named list: a selected method for creating the immunogenicity subset
#'
#' method$name
#'
#' \itemize{
#'   \item Full: subject level titer information is available for all diseased and all non-diseased subjects, i.e. immunogenicity subset is the full clinical trial
#'
#'   \item Fixed: subject level titer information is available for all diseased and some non-diseased subjects.
#' }
#'
#'
#' method$value
#'
#' \itemize{
#'   \item Full: value = NA; immunogenicity sample is the full clinical trial (non-diseased subset contains all non-diseased in the trial; diseased subset contains all disease cases in the trial)
#'
#'   \item Fixed: value = size of the immunogenicity subset, pre-defined number of subjects assayed for titers independently of their future disease status (non-diseased subset could rarely contain some diseased subjects, as the selection is done at the enrollment and prior the knowledge of future disease status; diseased subset contains all disease cases in the trial)
#' }
#'
#' @return
#' Immunogenicity subset with subject level information about vaccination status and disease status, provided in the form of \code{Population-class} object
#'
#' @usage
#' ImmunogenicitySubset(diseased, 
#'                      nondiseased, 
#'                      method = list(name = "Full", value = NA))
#'
#' @details
#' The total immunogenicity subset consists of the diseased immunogenicity subset and non-diseased immunogenicity subset. 
#' For all three methods implemented, we assume that the diseased immunogenicity subset contains all disease cases in the trial.
#' Based on the chosen method, the the size of the non-diseaded immunogenicity subset can be derived as follows:
#'
#' Size = number of subjects in the non-diseased immunogenicity subset
#'
#' Titers = values of titers from which we want to sample in order to simulate the non-diseased immunogenicity subset
#'
#' #Diseased = total number of diseased in the clinical trial
#'
#' #Nondiseased = total number of non-diseased in the clinical trial
#'
#'  \itemize{
#'   \item method$name = "Full"
#'   
#'   Size = #Nondiseased
#'   
#'   Titers = Nondiseased Titers 
#'
#'   \item method$name = "Fixed"
#'   
#'   Size = method$value
#'   
#'   Titers = Nondiseased Titers + Diseased Titers
#'   
#' }
#'
#' @examples
#' ## Example 1
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' ImmunogenicitySubset(diseased,
#'                      nondiseased,
#'                      method = list(name = "Fixed",
#'                                    value = 150))
#'
#' @export
ImmunogenicitySubset <- function(diseased,
                                 nondiseased,
                                 method = list(name = "Full",
                                               value = NA)) {
  if (!is(diseased, "Population")) {
    incorrectPopulationInput("diseased")
  }
  if (!is(nondiseased, "Population")) {
    incorrectPopulationInput("nondiseased")
  }
  
  if(any(is.na(match(names(method), c("name", "value") )))){
    stop(paste("The input value for method is incorrect. 'method' parameter has wrong names."))
  }
  
  if (is.na(
    match(method$name, c("Full", "Fixed")))
  ) { stop(paste("The input value for method$name is incorrect. 'method$name' needs to be either \"Full\" or \"Fixed\".")) }
  
  if (method$name == "Full") {
    forImmunogenicitySample <- mergePopulations(diseased, nondiseased)
    
  } else if (method$name == "Fixed") {
    immunogenicitySize <- method$value
    if (!is.numeric(immunogenicitySize) | immunogenicitySize <= 0 | is.na(immunogenicitySize)) { 
      stop(paste("The input value for method$value is incorrect. 'method$value' needs to be positive numeric value."))
    }
    forImmunogenicitySample <- mergePopulations(diseased, nondiseased)
    forImmunogenicitySample$N <- immunogenicitySize
  }
  
  immunogenicitySample <- samplePopulation(forImmunogenicitySample,
                                           round(forImmunogenicitySample$N),
                                           replace = FALSE)
  
  return(immunogenicitySample)
}

#' @title PoDBAY efficacy estimation; works with covariate data.
#'
#' @description 
#' Function calculates the PoDBAY efficacy based on the set of PoD curve 
#' parameters calculated in \code{PoDParamEstimation} function, vaccinated and 
#' control immunogenicity subset means and standard deviations.
#' 
#' Uses a weighted PoD curve along with user defined titer distributions
#' to calculate the CoP-based vaccine efficacy.
#'
#' @param podEstimatedParameters list of estimated PoD curve parameters
#' @param pdfSampledParametersCont list of sampled PDF parameters of control group
#' @param pdfSampledParametersVacc list of sampled PDF parameters of vaccinated group
#' @param vaccPOP \code{Population-class} object: vaccinated population 
#' @param contPOP \code{Population-class} object: control population
#' @param covInfo user supplied covariate information
#' @param podModelObject model object from \code{MakeModelObj} function that includes PoD function
#' @param pdfContModelObject model object from \code{MakeModelObj} function for control PDF
#' @param pdfVaccModelObject model object from \code{MakeModelObj} function for vaccinated PDF
#' @param plotDistTF set to true and return plot of distributions
#' @param fromDataTF boolean
#' @param integrationLHL endpoint of the integration interval
#' @param integrationRHL endpoint of the integration interval
#' @param integrationTOL numeric: integration tolerance
#'
#' @return 
#' efficacySet, set of PoDBAY efficacies corresponding to estimated set of PoD curve parameters
#'
#' @examples
#' ## Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(PoDModel)
#' data(covariateInfo)
#'
#' ## Example 1
#' # Creating imunogenicity subset, method = "Fixed", value = 150
#' ImmunogenicitySample <- 
#'   BlindSampling(diseased, 
#'                 nondiseased, 
#'                 method = list(name = "Fixed", 
#'                               value = 150))
#'                               
#' # Estimating PoD curve parameters
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' estimatedParameters <- PoDParamEstimation(diseased,
#'                        ImmunogenicitySample$immunogenicityNondiseased,
#'                        nondiseasedGenerationCount,
#'                        PoDModel,
#'                        repeatCount = 3)
#'                        
#' # Estimating PoDBAY efficacy  
#' PoDBAYEfficacy_MV(podEstimatedParameters = estimatedParameters$results,
#'                   vaccPOP = ImmunogenicitySample$immunogenicityVaccinated,
#'                   contPOP = ImmunogenicitySample$immunogenicityControl,
#'                   covInfo = covariateInfo, 
#'                   podModelObject = PoDModel,
#'                   fromDataTF = TRUE)
#'
#' @export
PoDBAYEfficacy_MV <- function(podEstimatedParameters,
                              pdfSampledParametersCont = NULL,
                              pdfSampledParametersVacc = NULL,
                              vaccPOP,
                              contPOP,
                              covInfo = NULL,
                              podModelObject,
                              pdfContModelObject = NULL,
                              pdfVaccModelObject = NULL,
                              plotDistTF = FALSE,
                              fromDataTF = FALSE,
                              integrationLHL = -100,
                              integrationRHL = 100,
                              integrationTOL = 1e-8) {
  
  if (!is(vaccPOP, "Population")) {
    incorrectPopulationInput("vaccPOP")
  }
  if (!is(contPOP, "Population")) {
    incorrectPopulationInput("contPOP")
  }
  efficacySet <- numeric()
  
  # Each iteration must reset the titerDist
  #   formals are altered in loops below
  if (is.null(pdfSampledParametersVacc) & is.null(pdfVaccModelObject)) {
    if (is.null(pdfSampledParametersCont) & is.null(pdfContModelObject)) {
      for (i in 1:nrow(podEstimatedParameters)) {
        print(i)
        jitterVaccinated <- JitterMean(vaccPOP)
        jitterControl <- JitterMean(contPOP)
        means <- list(
          vaccinated = jitterVaccinated,
          control = jitterControl
        )
        stdDevs <- list(
          vaccinated = vaccPOP$stdDev,
          control = contPOP$stdDev
        )
        
        # Define temp/utility populations for use in efficacy computation
        vaccPOPUtil <- vaccPOP$copy()
        vaccPOPUtil[["mean"]] <- means[["vaccinated"]]
        vaccPOPUtil[["stdDev"]] <- stdDevs[["vaccinated"]]
        contPOPUtil <- contPOP$copy()
        contPOPUtil[["mean"]] <- means[["control"]]
        contPOPUtil[["stdDev"]] <- stdDevs[["control"]]
        
        efficacy <- EfficacyComputation_MV(vaccPOP = vaccPOPUtil,
                                           contPOP = contPOPUtil,
                                           covInfo = covInfo,
                                           paramIn = podEstimatedParameters[i,],
                                           podModelObject = podModelObject,
                                           integrationLHL = integrationLHL,
                                           integrationRHL = integrationRHL,
                                           integrationTOL = integrationTOL,
                                           fromDataTF = fromDataTF)
        
        efficacySet <- c(efficacySet, efficacy)
        
        # # For testing
        # distDF.util <- data.frame(vaccMean = means$vaccinated,
        #                           vaccStd = stdDevs$vaccinated,
        #                           contMean = means$control,
        #                           contStd = stdDevs$control)
        # if (!exists("distDF")){
        #   distDF <- distDF.util
        # } else {
        #   distDF <- rbind(distDF,distDF.util)
        # }
      }
    } 
  } else {
    vaccPOP$userDefinedTiterDist <- TRUE
    vaccPOP$titerDist <- pdfVaccModelObject$predictionFun
    for (paramName in pdfVaccModelObject$params) {
      # Change default argument values to sampled parameter values so that 
      #   the function in "titerDist" field is only dependent on titer, and
      #   all other parameters are those found from sampling parameter distribution.
      formalsStringUtil <- paste0("formals(vaccPOP$titerDist)['",
                                  paramName,
                                  "'] <- pdfSampledParametersVacc[",
                                  i,
                                  ",][['",
                                  paramName,
                                  "']]")
      eval(parse(text=formalsStringUtil))
    }
    contPOP$userDefinedTiterDist <- TRUE
    contPOP$titerDist <- pdfContModelObject$predictionFun
    for (paramName in pdfContModelObject$params) {
      formalsStringUtil <- paste0("formals(contPOP$titerDist)['",
                                  paramName,
                                  "'] <- pdfSampledParametersCont[",
                                  i,
                                  ",][['",
                                  paramName,
                                  "']]")
      eval(parse(text=formalsStringUtil))
    }
    efficacy <- EfficacyComputation_MV(vaccPOP = vaccPOP,
                                       contPOP = contPOP,
                                       covInfo = covInfo,
                                       paramIn = podEstimatedParameters[i,],
                                       podModelObject = podModelObject,
                                       integrationLHL = integrationLHL,
                                       integrationRHL = integrationRHL,
                                       integrationTOL = integrationTOL,
                                       fromDataTF = fromDataTF)
    efficacySet <- c(efficacySet, efficacy)
  }
  return(efficacySet = efficacySet)
}

#' @title PoDBAY efficacy equation; works with covariate data.
#'
#' @description
#' Function calculates the PoDBAY efficacy based on the PoD curve parameters and titer distribution 
#' parameters (mean, sd) for vaccinated and control groups.
#'
#' @param vaccPOP \code{Population-class} object: vaccinated population 
#' @param contPOP \code{Population-class} object: control population
#' @param covInfo user supplied covariate information
#' @param paramIn PoD curve parameters
#' @param podModelObject model object from \code{MakeModelObj} function that includes PoD function
#' @param integrationLHL endpoint of the integration interval
#' @param integrationRHL endpoint of the integration interval
#' @param integrationTOL relative tolerance for integration
#' @param fromDataTF boolean
#'
#' @return
#' efficacy: numeric value
#'
#' @usage
#' EfficacyComputation_MV(vaccPOP,
#'                        contPOP,
#'                        covInfo = NULL,
#'                        paramIn,
#'                        podModelObject,
#'                        integrationLHL = -100,
#'                        integrationRHL = 100,
#'                        integrationTOL = 1e-8,
#'                        fromDataTF = FALSE)
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#' data(control)
#' data(PoDModel)
#' data(covariateInfo)
#'
#' ## Example 1 
#' PoDParams <- list("param.pmax" = 0.05,
#'                   "param.pmax.covariateBinary" = 1,
#'                   "param.et50" = 5,
#'                   "param.et50.covariateBinary" = 2,
#'                   "param.gamma" = 7,
#'                   "param.gamma.covariateBinary" = 1)
#' 
#' EfficacyComputation_MV(vaccPOP = vaccinated, 
#'                        contPOP = control,
#'                        covInfo = covariateInfo, # NULL does not work now
#'                        paramIn = PoDParams,
#'                        podModelObject = PoDModel,
#'                        fromDataTF = TRUE)
#' @export
EfficacyComputation_MV <- function(vaccPOP,
                                   contPOP,
                                   covInfo = NULL,
                                   paramIn,
                                   podModelObject,
                                   integrationLHL = -100,
                                   integrationRHL = 100,
                                   integrationTOL = 1e-8,
                                   fromDataTF = FALSE) {
  
  # Expected probability of disease vacc
  ExpectedPoDVacc <- ExpectedPoD_MV(trialPOP = vaccPOP,
                                    vaccStat = 1,
                                    covInfo = covInfo,
                                    paramForModelObject = paramIn,
                                    modelObject = podModelObject,
                                    fromDataTF = fromDataTF, # if true data is taken from trialPOP
                                    integrationLHL = integrationLHL,
                                    integrationRHL = integrationRHL,
                                    integrationTOL = integrationTOL)
  # Expected probability of disease control
  ExpectedPoDCont <- ExpectedPoD_MV(trialPOP = contPOP,
                                    vaccStat = 0,
                                    covInfo = covInfo,
                                    paramForModelObject = paramIn,
                                    modelObject = podModelObject,
                                    fromDataTF = fromDataTF, # if true data is taken from trialPOP
                                    integrationLHL = integrationLHL,
                                    integrationRHL = integrationRHL,
                                    integrationTOL = integrationTOL)
  efficacy <- 1-ExpectedPoDVacc/ExpectedPoDCont
  return(efficacy)
}


#' @title Expected probability of disease; works with covariate data. 
#'
#' @description Function returns expected probability of disease for a trial population given the 
#   PoD curve model, PoD curve parameters, and
#   the covariate information (either supplied in covInfo or
#   calculated from the trialPOP$subjectAttributes)
#'
#' @param trialPOP \code{Population-class} object: trial population (vaccinated or control)
#' @param vaccStat vaccination status; 1 for vaccinated population, 0 for control population
#' @param covInfo user supplied covariate information
#' @param paramForModelObject PoD curve parameters
#' @param modelObject model object from \code{MakeModelObj} function that includes PoD function
#' @param fromDataTF boolean: set to TRUE if data is taken from trialPOP
#' @param integrationLHL endpoint of the integration interval
#' @param integrationRHL endpoint of the integration interval
#' @param integrationTOL relative tolerance for integration
#'
#' @return 
#' Result of numerical integration of the PoD curve and population titer distribution
#'
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#' data(covariateInfo)
#' data(PoDModel)
#' 
#' ## Example 1
#' PoDParams <- list("param.pmax" = 0.05,
#'                   "param.pmax.covariateBinary" = 1,
#'                   "param.et50" = 5,
#'                   "param.et50.covariateBinary" = 2,
#'                   "param.gamma" = 7,
#'                   "param.gamma.covariateBinary" = 1)
#' 
#' ExpectedPoD_MV(vaccinated,
#'                vaccStat = 1,
#'                covariateInfo,
#'                PoDParams,
#'                PoDModel,
#'                fromDataTF = TRUE, 
#'                integrationLHL = -100,
#'                integrationRHL = 100,
#'                integrationTOL = 1e-8)
#'                
#' @importFrom pracma quadgk               
#' @export
ExpectedPoD_MV <- function(trialPOP,
                           vaccStat = NULL,
                           covInfo,
                           paramForModelObject,
                           modelObject,
                           fromDataTF = FALSE, # if true data is taken from trialPOP
                           integrationLHL = -100,
                           integrationRHL = 100,
                           integrationTOL = 1e-8) {
  
  # Assign the data from the population class
  if (fromDataTF==TRUE) {
    dataInput <- trialPOP$subjectAttributes
  } else {
    dataInput <- NULL
  }
  
  # Reduce to one argument for use with integrator
  PoD_titer <- function(titer) {
    return(
      WeightedPoD(titer,
                  vaccStat=vaccStat,
                  covInfo = covInfo,
                  paramForModelObject = paramForModelObject,
                  modelObject = modelObject,
                  fromDataTF = fromDataTF,
                  dataInput = dataInput)
    )
  }
  
  # the product of the distribution and the PoD curve
  #   with only the integration variable as the argument
  #   for use in integrator
  if (trialPOP$userDefinedTiterDist) {
    titerDist <- function(titer){return(trialPOP$titerDist(titer))}
  } else (
    titerDist <- function(titer){return(dnorm(titer, mean = trialPOP$mean, sd = trialPOP$stdDev))}
  )
  
  productFunction <- function(x){
    return(titerDist(titer = x)*PoD_titer(titer = x))
  }
  res <- quadgk(productFunction, integrationLHL, integrationRHL, tol = integrationTOL)
  return(
    res
  )
}

#' @title Builds a weighted PoD curve
#'
#' @description 
#' Build the weighted PoD curve
#' output is numeric vector of PoD values that corresponds to
#' titer list
#'
#' @param titer numeric: titer value
#' @param vaccStat vaccination status; 1 for vaccinated population, 0 for control population
#' @param covInfo user supplied covariate information
#' @param paramForModelObject PoD curve parameters
#' @param modelObject model object from \code{MakeModelObj} function that includes PoD function
#' @param fromDataTF bolean: set to TRUE if data is taken from trialPOP
#' @param dataInput NULL
#'
#' @return 
#' weighted PoD curve
#'
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#' data(covariateInfo)
#' data(PoDModel)
#' 
#' ## Example 1
#' PoDParams <- list("param.pmax" = 0.05,
#'                   "param.pmax.covariateBinary" = 1,
#'                   "param.et50" = 5,
#'                   "param.et50.covariateBinary" = 2,
#'                   "param.gamma" = 7,
#'                   "param.gamma.covariateBinary" = 1)
#' 
#' WeightedPoD(titer = 5,
#'             vaccStat = 1,
#'             covInfo = covariateInfo,
#'             paramForModelObject = PoDParams,
#'             modelObject = PoDModel,
#'             fromDataTF = TRUE,
#'             dataInput = vaccinated$subjectAttributes)
#'             
#' @importFrom dplyr %>% 
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_if
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tidyr expand_grid  
#' @importFrom modelr model_matrix        
#' @export
WeightedPoD <- function(titer,
                        vaccStat,
                        covInfo,
                        paramForModelObject,
                        modelObject,
                        fromDataTF = FALSE,
                        dataInput = NULL) { 
  
  allVars <- modelObject$vars
  allParams <- modelObject$params
  
  if (fromDataTF==TRUE & is.null(covInfo)) { # should it be here fromDataTF == FALSE
    warning("Need to include covInfo. Weights will be calculated based on population proportions.")
    stop()
  }
  
  if (is.null(covInfo)){
    covInfo <- list("var.Dummy" = list("name" = "var.Dummy",
                                       "values" = c("Universal Trait"),
                                       "prob" = c(1.0)))
  } 
  
  # make a df with rows the combination of covariates
  covStringGrid <- makeGridString(modelObject = modelObject, 
                                  covInfo = covInfo)
  
  # Make prob vec but don't cbind until after one hot encoding
  # for joint prob from data 
  if (fromDataTF) {
    totalProbabilityList <- findJointProbability(covGridDF = covStringGrid,
                                                 dataDF = dataInput)
  } else {
    totalProbabilityList <- findProb(covStringGrid,covInfo)
  }
  
  # Now include one hot encoding for each variable-value
  # Can only use modelr if there are two or more "levels" of categorical variable
  #   The lvlsUtil==1 code formats to one hot encoding "by hand"
  for (i in colnames(covStringGrid)) {
    lvlsUtil <- covStringGrid[toString(i)] %>% unique() %>% nrow()
    if (lvlsUtil==1) {
      # If only one factor then every row should be a 1
      onesList <- seq(from=1,to=1,length.out=nrow(covStringGrid))
      covStringGrid <- cbind(covStringGrid,onesList)
      # Rename in the modelr standard
      colnames(covStringGrid)[[ncol(covStringGrid)]] <- paste0(i,toString(covStringGrid[toString(i)][[1,1]]))
    } else {
      #library(modelr)
      oneHotUtilString <- paste0("oneHotUtilMat <- model_matrix(~0+",toString(i),",data=covStringGrid)")
      eval(parse(text=oneHotUtilString))
      covStringGrid <- cbind(covStringGrid,oneHotUtilMat)
    }
  }
  
  # Now include the probability of each outcome in the df
  covStringGrid <- cbind(covStringGrid,totalProbabilityList)
  colnames(covStringGrid)[ncol(covStringGrid)] <- "weight"
  
  # Include the vacc status
  covStringGrid <- covStringGrid %>% mutate("var.vaccStatus" = vaccStat)
  
  # Add parameter columns to df
  if (length(paramForModelObject)>0) {
    for (i in seq(from=1,to=length(paramForModelObject),by=1)) {
      paramColUtil <- paste0("covStringGrid <- cbind(covStringGrid,",
                             names(paramForModelObject[i]),
                             "=paramForModelObject[[i]])")
      eval(parse(text=paramColUtil))
    }
  }
  
  # Include titer -- switch to name dfOut
  dfOut <- expand_grid(covStringGrid,titer)
  colnames(dfOut)[ncol(dfOut)] <- "var.titer"
  
  ##working only for binary Covariates
  for(j in names(covInfo)){
    dfOut[j] <- as.numeric(
      dfOut[[j]] == levels(dfOut[[j]])[1])
  }
  
  # select for correct variables and parameters -> calc PoD
  PoDValues <- dfOut %>%
    dplyr::select(c(allVars,allParams)) %>%
    mutate_if(is.factor, as.numeric)
  
  PoDValues <- PoDValues %>% 
    mutate("PoD" = do.call(modelObject$predictionFun,.)) %>%  
    dplyr::select("PoD")
  
  # use weights to adjust PoD, making weighted-PoD
  dfOut <- dfOut %>% 
    mutate("PoD" = PoDValues[["PoD"]]*weight)
  
  # sum the PoDs for each titer-val
  dfOut <- dfOut %>% 
    group_by(var.titer) %>% 
    summarise("PoDTotal" = sum(PoD)) %>% 
    dplyr::select(PoDTotal)
  
  # note that these correspond to the user-input argument "titer"
  return(as.numeric(unlist(dfOut[,1])))
}

#' @title Grid of strings 
#'
#' @description Makes the grid of strings-- all covariate combinations given a model object
#   and the random covariate info
#'
#' @param modelObject model object from \code{MakeModelObj} function that includes PoD function
#' @param covInfo user supplied covariate information
#'
#' @return 
#' grid of strings
#'
#' @examples
#' ## Data preparation
#' data(covariateInfo)
#' data(PoDModel)
#' 
#' ## Example 1
#' makeGridString(modelObject = PoDModel,
#'                covInfo = covariateInfo)
#'
#' @export
makeGridString <- function(modelObject,
                           covInfo) {
  allVars <- modelObject$vars
  allParams <- modelObject$params
  
  for (j in seq(length(covInfo))) {
    if (!exists("expandGridString")) {
      expandGridString <- paste0("expand.grid(")
      expandGridString <- paste0(expandGridString,"covInfo[[",j,"]][['values']],")
    } else {
      expandGridString <- paste0(expandGridString,"covInfo[[",j,"]][['values']],")
    }
  }
  expandGridString <- stringr::str_sub(expandGridString, end=-2)
  expandGridString <- paste0("covStringGrid <- ",expandGridString,")")
  eval(parse(text=expandGridString))
  
  for (j in seq(length(covInfo))) {
    colnames(covStringGrid)[[j]] <- covInfo[[j]][['name']]
  }
  return(covStringGrid)
}


#' @title Find probability
#'
#' @description  Function to return the probability of each combination of variables
#'
#' @param df an expanded grid
#' @param covInfo a list of lists giving variable names, values, and corresponding probabilites
#'
#' @return 
#' list of probabilities
#'
#' @details 
#' Function treats each covariate as independent
#   P(a,b)=P(a)*P(b). Dependent case is handled in 
#   function findJointProbability() found below
#'
#' @examples
#' ## Data preparation
#' data(covariateInfo)
#' 
#' ## Example 1
#' covGridString <- data.frame("var.covariateBinary" = factor(c('valueA', 'valueB')))
#' 
#' findProb(covGridString,
#'          covariateInfo)
#'
#' @export
findProb <- function(df,
                     covInfo) {
  
  totalProbabilityList <- seq(from=0,to=0,length.out=nrow(df))
  for (i in seq(nrow(df))) {
    totalProbability <- 1
    for (j in seq(ncol(df))) {
      varName <- colnames(df)[[j]]
      varVal <- toString(df[[i,j]])
      factorPositions <- covInfo[[varName]][['values']]==varVal
      probOfVal <- covInfo[[varName]][['prob']][factorPositions]
      totalProbability <- totalProbability*probOfVal
    }
    totalProbabilityList[i] <- totalProbability
  }
  return(totalProbabilityList)
}

#' @title Find joint probability 
#'
#' @description Function returns joint probability of the different variable combinations based
#   on the population proportion as sees in DataDF. The combinations are given in covGridDF.
#'
#' @param covGridDF grid of strings: each row represents a possible combination of categorical covariates
#' @param dataDF data frame
#'
#' @return 
#' list of joint probabilities
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#' 
#' ## Example 1
#' covGridString <- data.frame("var.covariateBinary" = factor(c('valueA', 'valueB')))
#' 
#' findJointProbability(covGridString,
#'                      vaccinated$subjectAttributes)
#'
#' @importFrom dplyr %>%
#' @export   
findJointProbability <- function(covGridDF,
                                 dataDF) {
  
  # covGridDF is a grid of strings. Each row represents a possible
  #   combination of categorical covariates, eg male, old,... male, child,... etc 
  
  # The following loop filters the dataset 
  #   based on the rows of covGridDF-- each possible
  #   combination of covariates, then counts how many examples
  #   of that combination and divides by total number of 
  #   subjects -- output vector should sum to one
  totalProbabilityList <- seq(from=0,to=0,length.out=nrow(covGridDF)) # initialize vector
  nTotal <- dataDF %>% nrow() # count total subjects
  for (i in seq(nrow(covGridDF))) {
    totalProbabilityUtil <- dataDF
    for (colStringUtil in colnames(covGridDF)) {
      # make string to help filter for specific named covariates
      totalProbabilityString <- paste0("totalProbabilityUtil <- totalProbabilityUtil %>% filter(",
                                       colStringUtil,
                                       "== covGridDF[['",
                                       colStringUtil,
                                       "']][[i]])")
      eval(parse(text=totalProbabilityString))
    }
    # count all examples of such covariate combinations
    totalProbabilityUtil <- totalProbabilityUtil %>% nrow()
    # return population proportion which acts as joint probability
    totalProbabilityList[i] <- totalProbabilityUtil / nTotal
  }
  # sum(totalProbabilityList)
  return(totalProbabilityList)
}

#' @title Population mean jittering
#'
#' @description
#' Function jitters the mean of the population.
#'
#' Jittering is adding noise to the mean. The jittered mean is sampled from the distribution with the population mean and population standard deviation standardized by the number of subjects in the population.
#'
#' \deqn{Mean_{jitter} \sim N(mean, \frac{sd}{N} )}{ MeanJitter ~ N(mean, sd/N)}
#'
#' @param blindPopulation population - population with N, mean, stdDev statistics
#'
#' @return
#' Jittered mean - numeric value
#'
#' @usage
#' JitterMean(blindPopulation)
#'
#' @examples
#' ## Data preparation
#' data(vaccinated)
#'
#' ## Example 1
#' vaccinated$mean
#' JitterMean(vaccinated)
#'
#' @export
JitterMean <- function(blindPopulation) {
  if (!is(blindPopulation, "Population")) {
    incorrectPopulationInput("blindPopulation")
  }
  blindSample <- rnorm(1,blindPopulation$mean,blindPopulation$stdDev/sqrt(round(blindPopulation$N)))
  return(blindSample)
}

#' @title PoDBAY efficacy summary (mean, median, confidence intervals)
#'
#' @description
#' Function summarizes PoDBAY efficacy statistics (mean, median, confidence intervals) based on the set of estimated efficacies and chosen condfidence interval.
#'
#' @param efficacySet numeric vector - vector of estimated PoDBAY efficacies from \code{PoDBAYEfficacy()} function.
#' @param ci numeric - required confidence level
#'
#' @return
#' named list - mean, median, CILow, CIHigh
#'
#' @usage
#' EfficacyCI(efficacySet, ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data(efficacySet)
#'
#' ## Example 1
#' EfficacyCI(efficacySet, ci = 0.95)
#'
#' @details
#' Confidence intervals are calculated using quantiles of estimated efficacies.
#'
#' @export
EfficacyCI <- function(efficacySet, ci = 0.95) {
  CILow <- quantile(efficacySet, (1 - ci) / 2, names = F)
  CIHigh <- quantile(efficacySet, ci + (1 - ci) / 2, names = F)
  return(
    list(
      mean = mean(efficacySet),
      median = median(efficacySet),
      CILow = CILow,
      CIHigh = CIHigh
    )
  )
}

#' @title PoD curve parameters estimation; works with covariate data.
#'
#' @description
#' Function estimates the PoD curve parameters (as specified in \code{PoDModel}) using \code{PoDMLE} function. Number of PoD curves estimated equals to the repeatCount input parameter.
#'
#' The estimation is performed using provided subject-level data in diseased and non-diseased \code{Population-class} objects.
#'
#' @param diseased \code{Population-class} object: diseased population, contains subject level titers
#' @param nondiseased \code{Population-class} object: non-diseased population from immunogenicity subset, contains subject level titers
#' @param nondiseasedGenerationCount numeric: total number of non-diseased subjects in the clinical trial
#' @param PoDModel list: information on the PoD curve function as well as user-defined boundaries for the optimization
#' @param repeatCount numeric: how many times is the dataset bootstrapped and the PoD curve parameter estimation performed
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#' @param MLEinitialGuess NULL
#'
#' @return
#' results: PoD curve parameters after resetting the disease status, named data.frame of estimated PoD curve parameters (as specified in \code{PoDModel}); see details for more information
#'
#' resultsPriorReset: PoD curve parameters prior to resetting the status, named data.frame of estimated PoD curve parameters (as specified in \code{PoDModel}); see details for more information
#'
#' failcount: number of iterations in which MLE failed to converge; see details for more information
#'
#' @usage
#' PoDParamEstimation(diseased,
#'                    nondiseased,
#'                    nondiseasedGenerationCount,
#'                    PoDModel,
#'                    repeatCount = 500,
#'                    adjustTiters = FALSE,
#'                    adjustFrom = log2(10),
#'                    adjustTo = log2(5),
#'                    MLEinitialGuess = NULL)
#'
#' @examples
#' ## Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(PoDModel)
#'
#' ## Example 1
#' # Creating imunogenicity subset, method = "Full"
#' ImmunogenicitySample <- 
#'     BlindSampling(diseased, 
#'                          nondiseased, 
#'                          method = list(name = "Full", 
#'                                        value = "NA"))
#'
#' # Number of all non-diseased subjects in the clinical trial
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' PoDParamEstimation(diseased,
#'                    ImmunogenicitySample$immunogenicityNondiseased,
#'                    nondiseasedGenerationCount,
#'                    PoDModel,
#'                    repeatCount = 3)
#'
#' ## Example 2
#' # Creating imunogenicity subset, method = "Fixed", value = 150
#' ImmunogenicitySample <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Fixed", 
#'                                        value = 150))
#'                                        
#' # Number of all non-diseased subjects in the clinical trial
#' nondiseasedGenerationCount <- nondiseased$N
#'
#' PoDParamEstimation(diseased,
#'                    ImmunogenicitySample$immunogenicityNondiseased,
#'                    nondiseasedGenerationCount,
#'                    PoDModel,
#'                    repeatCount = 3)
#'
#' @details
#'
#' diseased: \code{Population-class} object of diseased population, contains subject level titers
#'
#' nondiseasedTiters:  \code{Population-class} object of non-diseased population in the immunogenicity subset, contains subject level titers
#'
#' There are two possible scenarios
#'
#' \itemize{
#'   \item Full: Full information about non-diseased titers is available, i.e subject level data for all non-diseased subjects from the clinical trial (nondiseasedGenerationCount = number of all non-diseased subjects in the clinical trial).
#'
#'   \item Ratio or Fixed: Information about non-diseased titers is available only for the immunogenicity subset. In order to compensate for these missing titers we upsampling of this subset to the total number of non-diseased (nondiseasedGenerationCount) in the trial is needed. 
#'
#' }
#'
#' nondiseasedGenerationCount: number of all non-diseased subjects in the clinical trial
#'
#' NOTE: Number of estimated parameters can be lower than repeatCount as MLE does not necessary converge in all estimations; failcount (number of iterations in which MLE failed to converge) is also returned; for details see \code{MLE} function.
#'
#'
#' Function steps
#'
#' \itemize{
#'   \item Upsample non-diseased if needed (needed for method Fixed) - from immunogenicity subset size (N = NondiseasedImmunogenicitySubset$N) to the whole trial size (N = nondiseasedGenerationCount). For details see \code{GenerateNondiseased} function.
#'
#'   \item Estimate PoD curve: resultsPriorReset
#'
#'   \item Reset disease status: the purpose is to estimate the confidence intervals of the PoD curve and its parameters
#'
#'   Part of the reset disease status procedure is the non-parametric bootstrap: titers of diseased and non-diseased subjects are pooled, and associated PoDs are calculated using their titer values and estimated PoD curve. Based on the subject level probabilities (PoDs), the disease status is re-estimated.
#'
#'   \item Re-estimate PoD curve: new diseased and non-diseased titers are used to reestimate the PoD curve
#' }
#'
#' @export
PoDParamEstimation <- function(diseased,
                               nondiseased,
                               nondiseasedGenerationCount,
                               PoDModel,
                               repeatCount = 500,
                               adjustTiters = FALSE,
                               adjustFrom = log2(10),
                               adjustTo = log2(5),
                               MLEinitialGuess = NULL) {
  
  if (floor(nondiseasedGenerationCount) < length(nondiseased$titers)) {
    warning(paste("The input value for \"nondiseasedGenerationCount\" is lower than number of nondiseasedTiters"))
  }
  
  if (adjustTiters) {
    diseased$titers[diseased$titers<=adjustFrom] <- adjustTo
    nondiseased$titers[nondiseased$titers<=adjustFrom] <- adjustTo
  }
  
  
  results <- data.frame()
  resultsPriorReset <- data.frame()
  
  failCount <- 0
  estimateOnThresholdWarningPriorReset <- 0
  estimateOnThresholdWarning <- 0
  for (i in 1:repeatCount) {
    print(paste("iteration", i))
    
    # Upsample non-diseased immunogenicity subset to the size of the trial
    upsampledNondiseased <- GenerateNondiseased(nondiseased,
                                                nondiseasedGenerationCount) 
    
    # initialParams argument for PoDMLE allows for specified starting parameters
    # Estimate PoD curve parameters
    estimatedPoDModel <- PoDMLE( 
      nondiseased = upsampledNondiseased,
      diseased = diseased,
      PoDModel = PoDModel,
      adjustTiters = adjustTiters,
      adjustFrom = adjustFrom,
      adjustTo = adjustTo,
      initialParams = MLEinitialGuess
    )
    
    
    # if the estimation fails, note the failure and skip to the next iteration
    if (is.null(estimatedPoDModel)) {
      failCount <- failCount + 1
      next
    }
    
    if (checkBounds(values = estimatedPoDModel$valueParams,
                    ranges = PoDModel$rangeParams) > 0) {
      estimateOnThresholdWarningPriorReset <- estimateOnThresholdWarningPriorReset + 1
      warning("Estimation of one or more parameters converged to the limit. Consider changing the permitted ranges of your parameters. ")
    }
    
    ## Reset disease status ##
    
    # Pool non-diseased and diseased data
    completePopulation <- mergePopulations(diseased, upsampledNondiseased)
    
    # Bootstrap all data
    bootstrappedPopulation <- samplePopulation(completePopulation, 
                                               length(completePopulation$titers), 
                                               replace = T)
    
    # Assign to each subject a probability of disease
    bootstrappedPopulation$titers[bootstrappedPopulation$titers<0]<-0 # probably no longer necessary because of adjustTiter if-statement above
    bootstrappedPopulation$assignPoD(
      PoD(titers = bootstrappedPopulation$titers,
          attributes = bootstrappedPopulation$subjectAttributes,
          PoDModel = estimatedPoDModel,
          valueParams = estimatedPoDModel$valueParams,
          adjustTiters = adjustTiters,
          adjustFrom = adjustFrom,
          adjustTo = adjustTo)
    )
    
    # Re-assign the disease status based on the probabilities
    dStatus <- mapply(rbinom, n = 1, size = 1, prob = bootstrappedPopulation$PoDs)
    bootstrappedPopulation$diseaseStatus <- sapply(dStatus, numToBool)
    
    # New diseased and non-diseased populations
    newDiseased <- ExtractDiseased(bootstrappedPopulation) 
    newNondiseased <- ExtractNondiseased(bootstrappedPopulation) 
    
    # Create non-diseased immunogenicity sample of the new population
    # choose appropriate method for the immunogenicty sample creation = Full or Fixed based on the input to the function
    # if input FULL SAMPLE -> method = "FULL"
    # otherwise method = "Fixed", method value = determined based on the data inputs #nondiseased
    
    if (round(nondiseasedGenerationCount) == length(nondiseased$titers)) {
      method <- list(name = "Full",
                     value = NA)
    } else {
      method <- list(name = "Fixed",
                     value = length(nondiseased$titers))
    }
    
    newImmunogenicitySample <- ImmunogenicitySubset(newDiseased,
                                                    newNondiseased,
                                                    method = method)
    
    # Get nondiseased population from immunogenicitySample population
    newImmunogenicitySampleNondiseased <- ExtractNondiseased(newImmunogenicitySample)
    
    # Upsample the population to the total number of non-diseased subjects in newNoniseased
    newUpsampledNondiseased <- GenerateNondiseased(newImmunogenicitySampleNondiseased,
                                                   newNondiseased$N) 
    
    #initialParams argument allows to specify initial guess for optimization
    # estimate new PoD curve parameters
    estimatedPoDModel2 <- PoDMLE(
      nondiseased = newUpsampledNondiseased,
      diseased = newDiseased,
      PoDModel = PoDModel,
      adjustTiters = adjustTiters,
      adjustFrom = adjustFrom,
      adjustTo = adjustTo,
      initialParams = MLEinitialGuess
    )
    
    if (is.null(estimatedPoDModel2)) {
      failCount <- failCount + 1
      next
    }
    
    if (checkBounds(values = estimatedPoDModel2$valueParams,
                    ranges = PoDModel$rangeParams) > 0) {
      estimateOnThresholdWarning <- estimateOnThresholdWarning + 1
      warning("Estimation of one or more parameters converged to the limit. Consider changing the permitted ranges of your parameters. ")
    }
    
    resultsPriorReset <- rbind(resultsPriorReset,
                               unname(unlist(estimatedPoDModel$valueParams)))
    
    
    results <- rbind(results,
                     unname(unlist(estimatedPoDModel2$valueParams)))
    
  }
  
  if (nrow(resultsPriorReset)>0) {
    names(resultsPriorReset) <- PoDModel$params
  }
  
  if (nrow(results)>0) {
    names(results) <- PoDModel$params
  }
  
  return(list(results = results,
              resultsPriorReset = resultsPriorReset,
              failCount = failCount,
              estimateOnThresholdWarningPriorReset = estimateOnThresholdWarningPriorReset,
              estimateOnThresholdWarning = estimateOnThresholdWarning)
  )
}

#' @title Check if optimization landed on a predefined boundary.
#'
#' @description  Function returns number of optimization boundaries hit for given PoD curve parameters (as specified in \code{PoDModel}). 
#' For a 6-parameter model, the maximal output of this function is 6.
#'
#' @param values list: named list of PoD curve parameters, output of optimization
#' @param ranges list: named list of ranges (lower bound and upper bound) for each parameter,used in the optimization
#'
#' @return
#' count: numeric, value of number the boundary was hit 
#'
#' @usage
#' checkBounds(values, ranges)
#'
#' @examples
#' ## Data preparation
#' data(nondiseased)
#' data(diseased)
#' data(PoDModel)
#'
#' ## Example 1
#' # Optimization
#' set.seed(2)
#' PoDParamEst <- PoDParamEstimation(diseased,
#'                                   nondiseased,
#'                                   nondiseased$N,
#'                                   PoDModel,
#'                                   repeatCount = 1,
#'                                   adjustTiters = FALSE)
#'
#' # Check number boundariew were hit
#' checkBounds(values = PoDParamEst$results,
#'             ranges = PoDModel$rangeParams)
#' @export
checkBounds <- function(values,ranges){
  counter = 0
  
  for (paramName in names(values)) {
    atBoundMin <- values[[paramName]] == ranges[[paramName]][[1]]
    atBoundMax <- values[[paramName]] == ranges[[paramName]][[2]]
    
    if(sum(atBoundMin)+sum(atBoundMax)>0){
      counter = counter + sum(atBoundMin) + sum(atBoundMax)
    }
  }
  return(counter)
}

#' @title Generation of upsampled non-diseased population; works with covariate data.
#'
#' @description Function upsamples (by random sampling with replacement) the immunogenicity subset population to the required size.
#'
#' If the size of the immunogenicity subset matches the required size, nothing happens and the original immunogenicity subset population is returned.
#'
#' @param immunogenicityNondiseased \code{Population-class} object: immunogenicity non-diseased population, contains subject level titers
#' @param nondiseasedGenerationCount numeric: total number of non-diseased subjects, required size of the non-diseased population
#'
#' @return 
#' generatedNondiseased: \code{Population-class} object of the upsampled non-diseased population, contains subject-level titers 
#'
#' @usage
#' GenerateNondiseased(immunogenicityNondiseased, nondiseasedGenerationCount)
#'
#' @examples
#' ## Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(PoDModel)
#'
#' ## Example 1
#' # Creating imunogenicity subset, method = "Fixed"
#' ImmunogenicitySample <- 
#'     BlindSampling(diseased, 
#'                   nondiseased, 
#'                   method = list(name = "Fixed", 
#'                                 value = 100))
#'
#' # Number of all non-diseased subjects in the clinical trial
#' nondiseasedGenerationCount <- nondiseased$N
#' 
#' # Upsampling of non-diseased titers
#' GenerateNondiseased(ImmunogenicitySample$immunogenicityNondiseased, nondiseasedGenerationCount)
#'
#' @details
#' The inputs should come from immunogenicity subset. \code{nondiseasedGenerationCount} represents number of all non-diseased patients in the clinical trial.
#'
#' Immunogenicity subset population is obtained from function \code{BlindSampling}. Immunogenicity subset represents a sample from the non-diseased population.
#'
#' In this function, sampling with replacement to the required \code{nondiseasedGenerationCount} of the immunogenecity subset is performed. The function is used inside \code{PoDParamEstimation} function.
#'
#' @export
GenerateNondiseased <- function(immunogenicityNondiseased, nondiseasedGenerationCount) {
  
  generatedNondiseased <-  if (round(length(immunogenicityNondiseased$titers)) == floor(nondiseasedGenerationCount)) {
    immunogenicityNondiseased
  } else {
    samplePopulation(immunogenicityNondiseased, nondiseasedGenerationCount, replace = TRUE)
  }
  
  return(generatedNondiseased = generatedNondiseased)
}

#' @title Setup for the maximum likelihood estimation (MLE); works with covariate data.
#'
#' @description Function estimates the optimal PoD curve parameters (as specified in \code{PoDModel}) using diseased and non-diseased populations and their subject-level titers. 
#' Initial guess of the PoD curve parameters can be provided as an input to the optimization. Otherwise, it is dynamically calculated using the provided datasets.
#'
#' @param diseased \code{Population-class} object: diseased population, contains subject level titers
#' @param nondiseased \code{Population-class} object: non-diseased population from immunogenicity subset, contains subject level titers
#' @param PoDModel list: information on the PoD curve function as well as user-defined boundaries for the optimization
#' @param initialParams TBD
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#' @param paramH0 TBD
#'
#' @return
#' list("predictionFun","vars","params","valueParams","AIC"), named list of PoD parameters and AIC, if MLE converges.
#'
#' Null, if MLE does not converge.
#'
#' @usage PoDMLE(nondiseased,
#'               diseased,
#'               PoDModel,
#'               initialParams = NULL,
#'               adjustTiters = FALSE,
#'               adjustFrom = log2(10),
#'               adjustTo = log2(5),
#'               paramH0 = 0)
#'
#' @details
#'
#' PoDMLE function estimates the PoD curve parameters by maximizing the likelihood value (see \code{MLE} function for details) based on the provided titers for diseased and non-diseased subjects.
#'
#' The \code{optim} function is used for optimization with method = "L-BFGS-B", 500 maximum iterations, and user-defined (or dynamically balculated based on provided datasets) boundaries for optimization.
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#' data(PoDModel)
#'
#' # Example 1
#' PoDMLE(nondiseased, diseased, PoDModel)
#'
#' @export
PoDMLE <- function(nondiseased,
                   diseased,
                   PoDModel,
                   initialParams = NULL,
                   adjustTiters = FALSE,
                   adjustFrom = log2(10),
                   adjustTo = log2(5),
                   paramH0 = 0) {
  
  # Get the mins and maxs for optim from user input
  boundsIndex <- 1
  myLower <- c()
  myUpper <- c()
  for (paramName in PoDModel$params) {
    myLower[boundsIndex] <- PoDModel$rangeParams[[paramName]][[1]]
    myUpper[boundsIndex] <- PoDModel$rangeParams[[paramName]][[2]]
    boundsIndex <- boundsIndex + 1
  }
  
  ## Initial guess of pmax, et50, gamma: use the base (3-parameter) PoD model estimates
  if (is.null(initialParams)) {
    baseInitialParams <- BaseModelFit(nondiseased,
                                      diseased,
                                      adjustTiters = adjustTiters,
                                      adjustFrom = adjustFrom,
                                      adjustTo = adjustTo)
    if (is.null(baseInitialParams)) {
      # Use ranges? if optimization in BaseModelFit doesn't converge
      baseInitialParams <- list("param.et50"  = PoDModel$rangeParams$param.et50[1] ,
                                "param.gamma" = PoDModel$rangeParams$param.gamma[1],
                                "param.pmax"  = PoDModel$rangeParams$param.pmax[1])
      
    }
    covariateParams <- PoDModel$params[!grepl(paste0("param.pmax$|param.gamma$|param.et50$"), PoDModel$params)]
    initialParams <- baseInitialParams
    if (length(covariateParams)>0){
      for (i in covariateParams) {
        initialParams[[i]] <- 1
      }
    } else {
      initialParams <- baseInitialParams
    }
  } 
  
  # Sort the initialParams list by the PoDModel$params vector
  initialParams <- initialParams[order(match(names(initialParams),PoDModel$params))]
  
  fitOptionsList <- list(par = initialParams,
                         fn = NLLPoD_JD,
                         PoDModel = PoDModel,
                         nondiseased = nondiseased, 
                         diseased = diseased,
                         lower = myLower,
                         upper = myUpper,
                         method = "L-BFGS-B",
                         control = list(maxit=5000),
                         hessian = TRUE) 
  
  TryCatchError <- FALSE
  tryCatch({
    # Perform actual fit
    fit <- do.call(optim,fitOptionsList)
    # Calculate AIC
    aicFit <- 2*fit$value + 2*length(PoDModel$params)
    # Return the p-values using our wald test function
    # paramInfoDF <- WaldStat(fit$par,
    #                         fit$hessian,
    #                         paramH0 = paramH0)
    # Return parameter values
    params <- fit$par
  },
  error = function(x){
    TryCatchError <<- TRUE
  })
  
  if(TryCatchError){
    return (NULL)
  }
  
  paramsList <- as.list(params)
  
  return(list("predictionFun" = PoDModel$predictionFun,
              "vars" = PoDModel$vars,
              "params" = PoDModel$params,
              "valueParams" = paramsList,
              "AIC" = aicFit))
}

#' @title Base model fit.
#'
#' @description Function calculates initial parameter guesses for PoDMLE optimization using base PoD model fit.
#'
#' @param diseased \code{Population-class} object: diseased population, contains subject level titers
#' @param nondiseased \code{Population-class} object: non-diseased population from immunogenicity subset, contains subject level titers
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#' 
#' 
#' @return
#' list("param.et50", "param.gamma", "param.pmax") - named list of PoD parameters, if MLE converges.
#'
#' Null, if MLE does not converge.
#'
#' @usage BaseModelFit(nondiseased,
#'                     diseased,
#'                     adjustTiters = FALSE,
#'                     adjustFrom = log2(10),
#'                     adjustTo = log2(5))
#'
#' @examples
#' ## Data preparation
#' data(diseased)
#' data(nondiseased)
#' 
#' ## Example 1
#' BaseModelFit(nondiseased,
#'              diseased,
#'              adjustTiters = TRUE,
#'              adjustFrom = 0,
#'              adjustTo = 0)
#'
#' @importFrom stats runif
#' @export
BaseModelFit <- function(nondiseased,
                         diseased,
                         adjustTiters = FALSE,
                         adjustFrom = log2(10),
                         adjustTo = log2(5)) {
  
  # pmax for base model fitting
  titers <- c(nondiseased$titers, diseased$titers)
  diseaseStatus <- c(rep(0, length(nondiseased$titers)),
                     rep(1, length(diseased$titers)))
  newOrder <- order(titers)
  titers <- titers[newOrder]
  diseaseStatus <- diseaseStatus[newOrder]
  lowTiterPercent <- 0.2
  lowTitersStatus <- diseaseStatus[1:round(length(diseaseStatus) * lowTiterPercent)]
  initialPmax <- (sum(lowTitersStatus) + 0.5) / (length(lowTitersStatus) + 0.5)
  
  # gamma for base model fitting
  coeff <- 9.1902
  range <- max(titers) - min(titers)
  upperBoundGamma <- coeff * 50 / range
  # option A
  initialGamma_A <- 6
  # option B
  initialGamma_B <- initialGamma_A + runif(1, min=-3, max=15)
  
  # et50 for base model fitting
  # option A
  tryError_A <- F
  tryCatch(
    {
      RootFunction <- function(x) dnorm(x, 
                                        mean = mean(diseased$titers), 
                                        sd = sd(diseased$titers)) - dnorm(x, 
                                                                          mean = mean(nondiseased$titers), 
                                                                          sd = sd(nondiseased$titers))
      RootSearch <- uniroot(RootFunction, interval = c(mean(diseased$titers), mean(nondiseased$titers)))
      initialEt50_A <- RootSearch$root
    }, 
    error = function(e) {
      tryError_A <<- T
    }
  )
  # option B
  initialEt50_B <- median(titers[titers > 0])
  
  # Find parameters, option A
  tryCatch(
    {  
      init <- c(initialEt50_A, initialGamma_A, initialPmax)
      names(init) <- c("et50", "gamma", "pmax")
      
      estimatedParameters_A <- optim(
        par = init,
        fn = MLE,
        nondiseasedTiters = nondiseased$titers,
        diseasedTiters = diseased$titers,
        adjustTiters = adjustTiters,
        adjustFrom = adjustFrom,
        adjustTo = adjustTo,
        lower = c(1e-6, -upperBoundGamma, 1e-6),
        upper = c(Inf, upperBoundGamma, 1),
        method = "L-BFGS-B",
        control = list(fnscale = -1, maxit = 500))
    },
    error = function(e) {
      tryError_A <<- T
    }
  )
  
  tryError_B <- F
  if (tryError_A) {
    # Find parameters, option B 
    initialParams <- c(initialEt50_B, initialGamma_B, initialPmax)
    names(initialParams) <- c("et50", "gamma", "pmax")
    tryCatch(
      estimatedParameters_B <- optim(
        par = initialParams,
        fn = MLE,
        nondiseasedTiters = nondiseased$titers,
        diseasedTiters = diseased$titers,
        adjustTiters = adjustTiters,
        adjustFrom = adjustFrom,
        adjustTo = adjustTo,
        lower = c(1e-6, -upperBoundGamma, 1e-6),
        upper = c(Inf, upperBoundGamma, 1),
        method = "L-BFGS-B",
        control = list(fnscale = -1, maxit = 500)),
      error = function(e) {
        tryError_B <<- T
      }
    )
  }
  if (tryError_B) {return(NULL)}
  if (tryError_A&!tryError_B) {
    return(list("param.et50"  = estimatedParameters_B$par[[1]],
                "param.gamma" = estimatedParameters_B$par[[2]],
                "param.pmax"  = estimatedParameters_B$par[[3]]))
  } else {
    return(list("param.et50"  = estimatedParameters_A$par[[1]],
                "param.gamma" = estimatedParameters_A$par[[2]],
                "param.pmax"  = estimatedParameters_A$par[[3]]))
    
  }
}

#' @title Maximum Likelihood estimation; works with covariate data.
#'
#' @description Function calculates the log likelihood value for base model which is used after the initial guesses of the parameters are set in the \code{PoDMLE} function.
#'
#' @param params named numeric vector: PoD curve parameters (et50, gamma, pmax)
#' @param nondiseasedTiters numeric vector: non-diseased subjects titers
#' @param diseasedTiters numeric vector: diseased subjects titers
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return log likelihood, numeric value
#'
#' @usage
#' MLE(params,
#'     nondiseasedTiters,
#'     diseasedTiters,
#'     adjustTiters = FALSE,
#'     adjustFrom = log2(10),
#'     adjustTo = log2(5))
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#' PoDParams <- c("pmax" = 0.03,
#'                "et50" = 5.3,
#'                "gamma" = 15)
#'
#' # MLE calculation
#' MLE(PoDParams, nondiseased$titers, diseased$titers)
#'
#' @details MLE function is used inside of PoDMLE function and esimates the PoD curve parameters.
#'
#' Based on the provided titers for diseased and non-diseased subjects the PoD curve parameters which maximize the log likelihood are chosen as optimal estimates of parameters.
#'
#' @export
MLE <- function(params,
                nondiseasedTiters,
                diseasedTiters,
                adjustTiters = FALSE,
                adjustFrom = log2(10),
                adjustTo = log2(5)) {
  
  if(any(is.na(match(names(params), c("pmax", "et50", "gamma") )))){
    stop(paste("The input value for params is incorrect. 'params' parameter has wrong names."))
  }
  
  pmax <- as.double(params["pmax"])
  et50 <- as.double(params["et50"])
  gamma <- as.double(params["gamma"])
  
  probDiseased <- PoDBaseModel(titer = diseasedTiters,
                               pmax, et50, gamma,
                               adjustTiters = adjustTiters,
                               adjustFrom = adjustFrom,
                               adjustTo = adjustTo)
  
  probNondiseased <- PoDBaseModel(titer = nondiseasedTiters,
                                  pmax, et50, gamma,
                                  adjustTiters = adjustTiters,
                                  adjustFrom = adjustFrom,
                                  adjustTo = adjustTo)
  
  logLikelihood <- sum(log(probDiseased)) + sum(log(1 - probNondiseased))
  return(logLikelihood)
}

#' @title Confidence intervals of PoD curve parameters. 
#'
#' @description Function calculates confidence intervals (user-defined) of the PoD curve parameters 
#'
#' @param estimatedParameters output of \code{PoDParamEstimation} function
#' @param ci numeric: value from (0, 1) interval, confidence level of interest
#'
#' @return CI of all PoD curve parameters
#'
#' @usage
#' PoDParamsCI(estimatedParameters, ci = 0.95)
#'
#' @export
PoDParamsCI <- function(estimatedParameters, ci = 0.95) {
  lower <- lapply(estimatedParameters, quantile, (1 - ci) / 2, names = F)
  upper <- lapply(estimatedParameters, quantile, ci + (1 - ci) / 2, names = F)
  
  return(
    list(
      lowerBounds = lower,
      upperBounds = upper)
  )
}

#' @title PoD curve point estimate; works with covariate data.
#'
#' @description Function returns PoD curve parameters corresponding to the point estimate of PoD curve.
#'
#' @param resultsPriorReset named data frame: set of estimated PoD curve parameters before resetting the disease status; for further details see \code{PoDParamEstimation} function.
#' @param PoDModel list: information on the PoD curve function
#' @param titers numeric vector: a grid of titers for PoD curve point estimate calculation
#' @param optim_titers logical: TRUE for a predefined sequence of titers 
#'
#' @return
#' paramsPointEstimate: named data frame of PoD curve parameters corresponding to the PoD curve point estimate
#'
#' @usage
#' PoDParamPointEstimation(resultsPriorReset, 
#'                         PoDModel,
#'                         titers = seq(from = 0, to = 20, by = 0.01), 
#'                         optim_titers = FALSE)
#'
#' @details
#' For each of estimated PoD curves in resultsPriorReset, the function values (probabilities of disease, PoD) for provided grid of titers are calculated.
#'
#' Median of function values (PoDs) at each provided titer is calculated.
#'
#' Subsequently, the PoD curve model is fitted to the median datapoints using \code{fitPoD} function, in order to get PoD curve parameters close to this median curve.
#' 
#'
#' @examples
#' ## Data preparation
#' data(estimatedParameters)
#' data(PoDModel)
#'
#' ## Example 1
#' # titers for which we want to optimize the functional values
#' titers <- seq(from = 0, to = 20, by = 0.01)
#'
#' # Point estimate of PoD curve
#' PoDParamPointEstimation(estimatedParameters$resultsPriorReset, PoDModel, titers)
#'
#' @export
PoDParamPointEstimation <- function(resultsPriorReset,
                                    PoDModel,
                                    titers = seq(from = 0, to = 20, by = 0.01),
                                    optim_titers = FALSE){
  
  
  
  # Generate vector of titers
  if (optim_titers) {
    resultsMed <- lapply(resultsPriorReset, median)
    titers <- c(seq(0, resultsMed$param.et50 - 1, by = 0.1),
                seq(resultsMed$param.et50 - 0.99, resultsMed$param.et50 + 0.99, by = 0.01),
                seq(resultsMed$param.et50 + 1, 2 * resultsMed$param.et50, by = 0.1))
  }
  
  # Generate dataframe of attributes
  covariateVars <- PoDModel$vars[!grepl(paste0("var.titer$"), PoDModel$vars)]
  attributes <- data.frame(matrix(ncol = length(covariateVars), nrow = length(titers)))
  
  if (length(covariateVars)>0){
    colnames(attributes) <- covariateVars
    for (i in 1:length(covariateVars)) {
      covariateValues <- rbinom(length(titers),1,0.5) # sampling from binom distribution should be used only if BINARY covariate is studied
      attributes[,i] <- covariateValues
    }
  }
  
  # Number of estimated PoD curves
  nCurves <- nrow(resultsPriorReset)
  
  # For each estimated PoD curve, calculate function values
  functionValues <- matrix(NA, nrow = nCurves, ncol = length(titers))
  for (i in 1:nCurves) { 
    # Populate PoD model valueParams with parameter estimation results
    valueParams <- as.list(resultsPriorReset[i,])
    
    # Assign probabilities based on param values
    functionValues[i, ] <- PoD(titers = titers, 
                               attributes = attributes,
                               PoDModel = PoDModel,
                               valueParams = valueParams,
                               adjustTiters = FALSE) # we are fitting the curve so there is no need to adjust titers
  }
  
  # Median of function values
  medianCurve <- apply(functionValues, 2, median)
  
  # Median PoD curve parameters
  medianParams <- apply(resultsPriorReset, 2, median)
  
  # Set initial guess of PoD curve parameters for optimization
  #init <- as.list((1 - (medianParams["param.gamma"] / 100)) * medianParams)
  init <- as.list(medianParams[PoDModel$params])
  
  # Get the mins and maxs for optimization from user input
  boundsIndex <- 1
  myLower <- c()
  myUpper <- c()
  for (paramName in PoDModel$params) {
    myLower[boundsIndex] <- PoDModel$rangeParams[[paramName]][[1]]
    myUpper[boundsIndex] <- PoDModel$rangeParams[[paramName]][[2]]
    boundsIndex <- boundsIndex + 1
  }
  
  # Find optimal PoD curve
  paramsPointEstimate <- optim(
    par = init,
    fn = fitPoD,
    titersInput = titers,
    attributes = attributes,
    PoDModel = PoDModel,
    medianCurve = medianCurve,
    lower = myLower,
    upper = myUpper,
    method = "L-BFGS-B",
    control = list(maxit = 1000, factr = 1e4)
  )
  
  return(paramsPointEstimate = as.list(paramsPointEstimate$par))
}

#' @title PoD curve, fitting function; works with covariate data.
#'
#' @description Function calculates the root mean squared error (RMSE) 
#' between provided PoD values and calculated PoD values. 
#' The latter are calculated using for provided titers and provided PoD curve parameters.
#'
#' By using the input titers \code{PoDParamPointEstimation} function and median 
#' of the estimated set of PoD curve parameters (output of \code{PoDParamEstimation} function), 
#' the point estimate of PoD curve can be obtained 
#' (for details see \code{PoDParamPointEstimation} function).
#'
#' @param params named data frame ("pmax", "slope", "et50"): provided PoD curve parameters
#' @param titersInput numeric vector: provided titers 
#' @param attributes named data frame: subjectAttributes for a given population
#' @param PoDModel list: information on the PoD curve function
#' @param medianCurve numeric vector: provided PoD values 
#'
#' @return
#' negative RMSE
#'
#' @usage
#' fitPoD(params, titersInput, attributes, PoDModel, medianCurve)
#'
#' @details
#'
#' \deqn{RMSE = \sqrt{\frac{\sum_{i}^{N} (PoD_{median}(titers) - 
#' PoD_{optimized}(titers))^2}{N}}}{ RMSE = sqrt( mean( (PoDmedian(titers) - PoDoptimized(titers))2 ) )}
#'
#' @examples
#' ## Data preparation
#' data(PoDModel)
#' data(vaccinated)
#' data(estimatedParameters)
#'
#' ## Example 1
#'
#' # grid of titers
#' titersInput <- seq(from = 0.01, to = 10, by = 0.01)
#'
#'
#' # functional values corresponding to the median of the estimated PoD curve parameters
#' medianCurve <- PoD(titers = titersInput,
#'                    attributes = vaccinated$subjectAttributes,
#'                    PoDModel = PoDModel,
#'                    valueParams = PoDModel$valueParams)
#'
#' # squared error of CurveTitersMedian and functional values of "params" curve
#' fitPoD(estimatedParameters$results[1,], 
#'        titersInput, 
#'        vaccinated$subjectAttributes, 
#'        PoDModel, 
#'        medianCurve)
#'
#' @export
fitPoD <- function(params, titersInput, attributes, PoDModel, medianCurve) {
  
  
  # functional values
  fittedCurve <- PoD(titers = titersInput, 
                     attributes = attributes,
                     PoDModel = PoDModel,
                     valueParams = params)
  
  # RMSE
  error <- sqrt(mean( (medianCurve - fittedCurve) ^ 2) )
  
  return(error)
}

#' @title Negative log likelihood for a general PoD curve model.
#'
#' @description Function returns the value of negative log likelihood of data given a provided PoD curve model.
#'
#' @param params named data frame ("pmax", "slope", "et50"): provided PoD curve parameters
#' @param PoDModel list: information on the PoD curve function
#' @param nondiseased \code{Population-class} object: non-diseased population from immunogenicity subset, contains subject level titers
#' @param diseased \code{Population-class} object: diseased population, contains subject level titers
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return
#' negative log likelihood
#'
#' @usage 
#' NLLPoD_JD(params,
#'           PoDModel,
#'           nondiseased,
#'           diseased,
#'           adjustTiters = FALSE, 
#'           adjustFrom = log2(10), 
#'           adjustTo = log2(5))
#'
#' @examples
#' ## Data preparation
#' data(PoDModel)
#' data(nondiseased)
#' data(diseased)
#'
#' ## Example 1
#'
#' params <- PoDModel$valueParams
#' NLLPoD_JD(params,
#'           PoDModel,
#'           nondiseased,
#'           diseased,
#'           adjustTiters = FALSE)
#' @export
NLLPoD_JD <- function (params,
                       PoDModel,
                       nondiseased,
                       diseased,
                       adjustTiters = FALSE, 
                       adjustFrom = log2(10), 
                       adjustTo = log2(5)) {
  
  PoDModel$valueParams <- params
  probDiseased <- PoD(titers = diseased$titers, 
                      attributes = diseased$subjectAttributes, 
                      PoDModel = PoDModel, 
                      valueParams = PoDModel$valueParams,
                      adjustTiters, 
                      adjustFrom, 
                      adjustTo)
  probNondiseased <- PoD(titers = nondiseased$titers, 
                         attributes = nondiseased$subjectAttributes, 
                         PoDModel = PoDModel, 
                         valueParams = PoDModel$valueParams,
                         adjustTiters, 
                         adjustFrom, 
                         adjustTo )
  logLikelihood <- sum(log(probDiseased)) + sum(log(1 - probNondiseased), na.rm=TRUE)
  NLL <- -logLikelihood
  return(NLL)
}

#' @title PoD curve plot; works with covariate data.
#'
#' @description Supplementary function for plotting the PoD curve with the confidence ribbon (of a required level). Input values are related to PoDBAY package structure.
#' 
#' @param titers numeric vector: grid of titers at which the confidence ribbon should be calculated
#' @param covInfo list: name: character vector with the name of covariate, values: vector of factors with values of covariate, prob: vector of probability of each value
#' @param estimatedParameters estimatedParameters list: 2 named data frames needed ("results", "resultsPriorReset") needed; set of estimated PoD curve parameters, output of \code{PoDParamEstimation} function
#' @param PoDModel list: predictionFun: function call, vars: character vector with names of variables (var.titer and currently just 1 additional binary covariate), params: character vector with names of parameters 
#' @param ci numeric, required confidence level
#'
#' @return
#' PoD curve plot
#'
#' @usage 
#' PoDCurvePlot(titers,
#'              covInfo,
#'              estimatedParameters,
#'              PoDModel,
#'              ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data(estimatedParameters)
#' data(PoDModel)
#' data(covariateInfo)
#' 
#' ## Example 1
#' PoDCurvePlot(titers = seq(0, 20, by = 0.01),
#'              covInfo = covariateInfo,
#'              estimatedParameters,
#'              PoDModel,
#'              ci = 0.95)
#'
#' @export
PoDCurvePlot <- function(titers = seq(0, 20, by = 0.01),
                         covInfo,
                         estimatedParameters,
                         PoDModel,
                         ci = 0.95){
  
  # plot format setting
  size <- 12
  themePlot <-
    theme(
      plot.title = element_text(size=size,hjust = 0.5),
      axis.text.y = element_text(size=size),
      axis.text.x = element_text(size=size),
      axis.title.y = element_text(size=size),
      axis.title.x = element_text(size=size),
      legend.text = element_text(size = size),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
    )
  
  #calculation of median curve and the ribbon confidence interval in given points (titer values)
  PoDCIdfs <- PoDCIInPoints(titers = titers,
                            covInfo = covInfo,
                            estimatedParameters = estimatedParameters,
                            PoDModel = PoDModel,
                            ci = 0.95)
  # ggplot of the first Curve (value of binary covariate = 1)
  PoDCurvePlot1 <- ggplot(PoDCIdfs$PoDCICurve1df) +
    geom_line(aes(x = titers, y = lowerCI), linetype = 5) +
    geom_line(aes(x = titers, y = PointEst )) +
    geom_line(aes(x = titers, y = upperCI), linetype = 5) +
    ylab("PoD") +
    ggtitle(paste0("PoD curve (with fixed covariate value = ", covInfo[[1]]$values[1] ,") with ",ci*100,"% confidence intervals")) +
    themePlot
  
  # ggplot of the second Curve (value of binary covariate = 0)
  PoDCurvePlot2 <- ggplot(PoDCIdfs$PoDCICurve2df) +
    geom_line(aes(x = titers, y = lowerCI), linetype = 5) +
    geom_line(aes(x = titers, y = PointEst )) +
    geom_line(aes(x = titers, y = upperCI), linetype = 5) +
    ylab("PoD") +
    ggtitle(paste0("PoD curve (with fixed covarieate value = ", covInfo[[1]]$values[2] ,") with ",ci*100,"% confidence intervals")) +
    themePlot
  
  return(list(PoDCurvePlot1 = PoDCurvePlot1,
              PoDCurvePlot2 = PoDCurvePlot2))
}

#' @title PoD Confidence interval and median for specific titer values
#'
#' @description Supplementary function for the calculation of PoD curve with the confidence ribbon (of a required level). Input values are related to PoDBAY package structure.
#'
#' @param titers numeric vector: grid of titers at which the confidence ribbon should be calculated
#' @param covInfo list: name: character vector with the name of covariate, values: vector of factors with values of covariate, prob: vector of probability of each value
#' @param estimatedParameters estimatedParameters list: 2 named data frames needed ("results", "resultsPriorReset") needed: set of estimated PoD curve parameters, output of \code{PoDParamEstimation} function
#' @param PoDModel list: predictionFun: function call, vars: character vector with names of variables (var.titer and currently just 1 additional binary covariate), params: character vector with names of parameters 
#' @param ci numeric, required confidence level
#'
#' @return
#' PoD CI ribbon, median and point estimate of curve
#'
#' @usage
#' PoDCIInPoints(titers = c(1e-10,100),
#'               covInfo,
#'               estimatedParameters,
#'               PoDModel,
#'               ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data(estimatedParameters)
#' data(PoDModel)
#' data(covariateInfo)
#'
#' ## Exaple 1
#' PoDCIInPoints(titers = c(1e-10,100),
#'               covInfo = covariateInfo,
#'               estimatedParameters,
#'               PoDModel,
#'               ci = 0.95)
#'
#' @export
PoDCIInPoints <- function(titers = c(1e-10,100),
                          covInfo,
                          estimatedParameters,
                          PoDModel,
                          ci = 0.95){
  
  #for now works just for 1 attribute (e.i. truthBrazil$vaccinated$covariateInfo[[1]])
  
  if(is.list(covInfo)){
    covInfoDF <- data.frame(factor(covInfo[[1]]$values))
    colnames(covInfoDF) <- covInfo[[1]]$name
  }else{
    covInfoDF <- as.data.frame(covInfo)
  }
  
  #calculate the point estimate of parameters
  PoDParamsPointEst <- PoDParamPointEstimation(estimatedParameters$resultsPriorReset[PoDModel$params], PoDModel, titers)
  
  #calculate the PoD point estimate curve (median curve) for the first value (value == 1)
  PoDParamsPointEstPoDCurve1 <- PoD(titers = titers, 
                                    attributes = covInfoDF[1, , drop = FALSE], 
                                    PoDModel = PoDModel, 
                                    valueParams = PoDParamsPointEst)
  
  #calculate the PoD point estimate curve (median curve) for the second value (value == 0)
  PoDParamsPointEstPoDCurve2 <- PoD(titers = titers, 
                                    attributes = covInfoDF[2, , drop = FALSE], 
                                    PoDModel = PoDModel, 
                                    valueParams = PoDParamsPointEst)
  
  
  #prepare empty matrixes for the PoDs of individual titers -> will be filled in the later step
  PoDMatrixCurve1 <- matrix(0, ncol = length(titers), nrow = nrow(estimatedParameters$results[PoDModel$params]))
  PoDMatrixCurve2 <- matrix(0, ncol = length(titers), nrow = nrow(estimatedParameters$results[PoDModel$params]))
  
  for (i in 1:nrow(estimatedParameters$results[PoDModel$params])) {
    # calculate PoDs for all PoD curves in estimatedParameters$results for the first value (value == 1)
    PoDMatrixCurve1[i,] <- PoD(titers = titers, 
                               attributes = covInfoDF[1, , drop = FALSE], 
                               PoDModel = PoDModel, 
                               valueParams = estimatedParameters$results[PoDModel$params][i,])
    
    # calculate PoDs for all PoD curves in estimatedParameters$results for the second value (value == 0)
    PoDMatrixCurve2[i,] <- PoD(titers = titers, 
                               attributes = covInfoDF[2, , drop = FALSE], 
                               PoDModel = PoDModel, 
                               valueParams = estimatedParameters$results[PoDModel$params][i,])
  }
  
  # calculate CI of functional values
  PoDCIMatrixCurve1 <- matrix(0, ncol = 3, nrow = ncol(PoDMatrixCurve1))
  PoDCIMatrixCurve2 <- matrix(0, ncol = 3, nrow = ncol(PoDMatrixCurve2))
  
  for (i in 1:ncol(PoDMatrixCurve1)) {
    PoDCIMatrixCurve1[i,] <- unlist(PoDCI(PoDMatrixCurve1[,i], ci = ci))
    PoDCIMatrixCurve2[i,] <- unlist(PoDCI(PoDMatrixCurve2[,i], ci = ci))
  }
  
  # prepare data for ggplot
  PoDCICurve1df <- data.frame(PoDCIMatrixCurve1,PoDParamsPointEstPoDCurve1)
  colnames(PoDCICurve1df) <- c("lowerCI", "median", "upperCI", "PointEst")
  PoDCICurve1df$titers <- titers
  
  PoDCICurve2df <- data.frame(PoDCIMatrixCurve2,PoDParamsPointEstPoDCurve2)
  colnames(PoDCICurve2df) <- c("lowerCI", "median", "upperCI", "PointEst")
  PoDCICurve2df$titers <- titers
  
  return(list(
    PoDParamsPointEst = PoDParamsPointEst,
    PoDCICurve1df = PoDCICurve1df, 
    PoDCICurve2df = PoDCICurve2df))
}

#' @title PoD curve confidence ribbon
#'
#' @description Supplementary function for \code{PoDCurvePlot} function. Function calculates the confidence ribbon around the PoD curve. 
#' 
#' @param data numeric vector for which we the confidence intervals should be calculated
#' @param ci numeric: required confidence level
#'
#' @return
#' lower bound of CI
#' median value
#' upper bound of CI
#'
#' @usage
#' PoDCI(data, ci = 0.95)
#'
#' @examples
#' ## Data preparation
#' data(estimatedParameters)
#' data(PoDModel)
#' data(covariateInfo)
#'
#' ## Example 1
#' RibbonTest(titers = c(1e-10,100),
#'            covInfo = covariateInfo,
#'            estimatedParameters,
#'            PoDModel,
#'            ci = 0.95)
#'
#' @export
PoDCI <- function(data, ci = 0.95){
  lower <- quantile(data, (1-ci)/2, na.rm = TRUE)
  uppper <- quantile(data, ci + (1-ci)/2, na.rm  = TRUE )
  median <- median(data, na.rm  = TRUE)
  
  return(list(lower = lower,
              median = median,
              upper = uppper)
  )
}

#' @title Ribbon test
#'
#' @description Function for checking if the confidence ribbon does not allow for constant or increasing function. Input values are related to PoDBAY package structure.
#'
#' @param titers numeric vector: grid of titers at which the confidence ribbon should be calculated
#' @param covInfo list: name: character vector with the name of covariate, values: vector of factors with values of covariate, prob: vector of probability of each value
#' @param estimatedParameters estimatedParameters list: 2 named data frames needed ("results", "resultsPriorReset") needed: set of estimated PoD curve parameters, output of \code{PoDParamEstimation} function
#' @param PoDModel list: predictionFun: function call, vars: character vector with names of variables (var.titer and currently just 1 additional binary covariate), params: character vector with names of parameters 
#' @param ci numeric, required confidence level
#'
#' @return
#' TRUE/FALSE for both curves and total results (TRUE = test passed, FALSE = test failed)
#'
#' @usage
#' RibbonTest(titers = c(1e-10,100),
#'            covInfo = list(var.covariateBinary = list(name = "var.covariateBinary",
#'                                                        values = c(0, 1),
#'                                                        prob = c(0.5, 0.5))),
#'            estimatedParameters,
#'            PoDModel,
#'            ci = 0.95)
#'
#' @examples
#' ## TBD
#'
#' @export
RibbonTest <- function(titers = c(1e-10,100),
                       covInfo = list(var.covariateBinary = list(name = "var.covariateBinary",
                                                                 values = c(0, 1),
                                                                 prob = c(0.5, 0.5))),
                       estimatedParameters,
                       PoDModel,
                       ci = 0.95){
  
  titerCIs <- PoDCIInPoints(titers = c(1e-10,100),
                            covInfo,
                            estimatedParameters,
                            PoDModel,
                            ci = 0.95)
  
  curve1CITest <- titerCIs$PoDCICurve1df$lowerCI[1] > titerCIs$PoDCICurve1df$lowerCI[2] &
    titerCIs$PoDCICurve1df$upperCI[1] > titerCIs$PoDCICurve1df$upperCI[2]
  
  curve2CITest <- titerCIs$PoDCICurve2df$lowerCI[1] > titerCIs$PoDCICurve2df$lowerCI[2] &
    titerCIs$PoDCICurve2df$upperCI[1] > titerCIs$PoDCICurve2df$upperCI[2]
  
  #Is the test passed? (TRUE = YES, FALSE = NO)
  resultCITest <-  (curve1CITest + curve2CITest) == 2
  
  
  return(list(curve1CITest = curve1CITest,
              curve2CITest = curve2CITest,
              resultCITest = resultCITest))
}

#' @title Population generation
#'
#' @description Function generates the population class with population summary statistics.
#'
#' @param N numeric: number of subjects in the population
#' @param mean numeric: mean value of antibody titers defining the population
#' @param stdDev numeric: standard deviation of antibody titer defining the population
#' @param unknownDistribution logical: TRUE if we have unknown aspect in the antibody titer distribution - NOT USED
#' @param UDFunction function: function defining the unknown aspect in the antibody titer distribution - NOT USED
#'
#' @return
#' generated population: an \code{Population-class} object of population with all its characteristics defined in the input parameters
#'
#' @usage
#' generatePopulation(N, mean, stdDev, unknownDistribution = FALSE, UDFunction = NULL)
#'
#' @examples
#'
#' # Example 1 - empty population
#' population0 <- generatePopulation()
#'
#' # Example 2
#' population1 <- generatePopulation(N = 100,
#'                                   mean = 5,
#'                                   stdDev = 2)
#'
#' @export
generatePopulation <- function(N = 0, 
                               mean = NA_real_, 
                               stdDev = NA_real_, 
                               unknownDistribution = FALSE, 
                               UDFunction = NULL) {
  if (is.na(mean))          {mean <- NA_real_ }
  if (is.na(stdDev))        {stdDev <- NA_real_ }
  
  generated <- population$new()
  generated$N <- N
  generated$mean <- mean
  generated$stdDev <- stdDev
  generated$getTiters()
  generated$unknownDistribution <- unknownDistribution
  if (unknownDistribution) { generated$UDFunction <- UDFunction}
  
  return(generated)
}

#' @title Covariate data generation
#'
#' @description Function generates a data frame of subject level covariate data 
#' according to covariate informatin given in input.
#'
#' @param covariateInfo named list: information about covariate name, values and distribution in the population
#' @param Nrow numeric: target number of rows in the output data frame
#' @param vaccStatus string: "vaccinated" or "control" 
#
#' @return data frame of subject level covariate values
#'
#' @usage
#' covariateData(covariateInfo = NULL, 
#'               Nrow = NULL, 
#'               vaccStatus = NULL)
#'
#' @examples
#' ## Data preparation
#' data(covariateInfo)
#' 
#' ## Example 1
#' covariateData(covariateInfo = covariateInfo, 
#'               Nrow = 10, 
#'               vaccStatus = "vaccinated")
#'
#' @export
covariateData <- function(covariateInfo = NULL, 
                          Nrow = NULL, 
                          vaccStatus = NULL){
  
  if(!is.null(covariateInfo)){
    if("var.vaccStatus" %in% names(covariateInfo)){
      covariateInfoUsed <- covariateInfo[!(names(covariateInfo) == "var.vaccStatus")]
    }else{
      covariateInfoUsed <- covariateInfo
    }
    # Data frame will comprise of:
    # 1. User-defined categorical covariates
    for (j in seq(length(covariateInfoUsed))) {
      nameUtil <- covariateInfoUsed[[j]][["name"]]
      covValUtil <- covariateInfoUsed[[j]][["values"]]
      probUtil <- covariateInfoUsed[[j]][["prob"]]
      covColUtil <- sample(covValUtil,
                           size = Nrow,
                           replace = TRUE,
                           prob = probUtil)
      
      if(!is.numeric(covValUtil)){
        covColUtil <- factor(covColUtil,
                             covariateInfoUsed[[j]][["values"]])
      }
      
      if (j==1){
        df <- data.frame("dummy.name"=covColUtil)
      } else {
        df <- cbind(df,data.frame(
          "dummy.name" = covColUtil))
      }
      colnames(df)[j]<- nameUtil
    }
  }else{
    df <- NULL
  }
  
  if(!is.null(vaccStatus)){
    # 2. Vaccination status, by default
    nameUtil <- "var.vaccStatus"
    covValUtil <- vaccStatus
    probUtil <- 1
    covColUtil <- sample(covValUtil,
                         size = Nrow,
                         replace = TRUE,
                         prob = probUtil)
    
    covColUtil <- factor(covColUtil,
                         c("vaccinated","control"))
    
    if(!is.null(df)){
      df <- cbind(data.frame(
        "dummy.name" = covColUtil),df)
    }else{
      df <- data.frame(
        "dummy.name" = covColUtil)
    }
    
    colnames(df)[1]<- nameUtil
  }
  
  return(df)
}

#' @title Probability of disease calculation for base PoD model; works with covariate data.
#'
#' @description Function calculates probability of disease (PoD) corresponding to given titers according to base model (a sigmoid PoD curve with three parameters: et50, gamma, pmax).
#'
#' @param titer numeric vector: subject level titers
#' @param pmax numeric: maximum PoD
#' @param et50 numeric: titer values corresponding to pmax/2 value, PoD(et50) = pmax/2
#' @param gamma numeric: parameter proportional to -slope of the PoD curve
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return vector of PoDs
#'
#' @usage
#' PoDBaseModel(titer, pmax, et50, gamma, adjustTiters = FALSE, adjustFrom = 0, adjustTo = 0)
#'
#' @examples
#' PoDBaseModel(-10:10, 0.05, 5, 7, FALSE, log2(10), log2(5))
#'
#' @details PoD is calculated as: \deqn{ PoD = p_{max} \frac{ (\frac{et50}{titer})^{\gamma} }{ 1 + (\frac{et50}{titer})^{\gamma}}, \ for \ titers \ > 0}{ PoD = pmax * (et50/titer)^{\gamma} / (1+ (et50/titer)^{\gamma}, for titers > 0} and 
#'\deqn{ PoD = pmax, \ for \ titers \ <= 0}{PoD = pmax for titers <= 0}.
#'
#' @export
PoDBaseModel <- function(titer, pmax, et50, gamma, adjustTiters = FALSE, adjustFrom = 0, adjustTo = 0) {
  
  if (adjustTiters & (adjustFrom < adjustTo) ) {warning(paste("The input value for \"adjustFrom\" is lower than \"adjustTo\" "))}
  if (adjustTiters) {titer[titer < adjustFrom] <- adjustTo}
  
  probDisease <- ifelse( titer > 0, pmax - pmax / ( 1 + ( et50 / titer ) ^ gamma), pmax)
  return(probDisease)
}


#' @title Probability of disease calculation for general PoD model; works with covariate data.
#'
#' @description Function calculates probability of disease (PoD) corresponding to given titers according to given general PoD equation.
#'
#' @param titers numeric vector: subject level titers
#' @param attributes named data frame: subjectAttributes for a given population
#' @param PoDModel list: information on the PoD curve function
#' @param valueParams named list: values of the general model PoD curve parameters
#' @param pmaxVar string (for mapping): usually "param.pmax"
#' @param negTitersPoD NULL
#' @param adjustTiters boolean: set to TRUE if titer values should be adjusted, for details see \code{PoD} function
#' @param adjustFrom numeric: value specifying the detection limit, all values below the detection limit will be adjusted to adjustTo value
#' @param adjustTo numeric: value to which titers below the detection limit will be adjusted
#'
#' @return vector of PoDs
#'
#' @usage
#' PoD(titers, 
#'     attributes, 
#'     PoDModel, 
#'     valueParams, 
#'     adjustTiters = FALSE, 
#'     adjustFrom = log2(10), 
#'     adjustTo = log2(5))
#'
#' @examples
#' ## Example 1
#' data(vaccinated)
#' data(PoDModel)

#' titers <- seq(from = -10, to = 9.98, by = 0.02)
#' attributes <- vaccinated$subjectAttributes
#' valueParams <- PoDModel$valueParams
#' 
#' PoD(titers, 
#'     attributes, 
#'     PoDModel, 
#'     valueParams, 
#'     adjustTiters = FALSE)
#'
#' @importFrom dplyr %>%
#' @export
PoD <- function(titers, 
                attributes, 
                PoDModel, 
                valueParams, 
                adjustTiters = FALSE, 
                adjustFrom = log2(10), 
                adjustTo = log2(5)){
  
  
  #detection limit for the titer values
  if (adjustTiters & (adjustFrom < adjustTo) ) {warning(paste("The input value for \"adjustFrom\" is lower than \"adjustTo\" "))}
  if (adjustTiters) {titers[titers < adjustFrom] <- adjustTo}
  
  ### NOW WORKING: independent variables: vaccStatus, seroStatus, gender, covariateBinary (binary categorical variables)
  ## should be generalized, now the code assumes that user supplies these covariates in the truth definition
  dataAll <- cbind(titers, attributes)
  colnames(dataAll)[1] <- "var.titer"
  PoDdat <- dataAll
  
  if (!is.integer(unlist(attributes))) {
    # Convert categorical to dummy variables
    # For now, hard-coded for binary categorical variables
    for(i in colnames(attributes)){
      PoDdat[i] <- as.numeric(
        dataAll[[i]] == levels(dataAll[[i]])[1])
    }
  }
  
  # Select only user-defined attributes 
  PoDdat <- PoDdat %>% dplyr::select(PoDModel$vars)
  
  list <- valueParams
  l <- c(PoDdat,list)
  probDisease <- as.numeric(do.call(PoDModel$predictionFun,l))
  
  return(probDisease)
}
# AssignPoD - see RefClassPopulation


#' @title Diseased population extraction; works with covariate data.
#'
#' @description
#' Function extracts diseased population from vaccinated and control populations 
#' if populations have assigned disease status (for example using \code{ClinicalTrial} function).
#'
#' @param ... \code{Population-class} object: one or more populations with assigned disease status
#'
#' @return
#' diseased population in the form of \code{Population-class} object
#'
#' @usage
#' ExtractDiseased(...)
#'
#' @examples
#' ## Example 1
#' # Data preparation
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with CI
#' ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' # Extracting the disease cases
#' ExtractDiseased(vaccinated, control)
#'
#' @export
ExtractDiseased <- function(...) {
  diseased <- generatePopulation(0)
  args = list(...)
  numArgs <- length(args)
  if (numArgs == 1){
    populationX <- args[[1]]
  } else if (numArgs == 2) {
    population1 <- args[[1]]
    population2 <- args[[2]]
    populationX <- mergePopulations(population1, population2)
  }
  
  titers <- populationX$getDiseasedTiters()
  attributes <- populationX$getDiseasedAttributes()
  diseased$titers <- titers
  diseased$subjectAttributes <- attributes
  diseased$N <- length(titers)
  diseased$mean <- mean(titers)
  diseased$stdDev <- sd(titers)
  diseased$diseaseStatus <- rep(TRUE, length(titers))
  
  if (length(populationX$PoDs)>0){
    PoDs <- populationX$getDiseasedPoDs()
    diseased$PoDs <- PoDs
  }
  return(diseased)
}

#' @title Non-diseased population extraction; works with covariate data.
#'
#' @description
#' Function extracts non-diseased population from vaccinated and control populations
#'  if populations have assigned disease status (for example using \code{ClinicalTrial} function).
#'
#' @param ... \code{Population-class} object: one or more populations with assigned disease status
#'
#' @return
#' nondiseased population in the form of \code{Population-class} object
#'
#' @usage
#' ExtractNondiseased(...)
#'
#' @examples
#' ## Example 1
#' # Data preparation
#' data(vaccinated)
#' data(control)
#'
#' # Estimating the disease status and case-count efficacy with CI
#' ClinicalTrial(vaccinated, control, CI = 0.95)
#'
#' # Extracting the non-diseased subjects
#' ExtractNondiseased(vaccinated, control)
#'
#' @export
ExtractNondiseased <- function(...) {
  nondiseased <- generatePopulation(0)
  args = list(...)
  numArgs <- length(args)
  if (numArgs == 1){
    populationX <- args[[1]]
  } else if (numArgs == 2) {
    population1 <- args[[1]]
    population2 <- args[[2]]
    populationX <- mergePopulations(population1, population2)
  }
  
  titers <- populationX$getNondiseasedTiters()
  attributes <- populationX$getNondiseasedAttributes()
  nondiseased$titers <- titers
  nondiseased$subjectAttributes <- attributes
  nondiseased$N <- length(titers)
  nondiseased$mean <- mean(titers)
  nondiseased$stdDev <- sd(titers)
  nondiseased$diseaseStatus <- rep(FALSE, length(titers))
  
  if (length(populationX$PoDs)>0){
    PoDs <- populationX$getNondiseasedPoDs()
    nondiseased$PoDs <- PoDs
  }
  return(nondiseased)
}


#' @title Merge two populations
#'
#' @description Function merges two instances of \code{Population-class} object.
#'
#' @param population1 \code{Population-class} object: population to be merged
#' @param population2  \code{Population-class} object: population to be merged
#' 
#' @return
#' merged population in the form of \code{Population-class} object
#'
#' @usage
#' mergePopulations(population1, population2)
#'
#' @examples
#' # Data preparation
#' data(diseased)
#' data(nondiseased)
#'
#' # Merge populations
#' mergePopulations(diseased, nondiseased)
#'
#' @export
mergePopulations <- function(population1, 
                             population2) {
  
  merged <- generatePopulation(0)
  merged$titers <- c(population1$titers, population2$titers)
  merged$N <- length(merged$titers)
  merged$mean <- mean(merged$titers)
  merged$stdDev <- sd(merged$titers)
  merged$diseaseStatus <- c(population1$diseaseStatus, population2$diseaseStatus)
  merged$PoDs <- c(population1$PoDs, population2$PoDs)
  merged$subjectAttributes <- rbind(population1$subjectAttributes, population2$subjectAttributes)
  
  return(merged)
}


#' @title Sample population
#'
#' @description Function upsamples/downsamples with/without replacement population data in the form of \code{Population-class} object
#'
#' @param population1 \code{Population-class} object: population to be sampled
#' @param Ntarget numeric: target size of the new population
#' @param replace boolean: TRUE if sample with replacement
#' 
#' @return
#' sampled population in the form of \code{Population-class} object
#'
#' @usage
#' samplePopulation(population1, Ntarget, replace)
#'
#' @examples
#' # Data preparation
#' data(diseased)
#'
#' # Sample population
#' samplePopulation(diseased, 100, TRUE)
#'
#' @export
samplePopulation <- function(population1, 
                             Ntarget,
                             replace) {
  
  sampled <- generatePopulation(0)
  
  if (length(population1$PoDs)>0){
    df <- cbind(population1$titers, population1$diseaseStatus, population1$PoDs, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.diseaseStatus"
    colnames(df)[3] <- "y.PoD"
    sampledDf  <- df[sample(nrow(df),
                            Ntarget,
                            replace = replace),]
    sampledTiters <- sampledDf$var.titer
    sampledDf$var.titer <- NULL  
    sampledDS <- sampledDf$y.diseaseStatus
    sampledDf$y.diseaseStatus <- NULL
    sampledPoDs <- sampledDf$y.PoD
    sampledDf$y.PoD <- NULL
    sampledAttributes <- sampledDf
    
    sampled$PoDs <- sampledPoDs
    
  } else {
    df <- cbind(population1$titers, population1$diseaseStatus, population1$subjectAttributes)
    colnames(df)[1] <- "var.titer"
    colnames(df)[2] <- "y.diseaseStatus"
    sampledDf  <- df[sample(nrow(df),
                            Ntarget,
                            replace = replace),]
    sampledTiters <- sampledDf$var.titer
    sampledDf$var.titer <- NULL  
    sampledDS <- sampledDf$y.diseaseStatus
    sampledDf$y.diseaseStatus <- NULL
    sampledAttributes <- sampledDf
  }
  
  sampled$titers <- sampledTiters
  sampled$N <- length(sampled$titers)
  sampled$mean <- mean(sampled$titers)
  sampled$stdDev <- sd(sampled$titers)
  sampled$diseaseStatus <- sampledDS
  sampled$subjectAttributes <-sampledAttributes
  
  return(sampled)
}

#' @title Make model 
#'
#' @description Functions takes information about PoD curve model 
#' and returns it in a unified format to be used throughout the package.
#'
#' @param bodyString string: the general PoD model function
#' @param functionVars string vector: variable names; should match bodyString
#' @param rangeVars named list: range for each variable
#' @param functionParams lstring vector: parameter names; should match bodyString
#' @param rangeParams named list: range for each parameter
#' @param valueParams named list: value for each parameter
#'
#' @return 
#' PoD model: list(predictionFun, vars, params, rangeParams, rangeVars, valueParams)
#'
#' @usage
#' MakeModelObj(bodyString,
#'              functionVars,
#'              rangeVars,
#'              functionParams,
#'              rangeParams,
#'              valueParams = NA)
#'
#' @examples
#' # Example 1
#' variable <- list()
#' parameter <- list()
#' PoDfun <- "I(param.pmax*(param.pmax.covariateBinary^var.covariateBinary) /(
#'   1+(var.titer/(param.et50*(param.et50.covariateBinary^var.covariateBinary))
#'   )^(param.gamma*(param.gamma.covariateBinary^var.covariateBinary) )))"
#' 
#' variable$name <- c("var.titer", "var.covariateBinary")
#' variable$range <- list(var.titer = c(0,Inf),
#'                        var.covariateBinary = c(0,1))
#' parameter$name <- c("param.pmax",
#'                     "param.pmax.covariateBinary",
#'                     "param.et50",
#'                     "param.et50.covariateBinary",
#'                     "param.gamma",
#'                     "param.gamma.covariateBinary")
#' parameter$range <- list("param.pmax" = c(1e-6, 0.1),  
#'                         "param.pmax.covariateBinary" = c(1e-6, 10), 
#'                         "param.et50" = c(1e-6, 100),
#'                         "param.et50.covariateBinary" = c(1e-3, 10), 
#'                         "param.gamma" = c(0, 100), 
#'                         "param.gamma.covariateBinary" = c(1e-6, 10)) 
#' parameter$value <- list(param.pmax = 0.05,
#'                         param.pmax.covariateBinary = 1,
#'                         param.et50 = 5,
#'                         param.et50.covariateBinary = 2, 
#'                         param.gamma = 7,
#'                         param.gamma.covariateBinary = 1)
#' PoDModel <- MakeModelObj(PoDfun,
#'                          variable$name,
#'                          variable$range,
#'                          parameter$name,
#'                          parameter$range,
#'                          parameter$value)
#'
#' @export
# This function makes a model object to be used in fit function
MakeModelObj <- function(bodyString,
                         functionVars,
                         rangeVars,
                         functionParams,
                         rangeParams,
                         valueParams = NA) {
  
  # Make arguments list for fit function
  allArgs <- c(functionVars,functionParams)
  argStringUtil <- "args <- alist("
  stringUtilIndex <- 1
  for (stringUtil in allArgs) {
    if (stringUtilIndex==length(allArgs)) {
      # note that each argument is given a default value of NULL
      argStringUtil <- paste0(argStringUtil,stringUtil,"=NULL)")
    } else {
      argStringUtil <- paste0(argStringUtil,stringUtil,"=NULL,")
    }
    stringUtilIndex <- stringUtilIndex + 1
  }
  
  # define the args and body based on the respective strings
  eval(parse(text=argStringUtil))
  eval(parse(text=paste0("body <- quote(",bodyString,")")))
  
  # name the function f
  f <- function() {}
  # assign argStringUtil as the arguments list of the function
  formals(f) <- args
  # Assign the bodyString as the body of the function
  body(f) <- body
  # Assign standard environment
  environment(f) <- parent.frame()
  
  # return function, variable names, parameter names, and the parameter ranges (used in fit function)
  return(list(predictionFun = f,
              vars = functionVars,
              params = functionParams,
              rangeParams = rangeParams,
              rangeVars = rangeVars,
              valueParams = valueParams))
}

#' @title Numeric to boolean
#'
#' @description  Converts numeric format to boolean format.
#'
#' @details If the function is supposed to be used on a vector, the form \code{sapply("vector", numToBool)} needs to be applied.
#'
#' @param x numeric value (0, 1)
#'
#' @return boolean value (T, F)
#'
#' @usage
#' numToBool(x)
#'
#' @examples
#' dStatus <- c(0,0,1,1,0,1)
#' sapply(dStatus, numToBool)
#'
#' @export
numToBool <- function(x) {if (x != 0) { return(T) } else { return(F) }}

#' @title Error message
#'
#' @description  Error message: the input value for "name" is incorrect
#'
#' @param name name of the input value
#'
#' @return error message: "the input value for "name" is incorrect"
#'
#' @export
incorrectInput <- function(name) {
  stop(paste("the input value for", name, "is incorrect"))
}

#' @title Population error message
#'
#' @description  Error meassage: the input value for "name" is incorrect.
#'
#' @param name name of the input value
#'
#' @return error message: "the input value for "name" is incorrect"
#'
#' @export
incorrectPopulationInput <- function(name) {
  stop(paste("Input value for", name," is incorrect. Input needs to be a population class. For more details see 'population' vignette - vignette(\"population\", package = \"PoDBA\")"))
}

#' @importFrom methods is new
#' @importFrom stats dnorm integrate median optim quantile rbinom rnorm sd uniroot 
#' @importFrom ggplot2 aes element_blank element_line element_rect element_text geom_line ggplot ggtitle theme ylab
NULL

#' @title Population class; works with covariate data.
#'
#' @description Population reference class which provides summary and subject level information about the population
#'
#' @field N numeric: number of subjects in the population
#' @field mean numeric: mean value of titers 
#' @field stdDev numeric: standard deviation of titers 
#' @field unknownDistribution logical: TRUE if titer distribution is not normally /log-normally distributed; FALSE titer distribution function needs to be defined by user
#' @field UDFunction function: user-defined titer distribution
#' @field titers numeric: subject level titers, generated with \code{getTiters} method
#' @field PoDs numeric: subject level probability of disease (PoD), generated with \code{assginPoD} method
#' @field diseaseStatus logical: subject level disease status (TRUE if diseased), generated with \code{ClinicilaTrial} function
#' @field userDefinedTiterDist logical: TRUE if titer distribution is not normally /log-normally distributed; FALSE titer distribution function needs to be defined by user
#' @field titerDist function: user-defined titer distribution
#' @field subjectAttributes data.frame: subject level covariate/attribute data, generated with \code{addSubjectAttributes}.
#'
#' @exportClass Population
population <- setRefClass(
  "Population",
  fields = list(
    N = "numeric",
    mean = "numeric",
    stdDev = "numeric",
    unknownDistribution = "logical",
    UDFunction = "function", # function with one parameter that defines the length of generated sequence
    titers = "numeric", # should be generated automatically
    PoDs = "numeric", # probability of disease for each subject
    diseaseStatus = "logical",
    userDefinedTiterDist = "logical",
    titerDist = "function",
    subjectAttributes = "data.frame"
  ),
  methods = list(
    initialize = function() {
      mean <<- NA_real_
      stdDev <<- NA_real_
      unknownDistribution <<- FALSE
      userDefinedTiterDist <<- FALSE
    },
    getTiters = function() {
      if (length(titers) == 0) {
        titers <<- rnorm(N, mean, stdDev)
      }
      return(titers)
    },
    popFun = function() {
      if (unknownDistribution) {
        toRet <- approxfun(
          density(rnorm(1e6, mean = mean, sd = stdDev) + getUnknown(1e6)),
          yleft = 0,
          yright = 0)
        return(toRet)
      }
      else {
        toRet <- function(x) {
          dnorm(x, mean = mean, sd = stdDev)
        }
        return(toRet)
      }
    },
    popX = function() {
      if (unknownDistribution) {
        return(titers + getUnknown(N))
      } else {
        return(titers)
      }
    },
    getUnknown = function(n) {
      return(UDFunction(n))
    },
    assignPoD = function(x){
      PoDs <<- x
    },
    addSubjectAttributes = function(df){
      subjectAttributes <<- df
    },
    getDiseasedCount = function() {
      return(length(which(diseaseStatus)))
    },
    getNondiseasedCount = function() {
      return(length(which(!diseaseStatus)))
    },
    getDiseasedTiters = function() {
      return(titers[which(diseaseStatus)])
    },
    getNondiseasedTiters = function() {
      return(titers[which(!diseaseStatus)])
    },
    getDiseasedAttributes = function(){
      return(filter(subjectAttributes,diseaseStatus))
    },
    getNondiseasedAttributes = function(){
      return(filter(subjectAttributes,!diseaseStatus))
    },
    getDiseasedPoDs = function(){
      return(PoDs[which(diseaseStatus)])
    },
    getNondiseasedPoDs = function(){
      return(PoDs[which(!diseaseStatus)])
    },
    getSubAttributeKeys = function() {
      return(names(subjectAttributes))
    }
  )
)

#' @title Subject level titers
#'
#' @name getTiters
#'
#' @description Returns subject level titers. If titers are not yet generated, the function generates them based on \code{Population-class} object attributes: N, mean, stdDev.
#'
#' @details Inputs into the function (N, mean, stdDev) are taken from the \code{Population-class} object attributes.
#'
#' @return Subject level titers
NULL
population$methods(
  getTiters = function() {
    if (!userDefinedTiterDist) {
      if (length(titers) == 0) {
        titers <<- rnorm(N, mean, stdDev)
      }
    } else {
      if (length(titers) == 0) {
        titers <<- SamplePDF(PDF=titerDist,sampleSize=N)
      }
    }
    return(titers)
  }
)

#' @title Population function
#'
#' @name popFun
#'
#' @description Function describing the titer distribution of the population: mean, standard deviation and an additional unknown factor affecting the shape of the distribution (e.g. mixture of two normals or other shapes defined by user).
#'
#' @details Inputs into the function (mean, stdDev, Unknowndistribution) and getUnknown method are taken from the \code{Population-class} object.
#'
#' @return Titer distribution function
NULL
population$methods(
  popFun = function() {
    if (unknownDistribution) {
      toRet <- approxfun(
        density(rnorm(1e6, mean = mean, sd = stdDev) + getUnknown(1e6)),
        yleft = 0,
        yright = 0)
      return(toRet)
    }
    else {
      toRet <- function(x) {
        dnorm(x, mean = mean, sd = stdDev)
      }
      return(toRet)
    }
  }
)

#' @title Add noise (e.g. measurement error) to population titers
#'
#' @name popX
#'
#' @description Function adds noise to population titers accounting for an unknown factor affecting the titer distibution.
#'
#' @details Inputs into the function: N, unknownDistribution and getUnknown() method are taken from the \code{Population-class} object.
#'
#' @return subject level titers
NULL
population$methods(
  popX = function() {
    if (unknownDistribution) {
      return(titers + getUnknown(N))
    } else {
      return(titers)
    }
  }
)

#' @title Generate unknown
#'
#' @name getUnknown
#'
#' @description Function generates unknown part of the titers which is eventually added to the original titers in \code{popX} and to the original titer distribution in \code{popFun}.
#'
#' @param n numeric: number of subjects in the population
#'
#' @details Input into the function: UDFunction is taken from the \code{Population-class} object. UDFunction is used for generating the unknown part of the titer distribution.
#'
#' @return unknown part of the titers
NULL
population$methods(
  getUnknown = function(n) {
    return(UDFunction(n))
  }
)

#' @title Assign probability of disease (PoD)
#'
#' @name assignPoD
#'
#' @description Function assigns subject-level probability of disease based on PoD curve and subject level titer.
#'
#' @param x numeric vector - vector of estimated PoD values
#'
#' @details The input into the function is either calculated using \code{PoD} function or if the PoD curve is unknown the same arbitrary PoD can be assigned to the whole population.
#'
#' @return Subject level probability of disease for the population
NULL
population$methods(
  assignPoD = function(x){
    PoDs <<- x
  }
)

#' @title Add to the population a dataframe with subject attributes (covariates)
#'
#' @name addSubjectAttributes
#'
#' @description Function adds a dataframe with subject-level covariate values.
#' 
#' @param data.frame   dataframe with subject-level covariate values
#'
#' @details The input into the function is calculated using \code{covariateData()} function.
NULL
population$methods(
  addSubjectAttributes = function(df){
    subjectAttributes <<- df
  }
)

#' @title Diseased count
#'
#' @name getDiseasedCount
#'
#' @description Function calculates the number of diseased subjects (disease status = TRUE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric: number of the diseased subjects in the \code{Population-class} object
NULL
getDiseasedCount <- function() {
  return(length(which(diseaseStatus)))
}

#' @title Non-diseased count
#'
#' @name getNondiseasedCount
#'
#' @description Function calculates the number of non-diseased subjects (disease status = FALSE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric: number of the non-diseased subjects in the \code{Population-class} object
NULL
getNondiseasedCount <- function() {
  return(length(which(!diseaseStatus)))
}

#' @title Diseased titers
#'
#' @name getDiseasedTiters
#'
#' @description Function returns titers of diseased subjects (disease status = TRUE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric vector: titers of diseased subjects in the \code{Population-class} object
NULL
getDiseasedTiters <- function() {
  return(titers[which(diseaseStatus)])
}

#' @title Non-diseased titers
#'
#' @name getNondiseasedTiters
#'
#' @description Function returns titers of non-diseased subjects (disease status = FALSE) in the \code{Population-class} object.
#'
#' @details Input into the function, "diseaseStatus", is taken from the \code{Population-class} object attribute. Information about disease status is written into the \code{Population-class} object by the \code{ClinicalTrial()} function.
#'
#' @return numeric vector: titers of non-diseased subjects in the \code{Population-class} object
NULL
getNondiseasedTiters <- function() {
  return(titers[which(!diseaseStatus)])
}

#' @title Initialize Population fields based on subjectAttributes field
#'
#' @name initFromSubjectAttributes
#'
#' @description Uses the user provided subjectAttributes to complete other fields
#'
#' @return Population
NULL
population$methods(
  initFromSubjectAttributes = function() {
    titers <<- subjectAttributes$titer
    diseaseStatus <<- as.logical(subjectAttributes$diseaseStatus)
    PoDs <<- as.numeric(subjectAttributes$PoD)
    N <<- length(subjectAttributes$titer)
    mean <<- mean(subjectAttributes$titer)
    stdDev <<- sd(subjectAttributes$titer)
  }
)

#' @title Initialize Population fields for vaccinated based on subjectAttributes field
#'
#' @name initVaccFromSubjectAttributes
#'
#' @description Uses the user provided subjectAttributes to complete other fields
#'
#' @return Population
NULL
population$methods(
  initVaccFromSubjectAttributes = function() {
    vaccLoc <- subjectAttributes$vaccStatus==1
    titers <<- subjectAttributes$titer[vaccLoc]
    diseaseStatus <<- as.logical(subjectAttributes$diseaseStatus[vaccLoc])
    PoDs <<- as.numeric(subjectAttributes$PoD[vaccLoc])
    N <<- length(subjectAttributes$titer[vaccLoc])
    mean <<- mean(subjectAttributes$titer[vaccLoc])
    stdDev <<- sd(subjectAttributes$titer[vaccLoc])
    subjectAttributes <<- subjectAttributes[vaccLoc,]
  }
)

#' @title Initialize Population fields for control based on subjectAttributes field
#'
#' @name initContFromSubjectAttributes
#'
#' @description Uses the user provided subjectAttributes to complete other fields
#'
#' @return Population
NULL
population$methods(
  initContFromSubjectAttributes = function() {
    contLoc <- subjectAttributes$vaccStatus==0
    titers <<- subjectAttributes$titer[contLoc]
    diseaseStatus <<- as.logical(subjectAttributes$diseaseStatus[contLoc])
    PoDs <<- as.numeric(subjectAttributes$PoD[contLoc])
    N <<- length(subjectAttributes$titer[contLoc])
    mean <<- mean(subjectAttributes$titer[contLoc])
    stdDev <<- sd(subjectAttributes$titer[contLoc])
    subjectAttributes <<- subjectAttributes[contLoc,]
  }
)

#' @title Check to see if titers are consistent
#'
#' @name checkTiters
#'
#' @description Checks to see if the user supplied titers match the generated titers
#'
#' @return String with test result
NULL
population$methods(
  checkTiters = function() {
    if (isTRUE(all.equal(subjectAttributes$titer, titers))) {
      print("Titers in subjectAttributes are consistent with those in the field 'titers'")
    } else {
      warning("Titers in subjectAttributes are NOT consistent with those in the field 'titers'")
    }
  }
)