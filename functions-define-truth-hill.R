# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

defineTruth <- function() {
  variable <- list()
  parameter <- list()
  vaccinated <- list()
  control <- list()
  
  ## Name ----
  name <- "results-logistic-effect30-n15000"
  
  ## PoD model for data generation ----
  PoD <- "I(1/(1+exp(-(param.beta*(param.beta.covariateBinary^var.covariateBinary)*var.titer + param.delta))))"
  
  # Independent variables
  variable$name <- c("var.titer", "var.covariateBinary")
  
  # Parameters
  parameter$name <- c("param.beta",
                      "param.beta.covariateBinary",
                      "param.delta")
  parameter$value <- list(param.beta = -0.33,
                          param.beta.covariateBinary = 0.442, # 0.753; 0.442 # for 10%; 30% effect
                          param.delta = -2)
  
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