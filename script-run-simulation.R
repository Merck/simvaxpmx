# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Elucidating vaccine efficacy using a correlate of protection, demographics, and logistic regression"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

rm(list=ls())

# Load required packages
library(Rcpp)
library(dplyr)
library(tidyr)
library(pracma)
library(modelr)
library(rlang)
library(ggridges)
library(tidyverse)
library(MASS)
library(vaxpmx) # version 0.0.3

# Source simulation function
source("functions-simulation.R")
source("functions-simulation-supplementary.R")

# Load the truth for the simulation
source("functions-define-truth-logistic.R")
#source("functions-define-truth-hill.R") 

experimentName <- defineTruth()$name

# Run clinical trial simulation
result <- Simulation()
