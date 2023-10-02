# 2023ESNb
Supplemental codes for "Calibrated Forecasts of Quasi-Periodic Climate Processes with Deep Echo State Networks and Penalized Quantile Regression" by Matthew Bonas, Christopher K. Wikle and Stefano Castruccio

## Data
A dataset "SimulatedData.RData" with 10 variables (locations) and 500 time points. This data is to be used in conjunction with the R scripts.

## functions.R
R script containing the user created functions used in both longrangeforecasting.R and calibration.R. This script does not need to be run manually for the other files automatically import the functions from this file.

## longrangeforecasting.R
R script containing code and methods used to generate long-range forecasts with the Echo State Network (ESN) using SimulatedData.RData.

## calibration.R
R script containing the code and methods used to calibrate the forecasts from the ESN using penalized quantile regression. 
