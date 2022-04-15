#####################################
#####################################
### Long-Range Forecasting w/ ESN ###
#####################################
#####################################

#clear enviroment and load libraries
rm(list = ls())
library(tidyverse)
library(Matrix)
library(abind)


#parallel libraries
library(doParallel)
library(parallel)
library(foreach)

#load functions and data
source('functions.R')
load('SimulatedData.RData')

#specify cores for parallel
options(cores = 4)


##############################
### Long-Range Forecasting ###
##############################

#Parameter specification
#User should perform a CV to determine parameter values
n.h = 60
nu = 0.55
lambda.r = 0.001
m = 4
alpha = 0.0023

#Fixed parameters
pi.w = 0.1
eta.w = 0.1 #only needed if distribution = 'Unif'
pi.win = 0.1
eta.win = 0.1 #only needed if distribution = 'Unif'
iterations = 100
tau = 1
trainLen = 400
testLen = 1
future = 1
forward = 20
locations = 10
rawData = sim.dat

#Create training and testing sets
sets = cttv(rawData, tau, trainLen, forward)
new.train = sets$yTrain
testindex = sets$xTestIndex[1]
newRaw = rawData[1:testindex,]

#Preallocate empty results matrix
mean.pred = array(NaN, dim = c(locations, iterations, forward))


#Begin future forecasts
for(f in 1:forward)
{
  #Generating input data
  input.dat = gen.input.data(trainLen = trainLen,
                             m = m,
                             tau = tau,
                             yTrain = new.train,
                             rawData = newRaw,
                             locations = locations,
                             xTestIndex = testindex,
                             testLen = testLen)
  y.scale = input.dat$y.scale
  y.train = input.dat$in.sample.y
  designMatrix = input.dat$designMatrix
  designMatrixOutSample = input.dat$designMatrixOutSample
  addScaleMat = input.dat$addScaleMat
  
  
  #Begin ESN forecasting
  testing = ensemble.esn(y.train = y.train,
                         x.insamp = designMatrix,
                         x.outsamp = designMatrixOutSample,
                         y.test =  NULL,
                         n.h = n.h,
                         nu = nu,
                         pi.w = pi.w,
                         pi.win = pi.win,
                         eta.w = eta.w,
                         eta.win = eta.win,
                         lambda.r = lambda.r,
                         alpha = alpha,
                         m = m,
                         iter = iterations,
                         future = future,
                         startvalues = NULL,
                         activation = 'tanh',
                         distribution = 'Normal',
                         polynomial = 1,
                         scale.factor = y.scale,
                         scale.matrix = addScaleMat,
                         verbose = F,
                         parallel = T)
  
  #Save predictions for each ensemble iteration
  mean.pred[,,f] = testing$predictions
  
  
  
  #Append mean of ensembles to training data
  new.train = rbind(new.train, testing$forecastmean)
  newRaw = rbind(newRaw, testing$forecastmean)
  
  #Update all values for next future point
  trainLen = trainLen + 1
  testindex = testindex + 1
  
  #print progress
  print(f)
}


#Forecast averages & MSE
forc.mean = t(sapply(1:locations, function(x) colMeans(mean.pred[x,,])))
mse = sum((t(sets$yTest)-forc.mean)^2)/(forward*locations); mse


