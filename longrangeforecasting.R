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
first.n.h = 150
last.n.h = 40
nu = c(0.4,1.0,1.0)
lambda.r = 0.1
m = 4
alpha = 0.5
reduced.units = 10

#Fixed parameters
tau = 1
layers = 3
pi.w = rep(0.1, layers)
pi.win = rep(0.1, layers)
eta.w = rep(0.1, layers)
eta.win = rep(0.1, layers)
trainLen = 400
testLen = 1
future = 1
forward = 10
locations = 10
iterations = 100
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
  
  
  n.h = c(rep(first.n.h, layers-1), last.n.h)
  #Begin DESN forecasting
  testing = deep.esn(y.train = y.train,
                     x.insamp = designMatrix,
                     x.outsamp = designMatrixOutSample,
                     y.test = sets$yTest,
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
                     future = testLen,
                     layers = layers,
                     reduced.units = reduced.units,
                     startvalues = NULL,
                     activation = 'tanh',
                     distribution = 'Normal',
                     scale.factor = y.scale,
                     scale.matrix = addScaleMat,
                     logNorm = FALSE,
                     fork = FALSE,
                     parallel = FALSE,
                     verbose = TRUE)
  
  #Save predictions for each ensemble iteration
  mean.pred[,,f] = testing$predictions
  
  
  
  #Append mean of ensembles to training data
  new.train = rbind(new.train, testing$forecastmean)
  newRaw = rbind(newRaw, testing$forecastmean)
  
  #Update all values for next future point
  trainLen = trainLen + 1
  testindex = testindex + 1
  
  #print progress
  #print(f)
}


#Forecast averages & MSE
forc.mean = t(sapply(1:locations, function(x) colMeans(mean.pred[x,,])))
mse = sum((t(sets$yTest)-forc.mean)^2)/(forward*locations); mse



