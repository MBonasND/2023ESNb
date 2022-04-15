###########################################
###########################################
### Calibration of Long-Range Forecasts ###
###########################################
###########################################


#clear enviroment and load libraries
rm(list = ls())
library(tidyverse)
library(Matrix)
library(abind)
library(scam)

#parallel libraries
library(doParallel)
library(parallel)
library(foreach)

#load functions and data
source('functions.R')
load('SimulatedData.RData')

#specify cores for parallel
#options(cores = 4)


########################################
### 'Windows' Long-Range Forecasting ###
########################################

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
tau = 1
start.range = 360
testLen = 1
future = 1
forward = 10
locations = 10
iterations = 100
rawData = sim.dat

n.w = 5 #user specified number of 'windows'
WindowForcs = array(NaN, dim = c(locations, iterations, 1))


#### code below can be computationally expensive - DO NOT RUN ####
#### speed of code below can be increased by running parallel = T in ESN function ###
#### parallel = T will run 'iter' variable in parallel across cores in the machine ###


for(w in 1:n.w)
{
  
  trainLen = (w-1)*forward + start.range
  
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
                           verbose = T,
                           parallel = F)
    
    #Save predictions for each ensemble iteration
    mean.pred[,,f] = testing$predictions
    
    
    
    #Append mean of ensembles to training data
    new.train = rbind(new.train, testing$forecastmean)
    newRaw = rbind(newRaw, testing$forecastmean)
    
    #Update all values for next future point
    trainLen = trainLen + 1
    testindex = testindex + 1
    
  }
  
  WindowForcs = abind(WindowForcs, mean.pred, along = 3)
  
}
WindowForcs = WindowForcs[,,-1]

###################################################################
### Determine Optimum SD Vector for Each Location - Algorithm 1 ###
###################################################################

locations = 10
n.w = 5
rawData = sim.dat
tau = 1
start.range = 360
horizon = 10
true.range = (start.range + tau + 1):(start.range + n.w*horizon + tau)

optim.sd.mat = matrix(NaN, nrow = locations, ncol = horizon) #
for(l in 1:locations)
{
  location = l
  dat = WindowForcs[location,,]
  true.y = rawData[true.range, location]
  means = apply(dat, 2, mean)
  
  #Generate data windows for j-step ahead forecasts
  true.window = list()
  mean.window = list()
  for(i in 1:horizon)
  {
    index = seq(i, (n.w-1)*horizon+i, horizon)
    true.window[[i]] = true.y[index]
    mean.window[[i]] = apply(dat[,index], 2, mean)
  }
  
  
  ######################################
  ### Optimal SD w/ Monotonic Spline ###
  ######################################
  
  #Generate optimal SD
  optim.sd = rep(0, horizon)
  for(i in 1:horizon)
  {
    optim.sd[i] = sd(true.window[[i]] - mean.window[[i]])
  }
  
  #Monotonic Spline
  testy = optim.sd
  testx = 1:horizon
  fit = scam(testy~s(testx, k=-1, bs="mpi"), 
             family=gaussian(link="identity"))
  
  #set optimimum values from monotonic spline
  optim.sd.mat[l,] = fit$fitted.values
}
