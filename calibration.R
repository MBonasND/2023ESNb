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
options(cores = 4)




#Calibration Functions
qs_cont = function (x.inter = 10, p.inter = 10, degree = 3, logit = FALSE, 
                    order = 2, kappa = 0, n.cyc = 100, c.crit = 1e-05, plot = TRUE, 
                    power = NULL, ...) 
{
  if (x.inter <= 0) {
    warning("the value of x.inter supplied is less than 0, the value of 10 was used instead")
    x.inter <- 10
  }
  if (p.inter <= 0) {
    warning("the value of p.inter supplied is less than 0, the value of 10 was used instead")
    p.inter <- 10
  }
  if (degree <= 0) {
    warning("the value of degree supplied is less than zero or negative the default value of 3 was used instead")
    degree <- 3
  }
  if (order < 0) {
    warning("the value of order supplied is zero or negative the default value of 2 was used instead")
    order <- 2
  }
  if (kappa < 0) {
    warning("the value of kapa supplied is less than 0, the value of zero was used instead")
    kappa <- 0
  }
  if (n.cyc < 0) {
    warning("the value of n.cyc is less than zero the default value of 100 was used instead")
    n.cyc <- 100
  }
  if (c.crit <= 0) {
    warning("the value of c.crit is less or equal than zero the default value of 1e-5 was used instead")
    c.crit <- 1e-05
  }
  out <- list(x.inter = x.inter, p.inter = p.inter, degree = degree, 
              logit = as.logical(logit)[1], order = order, kappa = kappa, 
              n.cyc = n.cyc, c.crit = c.crit, plot = as.logical(plot)[1], 
              power = power)
}




qs_mod = function (y, x, x.lambda = 1, p.lambda = 1, data = NULL, cent = 100 * 
                     pnorm((-4:4) * 2/3), xsi = 1, totpoints, control = qs_cont(...), 
                   print = TRUE, ...) 
{
  tpower <- function(x, t, p = 1) (x - t)^p * (x > t)
  bbase <- function(x, xl = min(x), xr = max(x), ndx = 10, 
                    deg = 3) {
    dx <- (xr - xl)/ndx
    kts <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, kts, FUN = tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * 
                                          dx^deg)
    B <- (-1)^(deg + 1) * P %*% t(D)
    B
  }
  rowtens = function(X) {
    one <- matrix(1, nrow = 1, ncol = ncol(X))
    kronecker(X, one) * kronecker(one, X)
  }
  ptrans <- function(x, p){ if (abs(p) <= 1e-04){ 
    log(x)
  }else{ I(x^p)}}
  invptrans <- function(x, p){ if (abs(p) <= 1e-04){
    exp(x)
  }else{ x^(1/p)}}
  scall <- deparse(sys.call())
  ylab <- deparse(substitute(y))
  xlab <- deparse(substitute(x))
  y <- if (!is.null(data)) {
    get(deparse(substitute(y)), envir = as.environment(data))
  }else{ y}
  x <- if (!is.null(data)){ 
    get(deparse(substitute(x)), envir = as.environment(data))
  }else{ x}
  if (!is.null(control$power)) {
    ox <- x
    x <- ptrans(x, control$power)
  }
  m <- length(x)
  xl <- min(x)
  xr <- max(x)
  nsegx <- control$x.inter
  nsegp <- control$p.inter
  bdeg <- control$degree
  p <- cent/100
  n <- length(p)
  Bx <- bbase(x, xl, xr, nsegx, bdeg)
  if (control$logit) {
    logitp <- log(p/(1 - p))
    Bp <- bbase(logitp, -20, 20, nsegp, bdeg)
  }else {
    Bp <- bbase(p, 0, 1, nsegp, bdeg)
  }
  nbx <- ncol(Bx)
  nbp <- ncol(Bp)
  Tx <- rowtens(Bx)
  Tp <- rowtens(Bp)
  Dx <- diff(diag((1:nbx)/xsi), diff = control$order) 
  Dp <- diff(diag((1:nbp)/xsi), diff = control$order)
  Px <- x.lambda * t(Dx) %*% Dx
  Pp <- p.lambda * t(Dp) %*% Dp
  P <- kronecker(Pp, diag(nbx)) + kronecker(diag(nbp), Px)
  P <- P + control$kappa * diag(nrow(P))
  Y <- outer(y, rep(1, n))
  Z <- 0 * Y + mean(Y)
  OP <- outer(rep(1, m), p)
  b <- 0.001
  for (it in 1:control$n.cyc) {
    R <- Y - Z
    W <- ifelse(R > 0, OP, 1 - OP)/sqrt(b + R^2)
    Q <- t(Tx) %*% W %*% Tp
    dim(Q) <- c(nbx, nbx, nbp, nbp)
    Q <- aperm(Q, c(1, 3, 2, 4))
    dim(Q) <- c(nbx * nbp, nbx * nbp)
    r <- t(Bx) %*% (Y * W) %*% Bp
    dim(r) <- c(nbx * nbp, 1)
    A <- solve(Q + P, r)
    dim(A) <- c(nbx, nbp)
    Znew <- Bx %*% A %*% t(Bp)
    dz <- sum(abs(Z - Znew))
    if (dz < control$c.crit) 
      break
    Z <- Znew
  }
  xg <- seq(xl, xr, length = totpoints)
  Bg <- bbase(xg, xl, xr, nsegx, bdeg)
  Zg <- Bg %*% A %*% t(Bp)
  if (!is.null(control$power)) {
    x <- ox
    xg <- invptrans(xg, control$power)
  }
  if (control$plot) {
    plot(x, y, pch = 15, cex = 0.5, col = gray(0.7), ylab = ylab, 
         xlab = xlab)
    matlines(x[order(x)], Z[order(x), ], type = "l", 
             lty = 1, lwd = 1)
  }
  colnames(Zg) <- as.character(round(cent, 2))
  per <- rep(0, length(cent))
  quantFun <- list()
  for (i in 1:length(cent)) {
    quantFun[[i]] <- splinefun(xg, Zg[, i], method = "natural")
    ll <- quantFun[[i]](x)
    per[i] <- (1 - sum(y > ll)/length(y)) * 100
    if (print) 
      cat("% of cases below ", cent[i], "centile is ", 
          per[i], "\n")
  }
  names(quantFun) <- namesFun <- as.character(round(cent, 2))
  out <- list(y = y, x = x, knots = xg, fitted.values = Zg, 
              cent = cent, sample.perc = per, quantFun = quantFun, 
              call = scall, ylab = ylab, xlab = xlab, namesFun = namesFun, 
              noObs = length(y))
  class(out) <- "quantSheets"
  return(invisible(out))
}

########################################
### 'Windows' Long-Range Forecasting ###
########################################

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
start.range = 200 #trainLen of first window forecasts
testLen = 1
future = 1
forward = 10
locations = 10
iterations = 100
rawData = sim.dat

n.w = 5 #user specified number of 'windows'
WindowForcs = array(NaN, dim = c(locations, iterations, 1))


#### code below can be computationally expensive ####
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
  ensemb.pred = array(NaN, dim = c(locations, iterations, forward))
  
  
  #Begin future forecasts
  for(f in 1:forward)
  {
    #Generating input data
    input.dat = gen.input.data(trainLen = trainLen,
                               m = m,
                               tau = tau,
                               yTrain = sets$yTrain,
                               rawData = rawData,
                               locations = locations,
                               xTestIndex = sets$xTestIndex,
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
    ensemb.pred[,,f] = testing$predictions
    
    
    
    #Append mean of ensembles to training data
    new.train = rbind(new.train, testing$forecastmean)
    newRaw = rbind(newRaw, testing$forecastmean)
    
    #Update all values for next future point
    trainLen = trainLen + 1
    testindex = testindex + 1
    
  }
  
  WindowForcs = abind(WindowForcs, ensemb.pred, along = 3)
  #print(w)
}
WindowForcs = WindowForcs[,,-1]

################################################################
### Calibrated Uncertainty via Penalized Quantile Regression ###
################################################################

#Declare windows, locations, horizon
locations = 10
n.w = 5
horizon = 10
rawData = sim.dat
tau = 1
horizon = 10
begin = 360

#Declafre penalty vector, wanted quantiles, PI interval
temp_quant_dist = list()
window_num = 1
xsi.seq = seq(0.1, 4, 0.1)
range = 1:horizon
interval = 0.95
wanted_quants = c(2.5, 50, 97.5)


#Begin calibration for each window independently
for(window in 1:n.w)
{
  start.range = begin + (window-1)*horizon
  true.range = (start.range + tau + 1):(start.range + horizon + tau)
  
  q_dist= array(NaN, dim = c(length(range), 2, locations))
  
  #Optimize PI Bands for each location independently
  for(l in 1:locations)
  {
    location = l
    dat = WindowForcs[location,,((window-1)*horizon + range)]
    true.y = rawData[true.range, location]
    means = apply(dat, 2, mean)
    
    interations = dim(dat)[1]
    #Generate residual windows for j-step ahead forecasts
    resid.window = list()
    for(i in 1:horizon)
    {
      resid.window[[i]] = as.numeric(matrix(true.y[i], nrow = iterations, ncol = n.w, byrow = TRUE) - dat[,i])
    }
    
    #Calculate average residual from ensemble 
    temp = unlist(lapply(resid.window, mean))
    temp_df = data.frame('timeframe' = 1:length(range), 'resid' = temp[range])
    
    
    #Optimize xsi parameter in quantile regression
    qcross = rep(NaN, length(xsi.seq))
    for(iter in 1:length(xsi.seq))
    {
      q_mod = qs_mod(resid,
                     timeframe,
                     data = temp_df,
                     cent = wanted_quants,
                     x.lambda = 1, p.lambda = 1, xsi = xsi.seq[iter],
                     totpoints = length(range),
                     print = FALSE,
                     control = qs_cont(degree = 3,
                                       logit = F,
                                       kappa = 0,
                                       plot = FALSE,
                                       c.crit = 1e-5, n.cyc = 50))
      
      qcross[iter] = (sum(as.numeric(t(apply(q_mod$fitted.values, 1, diff))) <= 0) >= 1)
    }
    
    index = which(qcross == 1)[1]-1
    if(!is.na(index) & index !=0)
    {
      optim.xsi = xsi.seq[index]
    }else if(is.na(index)){
      optim.xsi = tail(xsi.seq, 1)
    }else if(index == 0){
      optim.xsi = xsi.seq[1]
    }
    
    
    #Use optim xsi to produce PI bands from penalized quantile regression
    q_mod = qs_mod(resid,
                   timeframe,
                   data = temp_df,
                   cent = wanted_quants,
                   x.lambda = 1, p.lambda = 1, xsi = optim.xsi, 
                   totpoints = length(range),
                   print = FALSE,
                   control = qs_cont(degree = 3,
                                     logit = F,
                                     kappa = 0,
                                     plot = F, c.crit = 1e-5, n.cyc = 50))
    
    
    q_dist[,1,l] = q_mod$fitted.values[,2] - q_mod$fitted.values[,1]
    q_dist[,2,l] = q_mod$fitted.values[,3] - q_mod$fitted.values[,2]
  }
  
  
  
  #Kappa expansion method
  index = window
  tau = 1
  trainLen = begin + (index-1)*horizon
  testLen = 1
  forward = horizon
  
  #Generate training/testing/valid sets
  sets = cttv(rawData, tau, trainLen, testLen = forward)
  
  ensemb.pred = WindowForcs[,,((index-1)*horizon + range)]
  mean.pred = sapply(1:locations, function(x) apply(ensemb.pred[x,,], 2, mean))

  #Optimize kappa to fit wanted interval on observed data
  kappa.grid = seq(0,0.5,0.01)
  o.kappa = rep(0,locations)
  for(l in 1:locations)
  {
    cover.vec = rep(0, length(kappa.grid))
    for(k in 1:length(kappa.grid))
    {
      upper = mean.pred[,l] + (q_dist[,2,l]) + kappa.grid[k]
      lower = mean.pred[,l] - (q_dist[,1,l]) - kappa.grid[k]
      
      cover.vec[k] = mean(sets$yTest[,l] <= upper & sets$yTest[,l] >= lower)
    }
    
    o.kappa[l] = kappa.grid[which.min(abs(cover.vec-interval))]
  }
  
  
  #add kappa expansion to quantile regression bands
  q_dist[,2,] = q_dist[,2,] + matrix(o.kappa, nrow = horizon, ncol = locations, byrow = T)
  q_dist[,1,] = q_dist[,1,] + matrix(o.kappa, nrow = horizon, ncol = locations, byrow = T)
  
  
  
  print(window)
  temp_quant_dist[[window_num]] = q_dist
  window_num = window_num + 1
  
}


#Take average across PI bands
temp_upper = array(NaN, dim = c(horizon,locations,1))
temp_lower = array(NaN, dim = c(horizon,locations,1))
for(i in 1:length(temp_quant_dist))
{
  phold = temp_quant_dist[[i]]
  temp_upper = abind(temp_upper, phold[,2,], along = 3)
  temp_lower = abind(temp_lower, phold[,1,], along = 3)
}
temp_upper = temp_upper[,,-1]
temp_lower = temp_lower[,,-1]

quant_dist = array(NaN, dim = c(horizon,2,locations))
for(l in 1:locations)
{
  quant_dist[,2,l] = apply(temp_upper[,l,], 1, mean)
  quant_dist[,1,l] = apply(temp_lower[,l,], 1, mean)
}

#Save calibrated PI distances
save(quant_dist, file = 'Calibrated_QDist.RData')
