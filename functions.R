#################
#################
### Functions ###
#################
#################

#load libraries
library(tidyverse)
library(Matrix)

set.seed(NULL)
##########################
### Echo State Network ###
##########################

ensemble.esn = function(y.train,
                        x.insamp,
                        x.outsamp,
                        y.test = NULL,
                        n.h,
                        nu,
                        pi.w,
                        pi.win,
                        eta.win,
                        eta.w,
                        lambda.r,
                        alpha,
                        m,
                        iter,
                        future,
                        startvalues = NULL,
                        activation = 'tanh',
                        distribution = 'Normal',
                        scale.factor,
                        scale.matrix,
                        polynomial = 1,
                        verbose = TRUE,
                        parallel = FALSE,
                        fork = FALSE)
{  
  if(!parallel)
  {
    if(verbose)
    {
      prog.bar = txtProgressBar(min = 0, max = iter, style = 3)
    }
  }
  
  #Set conditions for iterations
  cap.t = dim(y.train)[1]
  locations = dim(y.train)[2]
  if(is.null(cap.t) | is.null(locations))
  {
    cap.t = length(y.train)
    locations = 1
  }
  
  n.x = (locations * (m+1)) + 1
  samp.w = n.h*n.h
  samp.win = n.h*n.x
  ensemb.mat = array(0, dim = c(locations, future, iter))
  
  #Set starting values if not specified
  if(is.null(startvalues))
  {
    startvalues = rep(0, n.h)
  }
  
  #Set the activation function
  if(activation == 'identity')
  {
    g.h = function(omega)
    {
      return(omega)
    } 
  } else if(activation == 'tanh') {
    g.h = function(omega)
    {
      placeholder = tanh(omega)
      return(placeholder)
    } 
  }
  
  ###########################
  #Begin ensemble iterations
  ###########################
  
  if(!parallel)
  {
    
    #Non-parallel iterations
    for(k in 1:iter)
    {
      
      #####################################
      #Simulating W and U weight matrices
      gam.w = rbernoulli(samp.w, p = pi.w) 
      gam.win = rbernoulli(samp.win, p = pi.win) 
      
      #Set W
      if(distribution == 'Unif')
      {
        unif.w = runif(samp.w, min = -eta.w, max = eta.w)
        W = Matrix((gam.w == 1)*unif.w + (gam.w == 0)*0, nrow = n.h, ncol = n.h, sparse = T)
      } else if(distribution == 'Normal')
      {
        norm.w = rnorm(samp.w, 0, 1)
        W = Matrix((gam.w == 1)*norm.w + (gam.w == 0)*0, nrow = n.h, ncol = n.h, sparse = T)
      }
      
      #Set W^in
      if(distribution == 'Unif')
      {
        unif.win = runif(samp.win, min = -eta.win, max = eta.win)
        WIN = Matrix((gam.win == 1)*unif.win + (gam.win == 0)*0, nrow = n.h, ncol = n.x, sparse = T)
      } else if(distribution == 'Normal')
      {
        norm.win = rnorm(samp.win, 0, 1)
        WIN = Matrix((gam.win == 1)*norm.win + (gam.win == 0)*0, nrow = n.h, ncol = n.x, sparse = T)
      }
      
      lambda.w = max(abs(eigen(W)$values))
      
      #Set problem to quadratic 
      if(polynomial > 1)
      {
        quad.ridge.dim = polynomial*n.h
        h.prior = rep(0, quad.ridge.dim)
        Ident.Mat = diag(polynomial*n.h) 
        reservoir = matrix(NaN, nrow = quad.ridge.dim, ncol = cap.t)
        h.forc.prior = rep(0, quad.ridge.dim)
        forc.reservoir = matrix(NaN, nrow = quad.ridge.dim, ncol = future)
      } else {
        h.prior = rep(0, n.h)
        Ident.Mat = diag(n.h)
        reservoir = matrix(NaN, nrow = n.h, ncol = cap.t)
        h.forc.prior = rep(0, n.h)
        forc.reservoir = matrix(NaN, nrow = n.h, ncol = future)
      }
      h.prior[1:n.h] = startvalues
      
      #Generate hidden units for training
      t = 1
      while(t <=cap.t)
      {
        WIN.x.in.product  = tcrossprod(WIN, x.insamp)
        omega = (nu/abs(lambda.w)) * W %*% h.prior[1:n.h] + WIN.x.in.product[,t] 
        h.tild.t = g.h(omega)
        h.temp = (1-alpha) * h.prior[1:n.h] + alpha * h.tild.t
        h.prior[1:n.h] = h.temp
        if(polynomial == 2)
        {
          quad.ridge.dim = 2*n.h
          h.prior[(n.h+1):quad.ridge.dim] = h.temp^2
        } else if(polynomial == 3)
        {
          h.prior[(n.h+1):(2*n.h)] = h.temp^2
          h.prior[(2*n.h+1):(3*n.h)] = h.temp^3
        } else if(polynomial == 4)
        {
          h.prior[(n.h+1):(2*n.h)] = h.temp^2
          h.prior[(2*n.h+1):(3*n.h)] = h.temp^3
          h.prior[(3*n.h+1):(4*n.h)] = h.temp^4
        }
        reservoir[,t] = h.prior
        t = t+1
      }
      
      #Estimating V weights using ridge regression
      ridgeMat = lambda.r * Ident.Mat
      V = t(y.train) %*% t(reservoir) %*% solve(tcrossprod(reservoir, reservoir) + ridgeMat)
      
      #Generating out-of-sample hidden units
      f = 1
      h.forc.prior[1:n.h] = reservoir[1:n.h,cap.t]
      while(f <= future)
      {
        WIN.x.out.product  = tcrossprod(WIN, x.outsamp)
        omega.hat = (nu/abs(lambda.w)) * W %*% h.forc.prior[1:n.h] + WIN.x.out.product[,f]
        h.tild.hat = g.h(omega.hat)
        h.hat.temp = (1-alpha) * h.forc.prior[1:n.h] + alpha * h.tild.hat
        h.forc.prior[1:n.h] = h.hat.temp
        if(polynomial == 2)
        {
          quad.ridge.dim = 2*n.h
          h.forc.prior[(n.h+1):quad.ridge.dim] = h.hat.temp^2
        } else if(polynomial == 3)
        {
          h.forc.prior[(n.h+1):(2*n.h)] = h.hat.temp^2
          h.forc.prior[(2*n.h+1):(3*n.h)] = h.hat.temp^3
        } else if(polynomial == 4)
        {
          h.forc.prior[(n.h+1):(2*n.h)] = h.hat.temp^2
          h.forc.prior[(2*n.h+1):(3*n.h)] = h.hat.temp^3
          h.forc.prior[(3*n.h+1):(4*n.h)] = h.hat.temp^4
        }
        forc.reservoir[,f] = h.forc.prior
        f = f+1
      }
      
      #Produce forecasts using coefficient matrix V and forc.reservoir
      ensemb.mat[,,k] = (scale.factor * (V %*% forc.reservoir)) + scale.matrix
      
      if(verbose)
      {
        setTxtProgressBar(prog.bar, k)
      }
    }
    
    
  } else if(parallel) {
    
    
    #Parallel Iterations
    set.seed(NULL)
    #Specify number of clusters
    if(fork)
    {
      cl <- parallel::makeForkCluster(getOption('cores'))
    } else if(!fork) {
      cl <- parallel::makeCluster(getOption('cores'))
    }
    
    # Activate cluster for foreach library
    doParallel::registerDoParallel(cl)
    
    
    ensemb.mat = foreach::foreach(k = 1:iter,
                                  .combine = abind,
                                  .inorder = FALSE) %dopar%
      {
        set.seed(NULL)
        #####################################
        #Simulating W and U weight matrices
        gam.w = purrr::rbernoulli(samp.w, p = pi.w) 
        gam.win = purrr::rbernoulli(samp.win, p = pi.win) 
        
        #Set W
        if(distribution == 'Unif')
        {
          unif.w = runif(samp.w, min = -eta.w, max = eta.w)
          W = Matrix::Matrix((gam.w == 1)*unif.w + (gam.w == 0)*0, nrow = n.h, ncol = n.h, sparse = T)
        } else if(distribution == 'Normal')
        {
          norm.w = rnorm(samp.w, 0, 1)
          W = Matrix::Matrix((gam.w == 1)*norm.w + (gam.w == 0)*0, nrow = n.h, ncol = n.h, sparse = T)
        }
        
        #Set W^in
        if(distribution == 'Unif')
        {
          unif.win = runif(samp.win, min = -eta.win, max = eta.win)
          WIN = Matrix::Matrix((gam.win == 1)*unif.win + (gam.win == 0)*0, nrow = n.h, ncol = n.x, sparse = T)
        } else if(distribution == 'Normal')
        {
          norm.win = rnorm(samp.win, 0, 1)
          WIN = Matrix::Matrix((gam.win == 1)*norm.win + (gam.win == 0)*0, nrow = n.h, ncol = n.x, sparse = T)
        }
        
        lambda.w = max(abs(eigen(W)$values))
        
        #Set problem to quadratic 
        if(polynomial > 1)
        {
          quad.ridge.dim = polynomial*n.h
          h.prior = rep(0, quad.ridge.dim)
          Ident.Mat = diag(polynomial*n.h) 
          reservoir = matrix(NaN, nrow = quad.ridge.dim, ncol = cap.t)
          h.forc.prior = rep(0, quad.ridge.dim)
          forc.reservoir = matrix(NaN, nrow = quad.ridge.dim, ncol = future)
        } else {
          h.prior = rep(0, n.h)
          Ident.Mat = diag(n.h)
          reservoir = matrix(NaN, nrow = n.h, ncol = cap.t)
          h.forc.prior = rep(0, n.h)
          forc.reservoir = matrix(NaN, nrow = n.h, ncol = future)
        }
        h.prior[1:n.h] = startvalues
        
        #Generate hidden units for training
        t = 1
        while(t <=cap.t)
        {
          WIN.x.in.product  = Matrix::tcrossprod(WIN, x.insamp)
          omega = (nu/abs(lambda.w)) * W %*% h.prior[1:n.h] + WIN.x.in.product[,t] 
          h.tild.t = g.h(omega)
          h.temp = (1-alpha) * h.prior[1:n.h] + alpha * h.tild.t
          h.prior[1:n.h] = h.temp
          if(polynomial == 2)
          {
            quad.ridge.dim = 2*n.h
            h.prior[(n.h+1):quad.ridge.dim] = h.temp^2
          } else if(polynomial == 3)
          {
            h.prior[(n.h+1):(2*n.h)] = h.temp^2
            h.prior[(2*n.h+1):(3*n.h)] = h.temp^3
          } else if(polynomial == 4)
          {
            h.prior[(n.h+1):(2*n.h)] = h.temp^2
            h.prior[(2*n.h+1):(3*n.h)] = h.temp^3
            h.prior[(3*n.h+1):(4*n.h)] = h.temp^4
          }
          reservoir[,t] = h.prior
          t = t+1
        }
        
        #Estimating V weights using ridge regression
        ridgeMat = lambda.r * Ident.Mat
        V = t(y.train) %*% t(reservoir) %*% solve(Matrix::tcrossprod(reservoir, reservoir) + ridgeMat)
        
        #Generating out-of-sample hidden units
        f = 1
        h.forc.prior[1:n.h] = reservoir[1:n.h,cap.t]
        while(f <= future)
        {
          WIN.x.out.product  = Matrix::tcrossprod(WIN, x.outsamp)
          omega.hat = (nu/abs(lambda.w)) * W %*% h.forc.prior[1:n.h] + WIN.x.out.product[,f]
          h.tild.hat = g.h(omega.hat)
          h.hat.temp = (1-alpha) * h.forc.prior[1:n.h] + alpha * h.tild.hat
          h.forc.prior[1:n.h] = h.hat.temp
          if(polynomial == 2)
          {
            quad.ridge.dim = 2*n.h
            h.forc.prior[(n.h+1):quad.ridge.dim] = h.hat.temp^2
          } else if(polynomial == 3)
          {
            h.forc.prior[(n.h+1):(2*n.h)] = h.hat.temp^2
            h.forc.prior[(2*n.h+1):(3*n.h)] = h.hat.temp^3
          } else if(polynomial == 4)
          {
            h.forc.prior[(n.h+1):(2*n.h)] = h.hat.temp^2
            h.forc.prior[(2*n.h+1):(3*n.h)] = h.hat.temp^3
            h.forc.prior[(3*n.h+1):(4*n.h)] = h.hat.temp^4
          }
          forc.reservoir[,f] = h.forc.prior
          f = f+1
        }
        
        #Produce forecasts using coefficient matrix V and forc.reservoir
        (scale.factor * (V %*% forc.reservoir)) + scale.matrix
      }
    
  }
  
  #Close progress bar or clusters
  if(!parallel)
  {
    if(verbose)
    {
      close(prog.bar)
    }
  } else if(parallel) {
    parallel::stopCluster(cl)
  }
  
  #Calculate forecast mean
  if(!parallel)
  {
    if(future > 1)
    {
      forc.mean = sapply(1:locations, function(n) rowMeans(ensemb.mat[n,,]))
    } else if(locations > 1) {
      forc.mean = apply(ensemb.mat[,1,], 1, mean)
    } else {
      forc.mean = mean(ensemb.mat[1,1,])
    }
  } else if(parallel) {
    if(locations > 1 & future == 1)
    {
      forc.mean = apply(ensemb.mat, 1, mean)
    } else if(locations == 1 & future > 1){
      forc.mean = (sapply(1:future, function(x) mean(ensemb.mat[,seq(x, ncol(ensemb.mat), future)])))
    } else if(locations > 1 & future > 1) {
      forc.mean = t(sapply(1:future, function(x) rowMeans(ensemb.mat[,seq(x, ncol(ensemb.mat), future)])))
    } else if(locations == 1 & future == 1) {
      forc.mean = mean(as.numeric(ensemb.mat))
    } else {
      forc.mean = NULL
    }
  }
  
  #Calculate MSE
  if(!is.null(y.test))
  {
    MSE=sum((y.test-forc.mean)^2)/(locations*future)
  } else {
    MSE = NULL
  }
  
  #Set parallel output items
  if(parallel)
  {
    V = NULL
    reservoir = NULL
    forc.reservoir = NULL
  }
  
  #Compile results
  esn.output = list('predictions' = ensemb.mat,
                    'forecastmean' = forc.mean,
                    'y.test' = y.test,
                    'MSE' = MSE,
                    'coefficients' = V,
                    'finalreservoir' = reservoir,
                    'finalforecastreservoir' = forc.reservoir)
  return(esn.output)
}




###########################
### Generate Input Data ###
###########################


gen.input.data = function(trainLen,
                          m,
                          tau,
                          yTrain,
                          rawData,
                          locations,
                          xTestIndex,
                          testLen)
{
  in.sample.len = trainLen - (m * tau)
  
  in.sample.x.raw = array(NA, dim = c(in.sample.len, m+1, locations))
  
  for(i in 1:in.sample.len)
  {
    in.sample.x.raw[i,,] = rawData[seq(i, (m*tau + i), by=tau), ]
  }
  
  #Scale in-sample x and y
  in.sample.y.raw = yTrain[(m*tau + 1):trainLen,]
  y.mean = mean(in.sample.y.raw)
  y.scale = sd(in.sample.y.raw)
  
  in.sample.y = (in.sample.y.raw - y.mean)/y.scale
  
  
  mean.train.x = mean(rawData[1:trainLen,])
  sd.train.x = sd(rawData[1:trainLen,])
  
  
  in.sample.x=(in.sample.x.raw - mean.train.x)/sd.train.x
  
  
  designMatrix = matrix(1,in.sample.len, (m + 1)*locations + 1)
  for(i in 1:in.sample.len){
    designMatrix[i,2:((m + 1)*locations + 1)] = as.vector(in.sample.x[i,,])
  }
  
  
  #Out-Sample
  out.sample.x.raw = array(NA, dim = c(testLen, m + 1, locations))
  for(i in 1:testLen)
  {
    out.sample.x.raw[i,,] = rawData[seq(xTestIndex[i]-(m*tau), xTestIndex[i], by=tau),]
  }
  
  
  #Scale out-sample x and y
  out.sample.x = (out.sample.x.raw - mean.train.x)/sd.train.x
  
  designMatrixOutSample = matrix(1, testLen, (m + 1)*locations + 1)
  for(i in 1:testLen)
  {
    designMatrixOutSample[i,2:((m + 1)*locations + 1)] = as.vector(out.sample.x[i,,])
  }
  
  
  
  #Additive scale matric
  addScaleMat = matrix(y.mean, locations, testLen)
  
  input.data.output = list('y.mean' = y.mean,
                           'y.scale' = y.scale,
                           'in.sample.y' = in.sample.y,
                           'in.sample.x' = in.sample.x,
                           'out.sample.x' = out.sample.x,
                           'in.sample.len' = in.sample.len,
                           'designMatrix' = designMatrix,
                           'designMatrixOutSample' = designMatrixOutSample,
                           'testLen' = testLen,
                           'addScaleMat' = addScaleMat)
  return(input.data.output)
}


#######################################################
### Generate Training, Testing, and Validation Sets ###
#######################################################


cttv = function(rawData, tau, trainLen, testLen, validLen = NULL, valid.flag = FALSE)
{
  #Create training and testing sets
  totlength = trainLen + testLen + tau
  yTrain = rawData[(tau+1):(trainLen+tau),]
  yTest = rawData[(trainLen+tau+1):totlength,]
  xTestIndex = seq((trainLen+1), (totlength-tau), 1)
  
  #Create valid sets
  if(valid.flag)
  {
    xValTestIndex=(trainLen+1-validLen):(trainLen)
    yValTestIndex=(trainLen+tau+1-validLen):(trainLen+tau)
    yValid = rawData[yValTestIndex,]
  } else {
    yValid = NULL
    xValTestIndex = NULL
  }
  
  #Return list
  output = list('yTrain' = yTrain,
                'yTest' = yTest,
                'yValid' = yValid,
                'xTestIndex' = xTestIndex,
                'xValTestIndex' = xValTestIndex)
  return(output)
}



