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




###############################
### Deep Echo-State Network ###
###############################


deep.esn = function(y.train,
                    x.insamp,
                    x.outsamp,
                    y.test = NULL,
                    n.h,
                    nu,
                    pi.w, 
                    pi.win,
                    eta.w,
                    eta.win,
                    lambda.r,
                    alpha,
                    m,
                    iter,
                    future,
                    layers = 3,
                    reduced.units,
                    startvalues = NULL,
                    activation = 'tanh',
                    distribution = 'Normal',
                    logNorm = FALSE,
                    scale.factor,
                    scale.matrix,
                    parallel = F,
                    fork = F,
                    verbose = T)
{
  if(!parallel)
  {
    if(verbose)
    {
      prog.bar = txtProgressBar(min = 0, max = iter, style = 3)
    }
  }
  
  
  ###########################################
  ### Initial Conditions and Known Values ###
  ###########################################
  
  #Set training length and locations
  cap.t = dim(y.train)[1]
  locations = dim(y.train)[2]
  if(is.null(locations) | is.null(cap.t))
  {
    cap.t = length(y.train)
    locations = 1
  }
  
  #Get number of samples for weight matrices
  samp.w = list()
  samp.win = list()
  for(ell in 1:layers)
  {
    samp.w[[ell]] = n.h[ell] * n.h[ell]
    if(ell == 1)
    {
      n.x = (locations * (m+1)) + 1
      samp.win[[1]] = n.h[ell] * n.x
    } else {
      samp.win[[ell]] = n.h[ell] * (reduced.units+1)
    }
  }
  
  #Starting values of hidden units
  if(is.null(startvalues))
  {
    startvalues = list()
    for(ell in 1:layers)
    {
      startvalues[[ell]] = rep(0, n.h[ell])
    }
  }
  
  #Set the activation function
  if(activation == 'identity')
  {
    g.h = function(x)
    {
      return(x)
    } 
  } else if(activation == 'tanh') {
    g.h = function(x)
    {
      placeholder = tanh(x)
      return(placeholder)
    } 
  }
  
  #Set output array
  if(!parallel)
  {
    ensemb.mat = array(0, dim = c(locations, testLen, iter))
  }
  
  #########################
  ### Forecast Ensemble ###
  #########################
  set.seed(NULL)
  
  if(parallel)
  {
    
    #Specify Parallel clusters
    if(fork)
    {
      cl = parallel::makeForkCluster(getOption('cores')) 
    } else if(!fork)
    {
      cl = parallel::makeCluster(getOption('cores'))
    }
    
    #Activate clusters
    doParallel::registerDoParallel(cl)
    
    #Begin parallel iterations
    ensemb.mat = foreach::foreach(k = 1:iter,
                                  .combine = abind,
                                  .inorder = FALSE) %dopar%
      {
        set.seed(NULL)
        
        ##########################################
        ### Generate W and WIN weight matrices ###
        ##########################################
        W = list()
        WIN = list()
        lambda.w = c()
        for(ell in 1:layers)
        {
          #Set sparsity
          gam.w = purrr::rbernoulli(samp.w[[ell]], p = pi.w[ell])
          gam.win = purrr::rbernoulli(samp.win[[ell]], p = pi.win[ell])
          
          #Generate W
          if(distribution == 'Unif')
          {
            unif.w = runif(samp.w[[ell]], min = -eta.w[ell], max = eta.w[ell])
            W[[ell]] = Matrix::Matrix((gam.w == 1)*unif.w + (gam.w == 0)*0,
                                      nrow = n.h[ell], ncol = n.h[ell], sparse = T)
          } else if(distribution == 'Normal')
          {
            norm.w = rnorm(samp.w[[ell]], 0, 1)
            W[[ell]] = Matrix::Matrix((gam.w == 1)*norm.w + (gam.w == 0)*0,
                                      nrow = n.h[ell], ncol = n.h[ell], sparse = T)
          }
          
          #Generate W^in
          n.input = c(n.x, rep(reduced.units+1, (layers-1)))
          if(distribution == 'Unif')
          {
            unif.win = runif(samp.win[[ell]], min = -eta.win[ell], max = eta.win[ell])
            WIN[[ell]] = Matrix::Matrix((gam.win == 1)*unif.win + (gam.win == 0)*0,
                                        nrow = n.h[ell], ncol = n.input[ell], sparse = T)
          } else if(distribution == 'Normal')
          {
            norm.win = rnorm(samp.win[[ell]], 0, 1)
            WIN[[ell]] = Matrix::Matrix((gam.win == 1)*norm.win + (gam.win == 0)*0,
                                        nrow = n.h[ell], ncol = n.input[ell], sparse = T)
          }
          
          #Specify spectral radius
          lambda.w[ell] = max(abs(eigen(W[[ell]])$values))
          
        }
        
        
        ###############################
        ### Initialize Hidden Units ###
        ###############################
        h.prior = list()
        reservoir = list()
        h.forc.prior = list()
        forc.reservoir = list()
        Ident.Mat = diag(((layers-1)*reduced.units) + n.h[layers])
        
        for(ell in 1:layers)
        {
          h.prior[[ell]] = startvalues[[ell]]
          reservoir[[ell]] = matrix(NaN, nrow = n.h[ell], ncol = cap.t)
          h.forc.prior[[ell]] = rep(0, n.h[ell])
          forc.reservoir[[ell]] = matrix(NaN, nrow = n.h[ell], ncol = future)
        }
        
        
        ####################################
        ### Update Training Hidden Units ###
        ####################################
        input.data = list()
        input.data[[1]] = x.insamp
        output.data = list()
        output.data[[1]] = x.outsamp
        for(ell in 1:layers)
        {
          WIN.x.in.product = Matrix::tcrossprod(WIN[[ell]], input.data[[ell]])
          for(t in 1:cap.t)
          {
            omega = g.h(as.matrix((nu[ell]/lambda.w[ell]) * W[[ell]] %*% h.prior[[ell]] + WIN.x.in.product[,t]))
            h.prior[[ell]] = (1-alpha)*h.prior[[ell]] + alpha*omega
            reservoir[[ell]][,t] = as.numeric(h.prior[[ell]])
          } 
          
          h.forc.prior[[ell]] = reservoir[[ell]][,cap.t]
          WIN.x.out.product = Matrix::tcrossprod(WIN[[ell]], output.data[[ell]])
          for(fut in 1:future)
          {
            omega.hat = g.h(as.matrix((nu[ell]/lambda.w[ell]) * W[[ell]] %*% h.forc.prior[[ell]] + WIN.x.out.product[,fut]))
            h.forc.prior[[ell]] = (1-alpha)*h.forc.prior[[ell]] + alpha*omega.hat
            forc.reservoir[[ell]][,fut] = as.numeric(h.forc.prior[[ell]])
          } 
          
          
          #Dimension reduction to combine layers
          if(layers > 1)
          {
            placeholder = wql::eof(cbind(reservoir[[ell]], forc.reservoir[[ell]]), n = reduced.units, scale. = FALSE) 
            mean.pca = apply(placeholder$REOF[1:cap.t,], 2, mean)
            sd.pca = apply(placeholder$REOF[1:cap.t,], 2, sd)
            placeholder$REOF = (placeholder$REOF - matrix(mean.pca, nrow = cap.t+future, ncol = ncol(placeholder$REOF), byrow = TRUE)) / 
              matrix(sd.pca, nrow = cap.t+future, ncol = ncol(placeholder$REOF), byrow = TRUE)
            input.data[[ell+1]] = cbind(rep(1,cap.t), placeholder$REOF[1:cap.t,1:reduced.units])
            output.data[[ell+1]] = cbind(as.matrix(rep(1,future), ncol = 1),
                                         matrix(placeholder$REOF[(cap.t+1):(cap.t+future),1:reduced.units], nrow = future, ncol = reduced.units))
          } else {
            input.data[[ell+1]] = NULL
            output.data[[ell+1]] = NULL
          }
          
        } 
        
        
        ###################################
        ### Estimate Coefficient Matrix ###
        ###################################
        
        #Get dimension reduced data on same scale
        h.tild = matrix(NaN, nrow = cap.t, ncol = 1)
        for(ell in 2:layers)
        {
          h.tild = cbind(h.tild, g.h(input.data[[ell]][,-1]))
        }
        h.tild = h.tild[,-1]
        
        #Estimate coefficients
        final.design = rbind(reservoir[[layers]], t(h.tild))
        ridgeMat = lambda.r * Ident.Mat
        V = t(y.train) %*% t(final.design) %*% solve(Matrix::tcrossprod(final.design, final.design) + ridgeMat)
        
        
        ###########################
        ### Calculate Forecasts ###
        ###########################
        
        #Get dimension reduced data on same scale
        h.tild.out = matrix(NaN, nrow = future, ncol = 1)
        for(ell in 2:layers)
        {
          h.tild.out = cbind(h.tild.out, matrix(g.h(output.data[[ell]][,-1]), ncol = reduced.units))
        }
        h.tild.out = matrix(h.tild.out[,-1], ncol = ((layers-1)*reduced.units))
        
        #Create output design matrix
        final.design.out = rbind(forc.reservoir[[layers]], t(h.tild.out))
        
        #Generate forecasts
        if(logNorm)
        {
          exp((scale.factor * (V %*% final.design.out)) + scale.matrix)
        } else {
          (scale.factor * (V %*% final.design.out)) + scale.matrix
        }
        
        
      } 
    
    
  } else {
    
    
    
    #Begin non-parallel iterations
    for(k in 1:iter)
    {
      set.seed(NULL)
      
      ##########################################
      ### Generate W and WIN weight matrices ###
      ##########################################
      W = list()
      WIN = list()
      lambda.w = c()
      for(ell in 1:layers)
      {
        #Set sparsity
        gam.w = purrr::rbernoulli(samp.w[[ell]], p = pi.w[ell])
        gam.win = purrr::rbernoulli(samp.win[[ell]], p = pi.win[ell])
        
        #Generate W
        if(distribution == 'Unif')
        {
          unif.w = runif(samp.w[[ell]], min = -eta.w[ell], max = eta.w[ell])
          W[[ell]] = Matrix::Matrix((gam.w == 1)*unif.w + (gam.w == 0)*0,
                                    nrow = n.h[ell], ncol = n.h[ell], sparse = T)
        } else if(distribution == 'Normal')
        {
          norm.w = rnorm(samp.w[[ell]], 0, 1)
          W[[ell]] = Matrix::Matrix((gam.w == 1)*norm.w + (gam.w == 0)*0,
                                    nrow = n.h[ell], ncol = n.h[ell], sparse = T)
        }
        
        #Generate W^in
        n.input = c(n.x, rep(reduced.units+1, (layers-1)))
        if(distribution == 'Unif')
        {
          unif.win = runif(samp.win[[ell]], min = -eta.win[ell], max = eta.win[ell])
          WIN[[ell]] = Matrix::Matrix((gam.win == 1)*unif.win + (gam.win == 0)*0,
                                      nrow = n.h[ell], ncol = n.input[ell], sparse = T)
        } else if(distribution == 'Normal')
        {
          norm.win = rnorm(samp.win[[ell]], 0, 1)
          WIN[[ell]] = Matrix::Matrix((gam.win == 1)*norm.win + (gam.win == 0)*0,
                                      nrow = n.h[ell], ncol = n.input[ell], sparse = T)
        }
        
        #Specify spectral radius
        lambda.w[ell] = max(abs(eigen(W[[ell]])$values))
        
      }
      
      
      ###############################
      ### Initialize Hidden Units ###
      ###############################
      h.prior = list()
      reservoir = list()
      h.forc.prior = list()
      forc.reservoir = list()
      Ident.Mat = diag(((layers-1)*reduced.units) + n.h[layers])
      
      for(ell in 1:layers)
      {
        h.prior[[ell]] = startvalues[[ell]]
        reservoir[[ell]] = matrix(NaN, nrow = n.h[ell], ncol = cap.t)
        h.forc.prior[[ell]] = rep(0, n.h[ell])
        forc.reservoir[[ell]] = matrix(NaN, nrow = n.h[ell], ncol = future)
      }
      
      
      ####################################
      ### Update Training Hidden Units ###
      ####################################
      input.data = list()
      input.data[[1]] = x.insamp
      output.data = list()
      output.data[[1]] = x.outsamp
      for(ell in 1:layers)
      {
        WIN.x.in.product = Matrix::tcrossprod(WIN[[ell]], input.data[[ell]])
        for(t in 1:cap.t)
        {
          omega = g.h(as.matrix((nu[ell]/lambda.w[ell]) * W[[ell]] %*% h.prior[[ell]] + WIN.x.in.product[,t]))
          h.prior[[ell]] = (1-alpha)*h.prior[[ell]] + alpha*omega
          reservoir[[ell]][,t] = as.numeric(h.prior[[ell]])
        } 
        
        h.forc.prior[[ell]] = reservoir[[ell]][,cap.t]
        WIN.x.out.product = Matrix::tcrossprod(WIN[[ell]], output.data[[ell]])
        for(fut in 1:future)
        {
          omega.hat = g.h(as.matrix((nu[ell]/lambda.w[ell]) * W[[ell]] %*% h.forc.prior[[ell]] + WIN.x.out.product[,fut]))
          h.forc.prior[[ell]] = (1-alpha)*h.forc.prior[[ell]] + alpha*omega.hat
          forc.reservoir[[ell]][,fut] = as.numeric(h.forc.prior[[ell]])
        } 
        
        
        #Dimension reduction to combine layers
        if(layers > 1)
        {
          placeholder = wql::eof(cbind(reservoir[[ell]], forc.reservoir[[ell]]), n = reduced.units, scale. = FALSE) 
          mean.pca = apply(placeholder$REOF[1:cap.t,], 2, mean)
          sd.pca = apply(placeholder$REOF[1:cap.t,], 2, sd)
          placeholder$REOF = (placeholder$REOF - matrix(mean.pca, nrow = cap.t+future, ncol = ncol(placeholder$REOF), byrow = TRUE)) / 
            matrix(sd.pca, nrow = cap.t+future, ncol = ncol(placeholder$REOF), byrow = TRUE)
          input.data[[ell+1]] = cbind(rep(1, cap.t), placeholder$REOF[1:cap.t,1:reduced.units])
          output.data[[ell+1]] = cbind(as.matrix(rep(1,future), ncol = 1),
                                       matrix(placeholder$REOF[(cap.t+1):(cap.t+future),1:reduced.units], nrow = future, ncol = reduced.units))
        } else {
          input.data[[ell+1]] = NULL
          output.data[[ell+1]] = NULL
        }
        
      } 
      
      
      ###################################
      ### Estimate Coefficient Matrix ###
      ###################################
      
      #Get dimension reduced data on same scale
      h.tild = matrix(NaN, nrow = cap.t, ncol = 1)
      for(ell in 2:layers)
      {
        h.tild = cbind(h.tild, g.h(input.data[[ell]][,-1]))
      }
      h.tild = h.tild[,-1]
      
      #Estimate coefficients
      final.design = rbind(reservoir[[layers]], t(h.tild))
      ridgeMat = lambda.r * Ident.Mat
      V = t(y.train) %*% t(final.design) %*% solve(Matrix::tcrossprod(final.design, final.design) + ridgeMat)
      
      
      
      ###########################
      ### Calculate Forecasts ###
      ###########################
      
      #Get dimension reduced data on same scale
      h.tild.out = matrix(NaN, nrow = future, ncol = 1)
      for(ell in 2:layers)
      {
        h.tild.out = cbind(h.tild.out, matrix(g.h(output.data[[ell]][,-1]), ncol = reduced.units))
      }
      h.tild.out = matrix(h.tild.out[,-1], ncol = ((layers-1)*reduced.units))
      
      #Create output design matrix
      final.design.out = rbind(forc.reservoir[[layers]], t(h.tild.out))
      
      #Generate forecasts
      if(logNorm)
      {
        ensemb.mat[,,k] = exp((scale.factor * (V %*% final.design.out)) + scale.matrix)
      } else {
        ensemb.mat[,,k] = (scale.factor * (V %*% final.design.out)) + scale.matrix
      }
      
      #update progress bare
      if(verbose)
      {
        setTxtProgressBar(prog.bar, k)
      }
      
    } 
    
  }
  
  
  ########################
  ### Finalize Results ###
  ########################
  
  #Close parallel clusters
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
    if(testLen > 1)
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
  
  #Compile results
  esn.output = list('predictions' = ensemb.mat,
                    'forecastmean' = forc.mean,
                    'MSE' = MSE)
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



