

#' @title  Semi-varying coefficient model with jump points (SVCJP)
#'
#' @description  Jump detection and coefficient estimation for the semi-varying coefficient model with jump points.
#'
#' @param tin A matrix containing the observed values of the index variables. Each column corresponds to one index variable.
#' @param yin A vector containing the observed values of the response.
#' @param xin A matrix containing the observed values of the covariates with varying coefficients. Each column corresponds to one covariate.
#' @param zin A matrix containing the observed values of the covariates with linear (non-varying) coefficients. Each column corresponds to one covariate.
#' @param win A vector of weights for each observation. If win = NULL (Default), equal weights will be assigned.
#' @param zeta The threshold value for change point (jump) detection. If zeta = NULL (Default), the threshold function will be calculated based on theoretical results.
#' @param bw.seq1 A vector containing the candidate values for bandwidth selection in Part 1 (for h_1). If bw.seq1 = NULL (Default), the candidate values will be generated based on theoretical results.
#' @param bw.seq2 A vector containing the candidate values for bandwidth selection in Part 2 (for h_tau, h_d, and h_2). If bw.seq2 = NULL (Default), the candidate values will be generated based on theoretical results.
#' @param kernel The smoothing kernel used for estimation, including "epan" (Default, the Epanechnikov kernel), "rect", "gauss", "gausvar", and "quar".
#' @param NbGrid The number of grid points used for jump detection in each index dimension. Default is 101.
#' @param nRegGrid The number of grid points to generate output estimation, mainly for plotting purpose. Default is 101.
#' @param kFolds The number of folds for bandwidth cross validation selection. Default is 5.
#' @param npoly The degree of polynomial. Default is 1 for local linear smoothing.
#' @param nder The order of derivative, which should be smaller than npoly. Default is 0 for local linear smoothing.
#' @param hkappa A numeric value (> 1) for neighborhood width adjustment in jump detection to avoid repeated identification of the same jump. Default is 2.
#' @param Refine A boolean indicating whether to perform the Refining Stage (iterative refining) in jump detection.
#' @param refine_tol The tolerance value for the Refining Stage. Default is 1e-4.
#' @param refine_maxiter The maximum number of iterations allowed in the Refining Stage. Default is 10.
#' @param verbose A boolean indicating whether to output more information during function execution.
#'
#' @return A list containing the following fields: 
#' \item{beta_hat}{A vector of estimated linear coefficients.}
#' \item{alp_est}{A list of estimated varying coefficient functions, including detected change point locations (jumptime), estimated jump sizes (jumpsize), and the estimated coefficient function values on the inputs (alp.hat)}
#' \item{yhat}{A vector of estimated response given the input.}
#' \item{mse}{The mean squared error of response.}
#' \item{tuning_parameters}{The tuning parameters selected for our SVCJP method, including h_1, h_tau, h_d, h_2, zeta, and rho_d.}
#'
#' @examples
#' data(syn_data_2d)
#' res_2d <- SVCJP(tin = tin, yin = yin, xin = xin, zin = zin)
#' 

SVCJP <- function(tin, 
                  yin, 
                  xin, 
                  zin, 
                  win = NULL,
                  zeta = NULL,  
                  bw.seq1 = NULL, 
                  bw.seq2 = NULL,
                  kernel= c("epan", "rect", "gauss", "gausvar", "quar"),
                  NbGrid = 101, 
                  nRegGrid = 101,
                  kFolds = 5,
                  npoly = 1, 
                  nder = 0,  
                  hkappa = 2,
                  Refine = T,
                  refine_tol = 1e-4,
                  refine_maxiter = 10,
                  verbose = T){
  
  kernel = match.arg(kernel)
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  q = dim(zin)[2]
  m = dim(tin)[2]
  
  t_ord = order(tin[,1])
  
  # weight
  if (is.null(win)) {
    win = rep(1, n)
    win = win/sum(win)
  } else {
    if (length(win) != n){
      stop("The length of provided weight win does not match data dimension.")
    }
    if (sum(win>=0) < n){
      stop("Weight win needs to be a vector of nonnegative numbers.")
    }
    win = win[t_ord]
    win = win/sum(win)
  }

  yin = yin[t_ord]
  xin = xin[t_ord, , drop=F]
  zin = zin[t_ord, , drop=F]
  tin = tin[t_ord, , drop=F]
  
  ## Bandwidth selection range
  bw_rn = 6
  bw1_base = n^{-0.175}/15*m
  if (is.null(bw.seq1)) bw.seq1 = seq(from = bw1_base/2, to = bw1_base*2, length.out = bw_rn)
  bw2_base = max(bw1_base*2, n^{-0.225}/2.3*m)
  if (is.null(bw.seq2)) bw.seq2 = seq(from = bw2_base/3*2, to = bw2_base*3/2, length.out = bw_rn)
  
  ##### Part 1: Estimation of the linear part
  
  if (verbose) message('Part 1: Estimation of the linear part.')
  
  ## k-fold cross validation to choose h_1
  
  cv1 = array(Inf, dim = c(length(bw.seq1), kFolds));
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  for (j in 1:length(bw.seq1)) {
    
    for (i in 1:kFolds) {
      
      ttest <- tin[theFolds[[i]], , drop=F]
      ytest <- yin[theFolds[[i]]]
      xtest <- xin[theFolds[[i]], , drop=F]
      ztest <- zin[theFolds[[i]], , drop=F]
      ttrain <- tin[-theFolds[[i]], , drop=F]
      ytrain <- yin[-theFolds[[i]]]
      xtrain <- xin[-theFolds[[i]], , drop=F]
      ztrain <- zin[-theFolds[[i]], , drop=F]
      wtrain <- win[-theFolds[[i]]]
      
      coef_cv1 = tryCatch(
        CPPlwls2d_s1(bw = bw.seq1[j], kernel_type = kernel, win = wtrain,
                     tin = ttrain, yin = ytrain, xin = xtrain, zin = ztrain, tout = ttest, npoly = npoly),
        error=function(err) {
          return(Inf)
        })
      nan_rate <- length(which(is.nan(coef_cv1[,1]))) / length(ytest)
      if (nan_rate > 0.2) warning('Bandwidth h1 =', bw.seq1[j], 'too small for Part 1 CV.')
      tid <- which(!is.nan(coef_cv1[,1]))
      beta_cv1 = colMeans(coef_cv1[tid, -c(1:(p+1)), drop=F])
      alp_cv1 = coef_cv1[tid, 1:(p+1), drop=F]
      
      yhat_cv1 = rowMeans(cbind(1, xtest[tid, , drop=F]) * alp_cv1) + ztest[tid, , drop=F] %*% beta_cv1
      
      cv1[j,i] = sum((ytest[tid] - c(yhat_cv1))^2)
      if (nan_rate >= 0.5) cv1[j,i] = Inf
      # print(cv1)
      if(is.na(cv1[j,i]) || is.nan(cv1[j,i])){
        cv1[j,i] = Inf
      }
      
    }
  }
  
  if(min(cv1) == Inf){
    stop("All bandwidths result in infinite CV costs. (Part 1)")
  }
  
  cvMean = apply(cv1, c(1), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  bopt1 = bw.seq1[cvMeanid];
  names(bopt1) = c('h_1')
  
  ## Use the chosen bandwidth h_1 to estimate beta
  
  h_1= unname(bopt1)
  
  beta_est = CPPlwls2d_s1(bw = h_1, kernel_type = kernel, win = win,
                          tin = tin, yin = yin, xin = xin, zin = zin,
                          tout = tin, npoly = npoly)[,-c(1:(p+1)), drop=F]
  tid <- which(!is.nan(beta_est[,1]))
  beta_hat = colMeans(beta_est[tid, , drop=F])
  
  #
  yin_s2 = c(yin - zin %*% beta_hat)
  
  if (verbose) message('  Bandwidth for Part 1: h_1 = ', h_1)
  
  ##### Part 2: Estimation of the nonparametric part
  
  if (verbose) message('Part 2: Estimation of the nonparametric part.')
  
  if(length(bw.seq2) == 1) {
    h_tau = bw.seq2
    h_d = 2 * bw.seq2
    zeta = zetaFun(bw = bw.seq2, 
                   tin = tin, yin = yin_s2, xin = xin, win = win,
                   npoly = npoly, nder = nder, kernel = kernel,
                   nRegGrid = nRegGrid)
  } else {
    
    ## select by CV procedure
    if (verbose) message('  Selecting tuning parameters (h_tau, h_d, h_2) for Part 2 ...')
    tunings = CVbandwidth(bw.seq2 = bw.seq2, zeta = zeta, win = win,
                          tin = tin, yin = yin_s2, xin = xin,
                          npoly = npoly, nder = nder, kernel = kernel,
                          NbGrid = NbGrid, nRegGrid = nRegGrid,
                          jumpcheck = T, kFolds = kFolds, 
                          hkappa = hkappa)
    h_tau = unname(tunings$bopt['h_tau'])
    h_d = unname(tunings$bopt['h_d'])
    zeta = tunings$zeta
    
    # k-fold cross validation to choose h_2
    h2.seq = seq(from = h_tau+(h_d-h_tau)/20, to = h_d-(h_d-h_tau)/2,
                 length.out = length(bw.seq2))
    cv2 = array(Inf, dim = c(length(h2.seq), kFolds));
    theFolds = SimpleFolds(1:n, kFolds)
    
    for (j in 1:length(h2.seq)) {
      
      for (i in 1:kFolds) {
        
        ttest <- tin[theFolds[[i]], , drop=F]
        ytest <- yin_s2[theFolds[[i]]]
        xtest <- xin[theFolds[[i]], , drop=F]
        ttrain <- tin[-theFolds[[i]], , drop=F]
        ytrain <- yin_s2[-theFolds[[i]]]
        xtrain <- xin[-theFolds[[i]], , drop=F]
        wtrain <- win[-theFolds[[i]]]
        
        muout = tryCatch(
          CoefJump(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain,
                   tout = ttest, xout = xtest,
                   h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h2.seq[j],
                   jumpcheck = T, npoly=npoly, nder= nder,
                   kernel = kernel,  NbGrid = NbGrid,
                   hkappa = hkappa, verbose = F)$muout,
          error=function(err) {
            return(Inf)
          })
        nan_rate <- length(which(is.nan(muout))) / length(ytest)
        if (nan_rate > 0.2) warning('Bandwidth: h_2=', h2.seq[j], 'too small for varying coef CV.')
        tid <- which(!is.nan(muout))
        
        cv2[j,i] = sum((ytest[tid] - muout[tid])^2)
        # print(cv2)
        if(is.na(cv2[j,i]) || is.nan(cv2[j,i]) || nan_rate>0.2){
          cv2[j,i] = Inf;
        }
        
      }
    }
    
    if(min(cv2) == Inf){
      stop("All bandwidths resulted in infinite CV costs. (Part 2, h_2)")
    }
    
    cvMean = apply(cv2, c(1), mean)
    cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
    if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[floor(dim(cvMeanid)[1]/2),]
    bopt2 = h2.seq[cvMeanid];
    names(bopt2) = c('h_2')
    
    h_2 = unname(bopt2)
    
  }
  
  ##
  if (verbose) {
    message('  Bandwidths for Part 2: ', 'h_tau = ', h_tau, '; h_d = ', h_d, '; h_2 = ', h_2)
    message('  Threshold: ', 'zeta = ', zeta)
  }
  
  ###
  obsGrid = tin
  ## cut in the region [h, 1-h]
  regGrid = apply(obsGrid, 2, function(x){
    seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = nRegGrid)
  })
  regGrid = as.matrix(expand.grid(as.data.frame(regGrid)))
  
  # Get the mean function using the bandwidth estimated above:
  smcObj = CoefJump(tin = tin, yin = yin_s2, xin = xin, win = win, 
                    tout = regGrid, xout = NULL,
                    h_tau = h_tau, h_d = h_d, zeta = zeta, h_2 = h_2,
                    jumpcheck = T, npoly = npoly, nder = nder, 
                    kernel = kernel, NbGrid = NbGrid,
                    hkappa = hkappa,
                    Refine = T,
                    refine_tol = refine_tol,
                    refine_maxiter = refine_maxiter,
                    verbose = verbose)
  
  alpha_hat <- sapply(smcObj$alp_est, function(z) z$alp.hat)
  yhat <- c(zin %*% beta_hat + rowSums(alpha_hat * cbind(1, xin)))
  
  mse <- mean((yhat[!is.na(yhat)]-yin[!is.na(yhat)])^2)
  
  if (verbose) message('MSE: ', mse)
  
  # Correspond the outputs to the inputs (reverse the ordering)
  alp_est = smcObj$alp_est
  for (j in 1:(p+1)) {
    alp_est[[j]]$gamma.hat[t_ord] = alp_est[[j]]$gamma.hat
    alp_est[[j]]$alp.hat[t_ord] = alp_est[[j]]$alp.hat
  }
  yhat[t_ord] = yhat
  
  
  ##
  res <- list(beta_hat = beta_hat, 
              alp_est = alp_est, 
              yhat = yhat, 
              mse = mse,
              tunning_parameters = list(h_1 = h_1, 
                                        h_tau = h_tau, 
                                        h_d = h_d, 
                                        h_2 = h_2,
                                        zeta = zeta, 
                                        rho_d = unname(smcObj$rho_d))
  )
  
  return(res)
  
  
}





#' @title Varying coefficient model with jump points 
#'
#' @description Jump detection and coefficient estimation for the varying coefficient model with jump points with given tuning parameters.
#'
#' @param tin A matrix containing the observed values of the index variables. Each column corresponds to one index variable.
#' @param yin A vector containing the observed values of the response.
#' @param xin A matrix containing the observed values of the covariates. Each column corresponds to one covariate.
#' @param win A vector of weights for each observation. If win = NULL (Default), equal weights will be assigned.
#' @param tout A matrix containing the testing values of the index variables. Default is NULL.
#' @param xout A matrix containing the testing values of the covariates. Default is NULL.
#' @param h_tau The bandwidth for jump location estimation.
#' @param h_d The bandwidth for jump size estimation.
#' @param h_2 The bandwidth for local kernel smoothing.
#' @param zeta The threshold value for change point (jump) detection. If zeta = NULL (Default), the threshold function will be calculated based on theoretical results.
#' @param kernel The smoothing kernel used for estimation, including "epan" (Default, the Epanechnikov kernel), "rect", "gauss", "gausvar", and "quar".
#' @param NbGrid The number of grid points used for jump detection in each index dimension. Default is 101.
#' @param npoly The degree of polynomial. Default is 1 for local linear smoothing.
#' @param nder The order of derivative, which should be smaller than npoly. Default is 0 for local linear smoothing.
#' @param hkappa A numeric value (> 1) for neighborhood width adjustment in jump detection to avoid repeated identification of the same jump. Default is 2.
#' @param Refine A boolean indicating whether to perform the Refining Stage (iterative refining) in jump detection.
#' @param refine_tol The tolerance value for the Refining Stage. Default is 1e-4.
#' @param refine_maxiter The maximum number of iterations allowed in the Refining Stage. Default is 10.
#' @param jumpcheck A boolean indicating whether to double check if the sizes of detected jumps are above the threshold value (zeta).
#' @param verbose A boolean indicating whether to output more information during function execution.
#'
#' @return A list containing the tuning parameters used (h_tau, h_d, h_2, zeta, and rho_d) and the following fields: 
#' \item{alp_est}{A list of estimated varying coefficient functions, including detected change point locations (jumptime), estimated jump sizes (jumpsize), and the estimated coefficient function values on the inputs (alp.hat)}
#' \item{muout}{A vector of predicted response for the testing data (tout, xout).}
#' 

CoefJump <- function(tin, 
                     yin, 
                     xin, 
                     win = NULL, 
                     tout = NULL, 
                     xout = NULL,
                     h_tau, 
                     h_d, 
                     h_2,
                     zeta, 
                     kernel= c("epan", "rect", "gauss", "gausvar", "quar"),
                     NbGrid = 101,
                     npoly = 1, 
                     nder = 0, 
                     hkappa = 2, 
                     Refine = F, 
                     refine_tol = 1e-4, 
                     refine_maxiter = 10,
                     jumpcheck = T,
                     verbose = T) {
  
  kernel = match.arg(kernel)
  if (is.null(tout)) tout = tin
  
  ##### Searching Stage
  
  if (verbose) message('  Searching Stage starts...')
  
  n = dim(xin)[1]
  p = dim(xin)[2]
  m = dim(tin)[2]
  
  yyin <- yin
  
  res <- vector(mode = 'list', length = p+1)
  names(res) <- paste0('alpha', c(0:p))
  
  rho_d = NULL
  
  CP_set = matrix(data = NA, nrow = 0, ncol = 3)
  
  for (ell in 1:m) {
    
    tin_dm = tin
    tin_dm[, c(1, ell)] = tin[, c(ell, 1)]
    tord = order(tin_dm[, 1])
    tin_dm = tin_dm[tord, ,drop=F]
    xin_dm = xin[tord, ,drop=F]
    yin_dm = yin[tord]
    win_dm = win[tord]
    
    # Generate basic grids for jump detect:
    obsGrid = tin_dm;
    if(is.null(NbGrid)){
      jumpGrid = obsGrid
    } else {
      jumpGrid = apply(obsGrid, 2, function(x){
        seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
      })
      D_h <- jumpGrid[, 1]
      jumpGrid[, c(1, m)] = jumpGrid[, c(m, 1)]
      jumpGrid = as.matrix(expand.grid(as.data.frame(jumpGrid)))
      jumpGrid[, c(1, m)] = jumpGrid[, c(m, 1)] # to ensure the first column varies the slowest
    }
    
    ## Local linear estimate based on one-sided kernel
    ##
    alp_est_lr = CPPlwls2d_s2_LR(bw = h_tau, kernel_type = kernel, win = win_dm,
                                 tin = tin_dm, yin = yin_dm, xin = xin_dm, tout = jumpGrid, npoly = npoly)
    alp_est_left = alp_est_lr[,1:(p+1), drop=F]
    alp_est_right = alp_est_lr[,-c(1:(p+1)), drop=F]
    
    for (j in 1:(p+1)) {
      
      res[[j]][[ell]] = vector(mode = 'list', length = 3)
      names(res[[j]][[ell]]) = c('jumptime','jumpsize','jumpsize_h_tau')
      
      alpj_diff = alp_est_right[,j] - alp_est_left[,j]
      
      ###################
      ## Integrate over t\ell
      
      # alpj_diff_abs = sapply(1:NbGrid, function(ii) mean( alpj_diff_abs[jumpGrid[,1] == D_h[ii]] ))
      alpj_diff_abs = abs(sapply(1:NbGrid, function(ii) mean( alpj_diff[jumpGrid[,1] == D_h[ii]] )))
      alpj_diff_abs = abs(alpj_diff_abs)
      
      alpj_diff_time = D_h[order(-alpj_diff_abs)]
      alpj_diff_size = alpj_diff_abs[order(-alpj_diff_abs)]
      
      ## only focus on jumps size that are bigger than zeta
      timeindex = which(abs(alpj_diff_size) > zeta)
      ll=length(timeindex)
      
      if (ll > 0) {
        
        alpj_diff_time=alpj_diff_time[timeindex]
        alpj_diff_size=alpj_diff_size[timeindex]
        
        # find cluster centers
        alpj_jumptime=alpj_diff_time[1]
        alpj_jumpsize=alpj_diff_size[1]
        for (i in 2:ll){
          if (sum(abs(alpj_diff_time-alpj_jumptime[i-1]) > hkappa * h_tau) >0){ 
            index=which(abs(alpj_diff_time-alpj_jumptime[i-1]) > hkappa * h_tau)
            alpj_diff_time = alpj_diff_time[index]
            alpj_diff_size = alpj_diff_size[index]
            alpj_jumptime=append(alpj_jumptime, alpj_diff_time[1])
            alpj_jumpsize=append(alpj_jumpsize, alpj_diff_size[1])
          }
          else{
            break
          }
        }
        alpj_jumpsize_h_tau = alpj_jumpsize[order(alpj_jumptime)]
        alpj_jumptime = sort(alpj_jumptime)
        
        ## jumpcheck step to estimate jump size
        # rho_d = 1.1*h_tau^2
        rho_d = 0.25*(h_tau^2 + sqrt(log(n)/(n*h_tau))*h_tau)
        names(rho_d) = c("rho")
        jumpset_l = c(alpj_jumptime - rho_d)
        jumpset_r = c(alpj_jumptime + rho_d)
        othergrid = apply(obsGrid[, -c(1), drop=F], 2, function(x){
          seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
        }, simplify = F)
        jumpset_l_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_l))))
        jumpset_r_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_r))))
        if (m>1) {
          jumpset_l_grid = jumpset_l_grid[,c(m, 1:(m-1))]
          jumpset_r_grid = jumpset_r_grid[,c(m, 1:(m-1))]
        }
        
        ## jump size based on the bandwidth h_d
        alpj_est_lr_l = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                        tin = tin_dm, yin = yin_dm, xin = xin_dm, 
                                        tout = jumpset_l_grid,
                                        npoly = npoly)
        alpj_est_lr_r = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                        tin = tin_dm, yin = yin_dm, xin = xin_dm, 
                                        tout = jumpset_r_grid,
                                        npoly = npoly)
        alpj_jumpsize_h_d = alpj_est_lr_r[,p+1+j] - alpj_est_lr_l[,j]
        
        alpj_jumpsize_h_d = sapply(1:length(jumpset_l), 
                                   function(ii) mean( alpj_jumpsize_h_d[jumpGrid[,1] == D_h[ii]] ))
        
        
        ## jumpcheck step for the jump location
        if(jumpcheck){
          alpj_jumptime = alpj_jumptime[abs(alpj_jumpsize_h_d) > zeta]
          alpj_jumpsize_h_d = alpj_jumpsize_h_d[abs(alpj_jumpsize_h_d) > zeta]
        }
        
        if(length(alpj_jumptime)==0){
          alpj_jumptime = NULL
          # if (verbose) cat('alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
        } else {
          res[[j]][[ell]]$jumptime = alpj_jumptime
          res[[j]][[ell]]$jumpsize = alpj_jumpsize_h_d
          # if (verbose) cat('alpha_', j-1, 'no.', ell, 'dimension --- change point at:', 
          #                  alpj_jumptime, '  jump size:', alpj_jumpsize_h_d, '\n')
        }
        
      } else {
        alpj_jumptime = NULL
        # if (verbose) cat('alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
      }
      
      if (!is.null(alpj_jumptime)) {
        for (jumpid in 1:length(alpj_jumptime)) CP_set = rbind(CP_set, c(ell, j, alpj_jumptime[jumpid]))
      }
      
    }
    
  }
  
  ##### Refining Stage
  
  if (m == 1 && Refine) {
    Refine = F
    message('  Refining is not needed when m=1.')
  }
  
  if (Refine && nrow(CP_set) > 0) { 
    
    if (verbose) message('  Refining Stage starts...')
    
    for (refine_iter in 1:refine_maxiter) {
      
      if (verbose) message('    iteration ', refine_iter)
      
      CP_set_prev = CP_set
      
      for (ell in 1:m) {
        
        tin_dm = tin
        tin_dm[, c(1, ell)] = tin[, c(ell, 1)]
        tord = order(tin_dm[, 1])
        tin_dm = tin_dm[tord, ,drop=F]
        xin_dm = xin[tord, ,drop=F]
        yin_dm = yin[tord]
        win_dm = win[tord]
        
        CP_set = CP_set[CP_set[,1] != ell, , drop = F]
        CP_set_dm = CP_set
        CP_set_dm[CP_set_dm[,1] == 1, 1] = ell
        
        # Generate basic grids for jump detect:
        obsGrid = tin_dm;
        if(is.null(NbGrid)){
          jumpGrid = obsGrid
        } else {
          jumpGrid = apply(obsGrid, 2, function(x){
            seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
          })
          D_h <- jumpGrid[, 1]
          jumpGrid[, c(1, m)] = jumpGrid[, c(m, 1)]
          jumpGrid = as.matrix(expand.grid(as.data.frame(jumpGrid)))
          jumpGrid[, c(1, m)] = jumpGrid[, c(m, 1)] # to ensure the first column varies the slowest
        }
        
        ## Local linear estimate based on one-sided kernel
        ##
        alp_est_lr = CPPlwls2d_s2_LR(bw = h_tau, kernel_type = kernel, win = win_dm,
                                     tin = tin_dm, yin = yin_dm, xin = xin_dm, tout = jumpGrid, npoly = npoly)
        alp_est_left = alp_est_lr[,1:(p+1), drop=F]
        alp_est_right = alp_est_lr[,-c(1:(p+1)), drop=F]
        
        for (j in 1:(p+1)) {
          
          res[[j]][[ell]] = vector(mode = 'list', length = 2)
          names(res[[j]][[ell]]) = c('jumptime','jumpsize')
          
          alpj_diff = alp_est_right[,j] - alp_est_left[,j]
          
          ###################
          ## Integrate over t\ell except for the neighborhoods of change points on other dimensions
          
          alpj_diff_abs = sapply(1:NbGrid, function(ii) {
            
            jumpGrid_sub = jumpGrid[jumpGrid[,1] == D_h[ii],]
            alpj_diff_sub = alpj_diff[jumpGrid[,1] == D_h[ii]]
            
            for (ell_temp in 2:m) {
              if (length(which(CP_set_dm[,1] == ell_temp & CP_set_dm[,2] == j)) == 0) next
              jp_ell = CP_set_dm[CP_set_dm[,1] == ell_temp & CP_set_dm[,2] == j, 3]
              isext_fun = function(x) !any(x > jp_ell-h_tau & x < jp_ell+h_tau)
              sub_id = which(sapply(jumpGrid_sub[,ell_temp], isext_fun))
              jumpGrid_sub = jumpGrid_sub[sub_id,]
              alpj_diff_sub = alpj_diff_sub[sub_id]
            }
            
            if (length(alpj_diff_sub) < length(alpj_diff)/2) NaN
            else mean(alpj_diff_sub)
          })
          
          if (any(is.nan(alpj_diff_abs))) alpj_diff_abs = abs(sapply(1:NbGrid, function(ii) mean( alpj_diff[jumpGrid[,1] == D_h[ii]] )))
          
          alpj_diff_abs = abs(alpj_diff_abs)
          
          alpj_diff_time = D_h[order(-alpj_diff_abs)]
          alpj_diff_size = alpj_diff_abs[order(-alpj_diff_abs)]
          
          ## only focus on jumps size that are bigger than zeta
          timeindex = which(abs(alpj_diff_size) > zeta)
          ll=length(timeindex)
          
          if (ll > 0) {
            
            alpj_diff_time=alpj_diff_time[timeindex]
            alpj_diff_size=alpj_diff_size[timeindex]
            
            # find cluster centers
            alpj_jumptime=alpj_diff_time[1]
            alpj_jumpsize=alpj_diff_size[1]
            for (i in 2:ll){
              if (sum(abs(alpj_diff_time-alpj_jumptime[i-1]) > hkappa * h_tau) >0){ 
                index=which(abs(alpj_diff_time-alpj_jumptime[i-1]) > hkappa * h_tau)
                alpj_diff_time = alpj_diff_time[index]
                alpj_diff_size = alpj_diff_size[index]
                alpj_jumptime=append(alpj_jumptime, alpj_diff_time[1])
                alpj_jumpsize=append(alpj_jumpsize, alpj_diff_size[1])
              }
              else{
                break
              }
            }
            alpj_jumpsize_h_tau = alpj_jumpsize[order(alpj_jumptime)]
            alpj_jumptime = sort(alpj_jumptime)
            
            ## refine stage to estimate jump size
            # rho_d = 1.1*h_tau^2
            rho_d = 0.25*(h_tau^2 + sqrt(log(n)/(n*h_tau))*h_tau)
            names(rho_d) = c("rho")
            jumpset_l = c(alpj_jumptime - rho_d)
            jumpset_r = c(alpj_jumptime + rho_d)
            othergrid = apply(obsGrid[, -c(1), drop=F], 2, function(x){
              seq( max(min(x), h_tau), min(max(x), 1- h_tau), length.out = NbGrid)
            }, simplify = F)
            jumpset_l_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_l))))
            jumpset_r_grid = as.matrix(expand.grid(c(othergrid, list(jumpset_r))))
            if (m>1) {
              jumpset_l_grid = jumpset_l_grid[,c(m, 1:(m-1))]
              jumpset_r_grid = jumpset_r_grid[,c(m, 1:(m-1))]
            }
            
            ## jump size based on the bandwidth h_d
            alpj_est_lr_l = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                            tin = tin_dm, yin = yin_dm, xin = xin_dm, 
                                            tout = jumpset_l_grid,
                                            npoly = npoly)
            alpj_est_lr_r = CPPlwls2d_s2_LR(bw = h_d, kernel_type = kernel, win = win,
                                            tin = tin_dm, yin = yin_dm, xin = xin_dm, 
                                            tout = jumpset_r_grid,
                                            npoly = npoly)
            
            alpj_jumpsize_h_d = alpj_est_lr_r[,p+1+j] - alpj_est_lr_l[,j]
            
            alpj_jumpsize_h_d = sapply(1:length(jumpset_l),
                                       function(ii) mean( alpj_jumpsize_h_d[jumpGrid[,1] == D_h[ii]] ))
            
            ## refine the jump location
            if(jumpcheck){
              alpj_jumptime = alpj_jumptime[abs(alpj_jumpsize_h_d) > zeta]
              alpj_jumpsize_h_d = alpj_jumpsize_h_d[abs(alpj_jumpsize_h_d) > zeta]
            }
            
            if(length(alpj_jumptime)==0){
              alpj_jumptime = NULL
              res[[j]][[ell]]$jumptime = NULL
              res[[j]][[ell]]$jumpsize = NULL
              # if (verbose) cat('alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
            } else {
              res[[j]][[ell]]$jumptime = alpj_jumptime
              res[[j]][[ell]]$jumpsize = alpj_jumpsize_h_d
              # if (verbose) cat('alpha_', j-1, 'no.', ell, 'dimension --- change point at:', 
              #                  alpj_jumptime, '  jump size:', alpj_jumpsize_h_d, '\n')
            }
            
          } else {
            alpj_jumptime = NULL
            res[[j]][[ell]]$jumptime = NULL
            res[[j]][[ell]]$jumpsize = NULL
            # if (verbose) cat('alpha_', j-1, 'no.', ell, 'dimension --- no change point after refine.', '\n')
          }
          
          if (!is.null(alpj_jumptime)) {
            for (jumpid in 1:length(alpj_jumptime)) CP_set = rbind(CP_set, c(ell, j, alpj_jumptime[jumpid]))
          }
          
        }
        
      }
      
      if (pracma::hausdorff_dist(CP_set, CP_set_prev) < refine_tol) {
        message('    Refining stage successfully converges.')
        break
      }
      else if (refine_iter == refine_maxiter && verbose) message('  Refining stage fails to converge.')
      
    }
    
  }
  
  if (verbose) message('Jump detection results:')
  
  for (ell in 1:m){
    for (j in 1:(p+1)){
      
      alpj_jumptime = res[[j]][[ell]]$jumptime
      alpj_jumpsize_h_d = res[[j]][[ell]]$jumpsize
      
      if (is.null(alpj_jumptime)) {
        
        yyin <- yyin
        
        if (verbose) message('  alpha_', j-1, ' no. ', ell, ' dimension --- no jump.')
        
      } else{
        
        ysubstr <- sapply(tin[, ell], 
                          function(z) sum(alpj_jumpsize_h_d*(z >= alpj_jumptime)))
        yyin <- yyin - ysubstr * cbind(1,xin)[,j]
        
        if (verbose) message('  alpha_', j-1, ' no. ', ell, ' dimension --- jump location: ', 
                             paste(alpj_jumptime, collapse = ', '), '  jump size: ', paste(alpj_jumpsize_h_d, collapse = ', '))
      }
      
    }
  }
  
  
  # alp: the smoothed alpha curves evaluated at tin
  gamma = CPPlwls2d_s2(bw = h_2, kernel_type = kernel, win = win,
                       tin = tin, yin = yyin, xin = xin, tout = tin, npoly = npoly)
  
  # convert alpha to truncated tout
  gamma_out = CPPlwls2d_s2(bw = h_2, kernel_type = kernel, win = win,
                           tin = tin, yin = yyin, xin = xin, tout = tout, npoly = npoly)
  
  for (j in 1:(p+1)) {
    
    res[[j]]$gamma.hat <- gamma[,j]
    res[[j]]$gamma.hat_out <- gamma_out[,j]
    
    jump = rep(0, dim(tin)[1])
    jump_out = rep(0, dim(tout)[1])
    
    for (ell in 1:m) {
      if (!is.null(res[[j]][[ell]]$jumptime))
        jump = jump + sapply(tin[, ell], function(z) sum(res[[j]][[ell]]$jumpsize * (z >= res[[j]][[ell]]$jumptime)))
      if (!is.null(res[[j]][[ell]]$jumptime))
        jump_out = jump_out + sapply(tout[, ell], function(z) sum(res[[j]][[ell]]$jumpsize * (z >= res[[j]][[ell]]$jumptime)))
    }
    
    res[[j]]$alp.hat <- gamma[,j] + jump
    res[[j]]$alp.hat_out <- gamma_out[,j] + jump_out
    
  }
  
  
  if (is.null(xout)) {
    muout <- NULL
  } else {
    alp_out <- sapply(res, function(z) z$alp.hat_out)
    muout <- rowSums(alp_out * cbind(1, xout)) 
  }
  
  
  return(list(alp_est = res,
              muout = muout, 
              h_tau = h_tau, 
              h_d = h_d, 
              h_2 = h_2,
              zeta = zeta, 
              rho_d = rho_d))
}









######################## Create k folds

SimpleFolds <- function(yy, k=10) {
  if (length(yy) > 1)
    allSamp <- sample(yy)
  else
    allSamp <- yy
  
  n <- length(yy)
  nEach <- n %/% k
  samp <- list()
  length(samp) <- k
  for (i in seq_along(samp)) {
    if (nEach > 0)
      samp[[i]] <- allSamp[1:nEach + (i - 1) * nEach]
    else
      samp[[i]] <- numeric(0)
  }
  restSamp <- allSamp[seq(nEach * k + 1, length(allSamp), length.out=length(allSamp) - nEach * k)]
  restInd <- sample(k, length(restSamp))
  for (i in seq_along(restInd)) {
    sampInd <- restInd[i]
    samp[[sampInd]] <- c(samp[[sampInd]], restSamp[i])
  }
  
  return(samp)
}



######################## Calculate threshold for change point detection

zetaFun <- function(bw, 
                    alpha = 0.01, 
                    tin, 
                    yin_s2, 
                    xin, 
                    win, 
                    npoly = 1, 
                    nder = 0, 
                    kernel = 'epan',
                    nRegGrid = 101, 
                    cutoff = max,
                    test_point = NULL){
  
  n <- dim(xin)[1]
  m <- dim(tin)[2]
  
  ### epanechnikov kernel
  K_Epa = function(u) 0.75 * (1 - u^2) * (abs(u) <= 1);
  K_Epa_r = function(u, r) 0.75 * (1 - u^2) *u^r * (abs(u) <= 1);
  nu0 = stats::integrate(K_Epa_r, 0, 1, r=0)$value
  nu1 = stats::integrate(K_Epa_r, 0, 1, r=1)$value
  nu2 = stats::integrate(K_Epa_r, 0, 1, r=2)$value
  
  nu0a = stats::integrate(K_Epa_r, -1, 1, r=0)$value
  nu1a = stats::integrate(K_Epa_r, -1, 1, r=1)$value
  nu2a = stats::integrate(K_Epa_r, -1, 1, r=2)$value
  
  K_Epa2 = function(u) (0.75 * (1 - u^2))^2 * (abs(u) <= 1);
  K_Epa2_r = function(u, r) (0.75 * (1 - u^2))^2 *u^r * (abs(u) <= 1);
  mu0 = stats::integrate(K_Epa2_r, 0, 1, r=0)$value
  mu1 = stats::integrate(K_Epa2_r, 0, 1, r=1)$value
  mu2 = stats::integrate(K_Epa2_r, 0, 1, r=2)$value
  
  KK_star1 = function(u) ( (nu2*nu0a * nu0*nu2a - (nu1*nu1a)^2) * nu0a - 
                             (nu1*nu0a * nu0*nu2a - nu0*nu1a * nu1*nu1a) * nu0a * u - 
                             (nu0*nu1a * nu2*nu0a - nu1*nu0a * nu1*nu1a) * nu1a ) / 
    ( (nu2*nu0a * nu0*nu2a - (nu1*nu1a)^2) * nu0*nu0a - 
        (nu1*nu0a * nu0*nu2a - nu0*nu1a * nu1*nu1a) * nu1*nu0a - 
        (nu0*nu1a * nu2*nu0a - nu1*nu0a * nu1*nu1a) * nu0*nu1a ) *
    (0.75 * (1 - u^2) * (abs(u) <= 1))
  
  KK_star2 = function(u) (KK_star1(u))^2
  
  KK_term = stats::integrate(KK_star2, 0, 1)$value
  
  obsGrid = tin
  regGrid = seq( max(min(obsGrid[,1]), bw), min(max(obsGrid[,1]), 1-bw), length.out = nRegGrid)
  
  coef_hat <- CPPlwls2d_s2(bw = bw, 
                           tin = tin, yin = yin_s2, xin = xin, win = win, 
                           tout = obsGrid, kernel_type = kernel, 
                           npoly = npoly, nder = nder)
  yhat <- c(rowSums(coef_hat * cbind(1, xin)))
  sigma2 <- mean((yhat-yin_s2)^2)
  
  if (m == 1){
    df = stats::density(tin, kernel = 'rectangular')
    if (is.null(test_point)) {
      df_term =  cutoff( 1 / stats::approx(df$x, df$y, xout=regGrid)$y )
    } else {
      df_term = 1 / stats::approx(df$x, df$y, xout=test_point)$y
    }
  } else if (m == 2){
    df = MASS::kde2d(x = tin[,1], y = tin[,2],
                     h = c(bw, bw),
                     n = nRegGrid, lims = c(max(min(obsGrid[,1]), bw), min(max(obsGrid[,1]), 1-bw),
                                            max(min(obsGrid[,2]), bw), min(max(obsGrid[,2]), 1-bw)))
    df_inv = 1/df$z
    df_inv_int = apply(df_inv, 1, mean)
    if (is.null(test_point)) {
      df_term =  cutoff( df_inv_int )
    } else {
      df_term = stats::approx(df$x, df_inv_int, xout=test_point)$y
    }
  } else {
    df = stats::density(tin[,1], kernel = 'rectangular')
    if (is.null(test_point)) {
      df_term =  cutoff( 1 / stats::approx(df$x, df$y, xout=regGrid)$y )
    } else {
      df_term = 1 / stats::approx(df$x, df$y, xout=test_point)$y
    }
  }
  
  Gamma <- t(cbind(1,xin)) %*% cbind(1,xin) / n
  Gamma_term = cutoff(1/diag(Gamma))
  
  zeta <- stats::qnorm(1-alpha / 2) * sqrt( 2*sigma2 * df_term * Gamma_term * KK_term / (n*bw) + 4 * 1/n )  
  
  return(zeta)
  
}



######################## Cross validation bandwidth selection for Part 2

CVbandwidth <- function(bw.seq2 = NULL, 
                        zeta = NULL, 
                        win = win,
                        tin = tin, 
                        yin = yin, 
                        xin = xin,
                        npoly = 1, 
                        nder = 0, 
                        kernel = 'epan',
                        NbGrid, 
                        nRegGrid, 
                        jumpcheck = T, 
                        kFolds = 5, 
                        alpha = 0.01, 
                        cutoff = max,
                        hkappa = 2){
  
  m = dim(tin)[2]
  n = dim(xin)[1]
  p = dim(xin)[2]
  
  cv2 = array(Inf, dim = c(length(bw.seq2), length(bw.seq2), kFolds));
  zeta.seq = rep(NA, length(bw.seq2))
  
  theFolds = SimpleFolds(1:n, kFolds)
  
  # cat('\n\n#############', 'CV procedure: \n')
  
  for (j in 1:length(bw.seq2)) {
    
    # cat('h_tau = ', bw.seq2[j], '\n')
    
    if (is.null(zeta)) {
      ## obtain threshold value from the hypothesis test
      zeta.seq[j] = zetaFun(bw = bw.seq2[j], alpha = alpha,
                            tin = tin, yin = yin, xin = xin, win = win,
                            npoly = npoly, nder = nder, kernel = kernel,
                            nRegGrid = nRegGrid, cutoff = cutoff)
      
    } else zeta.seq[j] = zeta
    
    for (k in 1:length(bw.seq2)) {
      
      if (2*bw.seq2[k] <= 1.5*bw.seq2[j]) next
      # cat('----h_d = ', 2*bw.seq2[k], ': ')
      
      for (i in 1:kFolds) {
        
        ttest <- tin[theFolds[[i]], , drop=F]
        ytest <- yin[theFolds[[i]]]
        xtest <- xin[theFolds[[i]], , drop=F]
        ttrain <- tin[-theFolds[[i]], , drop=F]
        ytrain <- yin[-theFolds[[i]]]
        xtrain <- xin[-theFolds[[i]], , drop=F]
        wtrain <- win[-theFolds[[i]]]
        
        muout = tryCatch(
          CoefJump(tin = ttrain, yin = ytrain, xin = xtrain, win = wtrain, 
                   tout = ttest, xout = xtest,
                   h_tau = bw.seq2[j], h_d = 2*bw.seq2[k], zeta = zeta.seq[j], 
                   h_2 = sqrt(bw.seq2[j] * 2*bw.seq2[k]),  # bw.seq2[j] + 1/2 * (2*bw.seq2[k] - bw.seq2[j]),  # sqrt(bw.seq2[j] * 2*bw.seq2[k]),  
                   jumpcheck = jumpcheck, npoly=npoly, nder= nder, 
                   kernel = kernel,  NbGrid = NbGrid, 
                   hkappa = hkappa, verbose = F)$muout,
          error=function(err) {
            warning('Invalid bandwidth during stage 2 CV. Try enlarging the window size. h_tau=', bw.seq2[j], 
                    ' h_d=', bw.seq2[k], ' the ', i, '-th fold \n')
            return(Inf)
          })
        nan_rate <- length(which(is.nan(muout))) / length(ytest)
        if (nan_rate > 0.2) warning('Bandwidth: h_tau=', bw.seq2[j], ' h_d=', 2*bw.seq2[k], 'too small for varying coef CV.')
        tid <- which(!is.nan(muout))
        
        cv2[j,k,i] = sum((ytest[tid] - muout[tid])^2)
        # print(cv2)
        if(is.na(cv2[j,k,i]) || nan_rate > 0.5){
          cv2[j,k,i] = Inf;
        }
        
        # cat('MSE: ', mean(ytest[tid] - muout[tid])^2, '\n')
        # cat('==')
        
        # cat('\n')
      }
    }
  }
  
  
  if(min(cv2) == Inf){
    stop("All bandwidths resulted in infinite CV costs. (Part 2)")
  }
  
  cvMean = apply(cv2, c(1,2), mean)
  cvMeanid = which(cvMean == min(cvMean), arr.ind=TRUE)
  if (length(dim(cvMeanid))>1 && dim(cvMeanid)[1]>1) cvMeanid = cvMeanid[floor(dim(cvMeanid)[1]/2),] # cvMeanid[dim(cvMeanid)[1],]
  bopt2 = bw.seq2[cvMeanid[c(1,2)]];
  names(bopt2) = c('h_tau', 'h_d')
  bopt2['h_d'] = 2 * bopt2['h_d']
  zeta = zeta.seq[cvMeanid[1]]
  
  boptList <- list('bopt' = bopt2, 'zeta' = zeta, 'cvMean' = cvMean)
  
  return(boptList)
  
}

