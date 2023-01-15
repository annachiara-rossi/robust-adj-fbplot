# The minimum regularized covariance determinant (MRCD) is a robust
# estimator for multivariate location and scatter, which detects outliers by fitting a
# robust covariance matrix to the data. The MRCD assumes that the observations
# are elliptically distributed. However, this property does not always apply to modern datasets. 
# Together with the time criticality of industrial processing, small n,
# large p problems pose a challenging problem for any data analytics procedure. Both
# shortcomings are solved with the proposed kernel Minimum Regularized Covariance
# Determinant estimator, where we exploit the kernel trick to speed-up computations.
# More specifically, the MRCD location and scatter matrix estimate are computed in
# a kernel induced feature space, where regularization ensures that the covariance matrix is well-conditioned, 
# invertible and defined for any dataset dimension.

# inputs:
# - x:      a vector or matrix whose columns represent variables, and rows represent 
#           observations. Missing and infinite values are not allowed
#           and should be excluded from the computations by the user.
# - alpha:  (1-alpha) measures the fraction of outliers the algorithm should resist.
#           It should be in [0.5 <= alpha <= 1]. (default = 0.75)
# - kModel: the kernel transformation used (default = Linear -> kMRCD = MRCD)

# outputs:
# - solution.outlyingnessIndices: outlyingness weight of each observation according.
# - solution.name: the name of the best initial estimator which was used to construct 
#   the final solution.
# - solution.hsubsetIndices: the h-subset element indices after C-step convergence.
# - solution.obj: the MRCD objective value.
# - solution.smd: the Squared Mahanobis Distance values of each observation
# - solution.rho: regularisation factor used
# - solution.scfac: finite sample correction factor used
# - solution.rd: Malahanobis distance value of each observation.
# - solution.ld: rescaled Malahanobis distance values, defined as log(0.1 + solution.rd);
# - solution.cutoff: outlier flagging cut-off
# - solution.flaggedOutlierIndices: indices of the flagged outliers

# set working directory
# setwd("~/Documents/Thesis work")
source(here::here("kMRCD/utils.R"))

CovKmrcd <- function(x, kModel = 'Linear', alpha = 0.75){
  
  # Maximum number of CStep iterations allowed 
  cStepIterationsAllowed = 100
  
  # Condition number one wants to achieve
  maxcond = 50
  
  # be sure the value of alpha is acceptable, otherwise stop  
  testit::assert('The percentage of regular observations, alpha, should be in [0.5-1]', 
                 alpha<=1 && alpha>=0.5)
  
  print(paste("Using", kModel, "kernel"))
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  if (kModel == 'RbfKernel'){
    K <- rdetools::rbfkernel(x)        # kernel matrix nxn
  }
  else{
    K <- linearkernel(x)
  }
  
  # choose some initial estimators and compute rho from each one of those
  # these should be methods that can be kernelized
  
  solution = list(
    list(outlyingnessIndices = SDO(K,alpha)$h_idx,
         name = 'SDO'
    ),
    list(outlyingnessIndices = SpatialRank(K,alpha)$h_idx,
         name = 'SpatialRank'
    ),
    list(outlyingnessIndices = SpatialMedianEstimator(K,alpha)$h_idx,
         name = 'SpatialMedian'
    ),
    list(outlyingnessIndices = SSCM(K)$h_idx,
         name = 'SSCM'
    )
  )
  # solution[[1]]$..  to access elements
  
  scfac = MCDcons(p, alpha)
  
  n_sol = pracma::numel(solution)
  rhoL = rep(NaN,n_sol)
  # refine each initial estimate
  for (index in 1:n_sol){
    # Compute solution for each estimator
    solution[[index]]$hsubsetIndices = 
      solution[[index]]$outlyingnessIndices[1:ceiling(n*alpha)]
    
    # Determine rho for each estimator 
    
    if (kModel == 'RbfKernel'){
      s = svd(center( rdetools::rbfkernel(x[solution[[index]]$hsubsetIndices,]) 
      ))$d
    }
    else{
      s = svd(center( linearkernel(x[solution[[index]]$hsubsetIndices,]) 
      ))$d
    }
    
    nx = length(solution[[index]]$hsubsetIndices)                
    e_min = min(s) 
    e_max = max(s)
    
    fncond = function(rho) (nx*rho + (1-rho)*scfac%*%e_max)/(nx*rho + (1-rho)*scfac%*%e_min) - maxcond
    
    tryCatch(
      {rhoL[index] = pracma::fzero(fncond, c(10^(-6),0.99))$x},
      # if that doesn't work, use:
      error = function(e) {
                           grid = seq(0.000001,1-0.000001,length.out=1000)
                           objgrid = abs(sapply(grid,fncond))
                           rhoL[index] = min(grid[which(objgrid == min(objgrid))])
              }
    )

  }
  
  # Determine the final rho as in Boudt et al.
  # Set rho as max of the rho_i's obtained for each subset in previous step                        
  rho = max(rhoL[ rhoL <= max(0.1, median(rhoL)) ])
  
  # Refine each initial estimation with C-steps
  Ktt_diag = diag(K)
  for (index in (1:n_sol)){                           
    for (iteration in 1:cStepIterationsAllowed){
      
      hSubset = solution[[index]]$hsubsetIndices   
      Kx = if (kModel == 'RbfKernel') rdetools::rbfkernel(x[hSubset,]) else linearkernel(x[hSubset,])
      nx = dim(Kx)[1]
      Kt = if (kModel == 'RbfKernel') rdetools::rbfkernel(X = x, Y = x[hSubset,]) else linearkernel(X = x, Y = x[hSubset,])
      Kc = center(Kx)
      Kt_c = center(Kx,Kt)
      Kxx = Ktt_diag - (2/nx)*rowSums(Kt) + rep((1/nx^2)*sum(Kx),n)
      # TODO: it seems to me that the target matrix only appears here -> substitute it?
      smd = (1/rho)*(Kxx - (1-rho)*scfac*rowSums((pracma::mrdivide(Kt_c,((1-rho)*scfac*Kc + nx*rho*pracma::eye(nx)))*Kt_c)))                    
     
      indices = sort(smd, index.return= TRUE)$ix
      # Redefine the h-subset
      solution[[index]]$hsubsetIndices = indices[1:nx]                
      if (pracma::isempty(setdiff(hSubset, solution[[index]]$hsubsetIndices))){ # prob::
        print(paste('Convergence at iteration', iteration,',', solution[[index]]$name ))                                               
        sigma = svd(Kc)$d
        sigma = (1-rho)*scfac%*%sigma + pracma::numel(solution[[index]]$hsubsetIndices)*rho                        
        solution[[index]]$obj =  sum(log(sigma))
        solution[[index]]$smd = smd
        break
      }
      
    }
    
    testit::assert('no C-step convergence', iteration < cStepIterationsAllowed)
  }
  
  # Select the solution with the lowest objective function  
  mIndex = 1
  for (i in 2:n_sol){
    if (solution[[i]]$obj < solution[[mIndex]]$obj){
      mIndex = i
    }
  }
  
  # ...and remove the other solutions
  solution = solution[[mIndex]]                            
  print(paste('-> Best estimator is ',solution$name))
  
  solution$rho = rho
  solution$scfac = scfac
  
  # TODO: compute the regularized covariance matrix and test it
  hSubset = solution$hsubsetIndices
  x_c = scale(x[hSubset,], scale = FALSE)
  
  Kx = if (kModel == 'RbfKernel') rdetools::rbfkernel(x_c) else linearkernel(x_c)
  nx = dim(Kx)[1]     # this is h = ceiling(n*alpha)
  Kc = center(Kx)
  svd_decomp = svd(scfac*Kc)
  D = diag(svd_decomp$d)
  Vc = scale(svd_decomp$v, scale = FALSE)
  
  Sigma_reg = (1-rho)/(nx) * t(x_c)%*%Vc%*%D%*%t(Vc)%*%x_c + rho*pracma::eye(p)
  solution$cov = Sigma_reg
  
  # Outlier flagging procedure 
  solution$rd = sqrt(pmax(solution$smd,0))     # it was max(sqrt(solution$smd),0) but square root is always positive
  solution$ld = log(0.1 + solution$rd)
  
  # use this if decide to implement unimcd
  #standardization = unimcd(solution$ld, pracma::numel(solution$hsubsetIndices))
  #tmcd = standardization$tmcd
  #smcd = standardization$smcd
  estim = robustbase::covMcd(solution$ld, alpha)
  tmcd = estim$center
  smcd = estim$cov
  
  solution$cutoff = exp(tmcd + qnorm(0.995) * smcd) - 0.1
  solution$flaggedOutlierIndices = which(solution$rd > solution$cutoff[1])
  
  return (solution)
  
}








