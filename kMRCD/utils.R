
# Count the number of inputs to a function
nargin <- function() {
  if(sys.nframe()<2) stop("must be called from inside a function")
  length(as.list(sys.call(-1)))-1
}

# Centering of the kernel matrix
center <- function(omega, Kt){
  # omega needs to be a square matrix!
  nb_data = dim(omega)[1]
  meanvec = rowMeans(omega)
  MM = mean(meanvec)
  
  if (nargin()<2){
    Kc= omega - meanvec%*%t(rep(1,nb_data)) - rep(1,nb_data)%*%t(meanvec) + MM
  }
  else{
    nt = dim(Kt)[1]
    meanvecT= rowMeans(Kt)
    Kc= Kt - rep(1,nt)%*%t(meanvec) - meanvecT%*%t(rep(1,nb_data)) + MM
  }
  
  return (Kc)
}

# Scale computation
W_scale <- function(x){
  # x: 2D matrix
  n = dim(x)[1]
  p = dim(x)[2]
  
  # robust estimate of the standard deviation and mean of each column
  sigma0 = apply(x,FUN=mad,MARGIN=2, constant = 1)
  med0 = apply(x,FUN=median,MARGIN=2)
  
  # create anonymous function
  Wc = function(x) (1-(x/4.5)^2)^2 * (abs(x)<4.5)
  w = Wc(x - pracma::repmat(med0,n,1)) / pracma::repmat(sigma0,n,1)
  loc=t(diag(t(x)%*%w))/apply(w,FUN = sum, MARGIN = 2)
  
  rc= function(x) (min(x^2,3^2)) 
  b=3*qnorm(3/4)
  
  nes=n*(2*((1-b^2)*pnorm(b)-b*dnorm(b)+b^2)-1)
  val = apply( (x-pracma::repmat(loc,n,1))/pracma::repmat(sigma0,n,1), rc, MARGIN = c(1,2) )
  scale=sigma0^2/nes*apply(val,FUN = sum, MARGIN = 2)
  return(sqrt(scale))
}


MCDcons <- function(p, alpha){ 
  qalpha = qchisq(alpha, df = p)
  caI = pgamma(qalpha/2, shape = p/2 + 1, scale = 1) / alpha
  
  return(1/caI)
}

# spatial median computation
SpatialMedian <- function(K){
  testit::assert('The kernel matrix should be a square matrix', 
                 dim(K)[1] == dim(K)[2])
  n = dim(K)[2]
  testit::assert('error: no samples', n > 0)
  
  gamma = rep(1/n,n)
  nn = rep(1,n)
  for (i in 1:10){
    w = nn / sqrt( diag(K) - 2*t(K)%*%gamma + rep(t(gamma)%*%K%*%gamma,n) )   # elementwise
    gamma = w/sum(w)  # elementwise
  }
  return (gamma)
}

##### unimcd <- function(y,h,varargin){####
#   
#   if (nargin() > 2){
#     centered=varargin[1]
#   }
#   else
#     centered=0
#   
#   ncas=length(y)
#   len=ncas-h+1
#   xorig=y
#   
#   if (len==1){
#     if (centered==1){
#       tmcd=0
#     }
#     else{
#       tmcd=mean(y)
#     }
#     smcd=sqrt(var(y))
#     weights=ones(length(y),1)
#   }
#   
#   else{
#     tmcd = c()
#     smcd = c()
#     for (j in 1:dim(y)[2]){
#       est = robustbase::covMcd(y[,j], h/n)
#       tmcd = c(tmcd, est$center)
#       smcd = c(smcd, est$cov)
#     }
#     
#   }
#   return (list(tmcd = tmcd, smcd = smcd))
# }


################# REFINEMENT STEP ##################
refinement_step <- function(h_idx, K){
  
  n = dim(K)[1]
  n_h = length(h_idx)
  K_h = K[h_idx,h_idx]
  Kt = K[ ,h_idx]
  
  # Compute covariance matrix
  K_tilde = center(K_h)
  decomp = svd(K_tilde)
  U = decomp$u
  S_F = decomp$d
  
  eps = pracma::eps(1.0)  # Numerical stability
  mask = S_F > 1000*eps 
  
  U = U[,mask]
  S_F = S_F[mask]
  
  if (is.null(dim(U))){
    U_scaled = U / pracma::repmat(sqrt(t(S_F)),length(U),1 )
  }else{
    U_scaled = U / pracma::repmat(sqrt(t(S_F)),dim(U)[1],1 )
  }
  
  
  # Step 1
  o = rep(1,n_h)
  gamma = o/n_h
  K_Phi_PhiTilde = (Kt -  Kt %*% gamma %*% t(o))
  B_F = K_Phi_PhiTilde %*% U_scaled
  lambda_F = W_scale(B_F)
  
  # Step 2: estimate the center
  K_Adapted = K_Phi_PhiTilde%*%U_scaled%*%pracma::Diag(lambda_F^(-1))%*%t(U_scaled)%*%t(K_Phi_PhiTilde)
  gamma_c = SpatialMedian(K_Adapted)
  
  # Step 3: Calculate Mahalanobis
  on = rep(1,n)
  Kt_cCov = Kt - on %*% t(gamma_c) %*% Kt - Kt %*% gamma %*% t(o) + 
            pracma::repmat(t(gamma_c) %*% Kt %*% gamma,n,n_h) * (on %*% t(o))
  mahal_F = rowSums((Kt_cCov%*%U_scaled%*%pracma::Diag(lambda_F^(-2))) * (Kt_cCov%*%U_scaled))
  
  h_idx <- sort(mahal_F, index.return= TRUE)$ix
  
  return (h_idx)
}


################ INITIAL ESTIMATORS ################

SpatialMedianEstimator <- function(K,h){
  testit::assert('The kernel matrix should be a square matrix', 
                 dim(K)[1] == dim(K)[2])
  n = dim(K)[1]
  
  gamma = SpatialMedian(K)
  
  # Calculate Euclidean distance of the elements of K to spatial median
  dist = diag(K) + pracma::repmat(t(gamma)%*%K%*%gamma,n,1) - 
          2*rowSums(K*pracma::repmat(t(gamma),n,1))
  # take the h indexes with the lowest dist
  # but I want the indexes in the original data corresponding to these values
  h_idx <- sort(dist, index.return= TRUE)$ix
  
  h_idx <- refinement_step(h_idx[1:ceiling(n*h)], K)
  
  return (list(h_idx = h_idx, 
               dist = dist, 
               gamma = gamma))
}

SDO <- function(K,h){
  testit::assert('The kernel matrix should be a square matrix', 
                 dim(K)[1] == dim(K)[2])
  n = dim(K)[1]
  gamma = rep(0,n)
  
  for (direction in 1:500){
    rindices = rep(0,2)
    
    while (rindices[1] == rindices[2]) 
      rindices  = pracma::randperm(n, 2)
    
    lambda = rep(0,n)
    lambda[rindices[1]]=1
    lambda[rindices[2]]=-1
    
    a = K %*% lambda / pracma::repmat(sqrt(t(lambda)%*%K%*%lambda),n,1)
    sdo = abs(a - median(a)) / DescTools::MeanAD(a)
    mask = sdo>gamma
    gamma[mask] = sdo[mask]
  }                
  
  h_idx <- sort(gamma, index.return= TRUE)$ix
  h_idx <- refinement_step(h_idx[1:ceiling(n*h)], K)
  
  return(list(h_idx = h_idx, 
              gamma = gamma))
}

SpatialRank <- function(K,h){
  testit::assert('The kernel matrix should be a square matrix', 
                 dim(K)[1] == dim(K)[2])
  n = dim(K)[1]
  ook = rep(0,n)
  
  for (k in 1:n){
    tmpA = K[k,k] - t(pracma::repmat(K[, k], n,1) ) - pracma::repmat(K[k, ], n, 1) + K
    tmpB = sqrt(K[k,k] + t(diag(K)) - 2*K[k,])
    tmpC = pracma::repmat(tmpB, n, 1) * pracma::repmat(t(tmpB), 1, n)               
    mask = matrix(rep(TRUE,n*n),nrow = n)
    mask[k,] = FALSE
    mask[,k] = FALSE       
    ook[k] = sum(tmpA[mask]/tmpC[mask])
  }
  
  ook = (1/n)*sqrt(ook)
  
  h_idx <- sort(ook, index.return= TRUE)$ix
  h_idx <- refinement_step(h_idx[1:ceiling(n*h)], K)
  
  return(list(h_idx = h_idx, 
              ook = ook))
}

SSCM <- function(K){
  testit::assert('The kernel matrix should be a square matrix', 
                 dim(K)[1] == dim(K)[2])
  n = dim(K)[1]
  # Center kernel with spatial median
  gamma = SpatialMedian(K)
  
  # Spatial Sign Covariance matrix
  o = rep(1,n)           
  Kc = K - o %*% t(gamma) %*% K - K %*% gamma %*% t(o) + pracma::repmat(t(gamma) %*% K %*% gamma,n,n) * (o %*% t(o))
  
  D = diag(as.vector((diag(K) - 2*rowSums(K*pracma::repmat(t(gamma),n,1)) + pracma::repmat(t(gamma)%*%K%*%gamma,n,1))^(-1)))
  
  K_tilde = D^(1/2)%*%Kc%*%D^(1/2)   # until here everything's fine
  decomp = svd(K_tilde)
  U = decomp$u            # from the 1'th column on all signs are swapped
  S_F = decomp$d

  eps = pracma::eps(1.0)
  mask = S_F > 1000*eps   # here I get 116 elements instead of 111
  
  U = U[,mask]
  S_F = S_F[mask]
  U_scaled = U / pracma::repmat(sqrt(t(S_F)), dim(U)[1], 1)
            
  # Step 1
  K_Phi_PhiTilde = (K - K %*% gamma %*% t(o))%*%D^(1/2)
  
  B_F = K_Phi_PhiTilde %*% U_scaled
  lambda_F= W_scale(B_F)
  
  # Step 2: Estimate the center
  K_Adapted = K_Phi_PhiTilde %*% U_scaled %*% diag(lambda_F^(-1))%*% 
              t(U_scaled) %*% t(K_Phi_PhiTilde)
  gamma_c = SpatialMedian(K_Adapted)
  
  # Step 3: Calculate Mahalanobis
  K_cCov = K - o %*% t(gamma_c) %*% K - K %*% gamma %*% t(o) + 
            pracma::repmat(t(gamma) %*% K %*% gamma_c,n,n) * (o %*% t(o))
  
  K_cCov = K_cCov %*% D^(1/2)

  mahal_F = rowSums((K_cCov %*% U_scaled %*% diag(lambda_F^(-2))) * (K_cCov %*% U_scaled))
  # slightly different values
  
  h_idx <- sort(mahal_F, index.return= TRUE)$ix
  
  return(list(h_idx = h_idx, 
              gamma = gamma))
  
}


linearkernel <- function(X, Y = NULL){
  if(!is.matrix(X))
  {
    print("X will be transformed to a matrix containing samples in its rows")
    X = as.matrix(X)
  }
  
  if(is.null(Y)){
    return(X%*%t(X))
  }
  else{
    return(X%*%t(Y))
  }
  
}
