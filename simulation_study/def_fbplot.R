####________________ FINAL fbplot IMPLEMENTATION _______________####

# import Tarabelloni's code
source( 'rfda/basisOp/functional_basis.R')
source( 'rfda/geometry/geometry.R')
source( 'rfda/fdata/fData.R')
source( 'rfda/stats/stats.R' )
source( 'rfda/stats/generateFD.R')
source( 'rfda/stats/expCovariance.R')
source( 'rfda/stats/robustEstimators.R' )

# import roahd depths computation
source('roahd/band_depths.R')

# needed for KFPCA
require(interp)

# needed for kMRCD
source('kMRCD/kMRCD.R')

### computation of fences - same for all Covariance estimators
.my_fbplot_fData = function( time_grid, Data, Depths = 'MBD', Fvalue = 1.5 )
{
  # NB: takes in input the coefficients of fData
  
  # Number of observations
  N = nrow( Data )
  
  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definiton of
    # depth inside the vector that supposedly should contain depth values
    Depths = eval( parse( text = paste( Depths, '( Data )', sep = '' ) ) )
  } else {
    stopifnot( length( Depths ) == N )
  }
  
  Data_center = Data[ which.max( Depths ), ]
  
  id_central_region = which( Depths >= stats::quantile( Depths, prob = 0.5 ) )
  
  max_envelope_central = apply( Data[ id_central_region, ], 2, max )
  min_envelope_central = apply( Data[ id_central_region, ], 2, min )
  
  fence_upper = ( max_envelope_central - min_envelope_central ) * Fvalue + max_envelope_central
  fence_lower = ( min_envelope_central - max_envelope_central ) * Fvalue + min_envelope_central
  
  ID_outlying = which( apply( Data, 1, function(x)( any( x > fence_upper  ) |
                                                      any ( x < fence_lower ) ) ) )
  
  return( list( ID_out = ID_outlying,
                min_envelope_central = min_envelope_central,
                max_envelope_central = max_envelope_central,
                fence_lower = fence_lower,
                fence_upper = fence_upper ) )
}


### F* tuning using two different approaches: 
# 1) multivariate estimators of the covariance structure
# 2) functional operators whose eigenfunctions estimate the ones of the original
#    covariance structure

my_fbplot.fData = function( Data,
                            Cov_estimator,
                            mcd_control = list(TRUE,FALSE),
                            mrcd_target = "equicorrelation",
                            h = NULL,
                            Depths = 'MBD',
                            L = 10, # number of components to retain
                            Fvalue = 1.5,
                            adjust = FALSE,
                            display = TRUE,
                            xlab = NULL,
                            ylab = NULL,
                            main = NULL,
                            ... ){
  # Data: an f_data object from rfda package
  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definition of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths
    
    if( Depths_spec == 'MBD' )
    {
      Depths = MBD( expand(Data), manage_ties = TRUE )
      #Depths = MBD( Data$coefficients, manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( expand(Data) )',
                                          sep = '' ) ) )
      #Depths = eval( parse( text = paste( Depths, '( Data$coefficients )',
      #                                    sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == Data$N )
  }
  
  if( ! is.list( adjust ) )
  {
    # Plain functional boxplot with default F value: F = 1.5
    out = .my_fbplot_fData( time_grid, expand(Data), Depths, Fvalue )
    
  } else {
    
    nodenames = c( 'N_trials', 'trial_size', 'TPR', 'F_min', 'F_max',
                   'tol', 'maxiter', 'VERBOSE' )
    unused = setdiff( names( adjust ), nodenames )
    
    # Checking for unused parameters
    if( length( unused ) > 0 )
    {
      for( i in unused )
        warning( 'Warning: unused parameter ', i, ' in adjust argument of fbplot' )
    }
    
    N_trials = ifelse( is.null( adjust$N_trials ),
                       20,
                       adjust$N_trials )
    
    trial_size = ifelse( is.null( adjust$trial_size ),
                         8 * Data$N,
                         adjust$trial_size )
    
    TPR = ifelse( is.null( adjust$TPR ),
                  2 * stats::pnorm( 4 * stats::qnorm( 0.25 ) ),
                  adjust$TPR )
    
    F_min = ifelse( is.null( adjust$F_min ),
                    0.5,
                    adjust$F_min )
    
    F_max= ifelse( is.null( adjust$F_max ),
                   5,
                   adjust$F_max )
    
    tol = ifelse( is.null( adjust$tol ),
                  1e-3,
                  adjust$tol )
    
    maxiter = ifelse( is.null( adjust$maxiter ),
                      100,
                      adjust$maxiter )
    
    VERBOSE = ifelse( is.null( adjust$VERBOSE ),
                      FALSE,
                      adjust$VERBOSE )
    
    # compute Covariance estimator
    switch(Cov_estimator,
           
           # multivariate estimators
           Ledoit_Wolf = {
             Cov = cvCovEst::linearShrinkLWEst(expand(Data))
           },
           OGK_mad = {
             Cov = robustbase::covOGK(expand(Data), sigmamu = robustbase::s_mad )$cov
           },
           OGK_Qn = {
             Cov = robustbase::covOGK(expand(Data), sigmamu = robustbase::s_Qn )$cov
           },
           MCD = {
             # by default this uses corrections for finite samples
             Cov = robustbase::covMcd(x = expand(Data), alpha = alpha, nsamp = "best", 
                                      use.correction = mcd_control[[1]], 
                                      raw.only =  mcd_control[[2]])$cov
           },
           MRCD_linear = {
             # identity target
             Cov = rrcov::CovMrcd(x = expand(Data), alpha = 0.75)$cov
           },
           MRCD_0.5 = {
             # as a default use finite sample corrections and reweighting
             if (is.null(h)){
               Cov = rrcov::CovMrcd(x = expand(Data), target=mrcd_target, alpha = 0.5)$cov
             }
             else{
               Cov = rrcov::CovMrcd(x = expand(Data), target=mrcd_target,  h = h)$cov
             }
           },
           MRCD_0.75 = {
             # as a default use finite sample corrections and reweighting
             if (is.null(h)){
               Cov = rrcov::CovMrcd(x = expand(Data), target=mrcd_target, alpha = 0.75)$cov
             }
             else{
               Cov = rrcov::CovMrcd(x = expand(Data), target=mrcd_target,  h = h)$cov
             }
           },
           kMRCD = {
           Cov = CovKmrcd(x = expand(Data), alpha = .75)$cov
            },
           
           # functional operators
           Spherical = {
             Cov = covSph( Data, ASG = TRUE )   
           },
           Median = {
             Cov = covMed( Data, ASG = TRUE )   
           },
           Kendall = {
             # in this case I don't have an estimate of the covariance matrix,
             # but directly of the eigenfunctions
             Cov = NULL
           } 
    )
    
    # check if the estimator is multivariate or functional
    # because the two procedures are different
    
    # TODO: add more estimators when they are ready
    multivariate_est <- list("Ledoit_Wolf", "OGK_mad", "OGK_Qn", "MCD", 
                             "MRCD_linear", "kMRCD", "MRCD_0.5", "MRCD_0.75")
    functional_est <- list("Spherical", "Median", "Kendall")
    
    # starting the tuning procedure. The quantity to be minimized is the same
    N = nrow( Data$coefficients )
    
    Fvalues = rep( 0, N_trials )
    
    # TODO: in the cost functional Data_gauss should also be expanded?
    cost_functional = function( F_curr )( length(
      .my_fbplot_fData( time_grid,
                        Data_gauss,
                        Depths = Depths_spec,
                        Fvalue = F_curr )$ID_out ) /
        trial_size - TPR )
    
    
    # do this before the for loop, since they're fixed
    if( Cov_estimator %in% multivariate_est ){
      # Cholesky factor
      CholCov <- chol( Cov )
      
      # Centerline of the dataset
      centerline = expand(Data)[ which.max( Depths ), ]
    }
    
    
    if( Cov_estimator %in% functional_est ){
      ecouples = eigencouples_estimation(Data, Cov, L, N)
    }
    
    
    for( iTrial in 1 : N_trials )
    {
      if( VERBOSE > 0 )
      {
        cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
      }
      
      # gaussian data generation
      if( Cov_estimator %in% multivariate_est ){
        Data_gauss = roahd::generate_gauss_fdata( trial_size, centerline, 
                                                  CholCov = CholCov )
      }
      if( Cov_estimator %in% functional_est ){
        Data_gauss = expand(KL_data_generation(Data, ecouples, L, N, trial_size))
      }
      
      if( VERBOSE > 0 )
      {
        cat( ' * * * * beginning optimization\n' )
      }
      
      opt = stats::uniroot( cost_functional,
                            interval = c( F_min, F_max ),
                            tol = tol,
                            maxiter = maxiter, extendInt = "yes" )
      if( VERBOSE > 0 )
      {
        cat( ' * * * * optimization finished.\n')
      }
      
      Fvalues[ iTrial ] = opt$root
    }
    
    Fvalue = mean( Fvalues )
    
    # TODO: check this!!
    out = .my_fbplot_fData( time_grid, expand(Data), Depths, Fvalue = Fvalue  )
  }
  
  ID_out = out$ID_out
  
  # Plotting part
  if( is.numeric( display ) )
  {
    grDevices::dev.set( display )
  }
  
  if( ! display == FALSE )
  {
    # Creating color palettes
    col_non_outlying = scales::hue_pal( h = c( 180, 270 ),
                                        l = 60 )( Data$N - length( ID_out ) )
    
    col_non_outlying = roahd::set_alpha( col_non_outlying, 0.5 )
    
    if( length( ID_out ) > 0 )
    {
      col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                      c = 150 )( length( ID_out ) )
    } else {
      col_outlying = scales::hue_pal( h = c( - 90, 180  ),
                                      c = 150 )( 1 )
    }
    
    col_envelope = roahd::set_alpha( 'blue', alpha = 0.4 )
    col_center = roahd::set_alpha( 'blue', alpha = 1 )
    col_fence_structure = 'darkblue'
    
    time_grid = Data$geometry$basis$timeStr$grid
    
    xlab = ifelse( is.null( xlab ), '', xlab )
    ylab = ifelse( is.null( ylab ), '', ylab )
    main = ifelse( is.null( main ), '', main )
    
    if( length( ID_out ) > 0 )
    {
      # Plotting non-outlying data
      graphics::matplot( time_grid,
                         t( expand(Data)[ - ID_out, ] ), lty = 1, type = 'l',
                         col = col_non_outlying,
                         ylim = range( expand(Data) ),
                         xlab = xlab, ylab = ylab, main = main, cex.main = 2, ... )
      
      # Computing maximum and minimum envelope
      max_envelope_limit = apply( expand(Data)[ - ID_out, ], 2, max )
      min_envelope_limit = apply( expand(Data)[ - ID_out, ], 2, min )
    } else {
      # Plotting all data
      graphics::matplot( time_grid,
                         t( expand(Data) ), lty = 1, type = 'l',
                         col = col_non_outlying,
                         ylim = range( expand(Data) ),
                         xlab = xlab, ylab = ylab, main = main, cex.main = 2, ... )
      
      # Computing maximum and minimum envelope
      max_envelope_limit = apply( expand(Data), 2, max )
      min_envelope_limit = apply( expand(Data), 2, min )
    }
    
    # Filling in the central envelope
    
    graphics::polygon( c(time_grid, rev( time_grid) ),
                       c( out$min_envelope_central, rev( out$max_envelope_central ) ),
                       col = col_envelope, border = NA)
    graphics::lines( time_grid, out$max_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    graphics::lines( time_grid, out$min_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    
    # Plotting the sample median
    graphics::lines( time_grid, expand(Data)[ which.max( Depths ), ], lty = 1, type = 'l',
                     col = col_center, lwd = 3)
    
    graphics::lines( time_grid, max_envelope_limit, lty = 1,
                     col = col_fence_structure, lwd = 3 )
    graphics::lines( time_grid, min_envelope_limit, lty = 1,
                     col = col_fence_structure, lwd = 3 )
    
    # Plotting vertical whiskers
    half.time_grid = which.min( abs( time_grid - 0.5 ) )
    graphics::lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
                     c( out$max_envelope_central[ half.time_grid ],
                        max_envelope_limit[ half.time_grid ] ),
                     lty = 1, col = col_fence_structure, lwd = 3 )
    
    graphics::lines( c( time_grid[ half.time_grid ], time_grid[ half.time_grid ] ),
                     c( out$min_envelope_central[ half.time_grid ],
                        min_envelope_limit[ half.time_grid ] ),
                     lty = 1, col = col_fence_structure, lwd = 3 )
    
    # Plotting outlying data
    if( length( ID_out ) > 0 )
    {
      #graphics::matplot( time_grid, t( roahd::toRowMatrixForm( Data$coefficients[ ID_out, ] ) ),
      #                   lty = 1, type = 'l', col = col_outlying, lwd = 3, add = T )
      graphics::matplot( time_grid, t( roahd::toRowMatrixForm( expand(Data)[ ID_out, ] ) ),
                         lty = 1, type = 'l', col = col_outlying, lwd = 2, add = T )
    }
  }
  
  return( list( Depth = Depths,
                Fvalue = Fvalue,
                ID_outliers = ID_out,
                Cov_estimator = Cov_estimator) )

}


eigencouples_estimation = function(Data, Cov, L, N){
  
  if( is.null(Cov) ){
    # Kendall's Tau function
    
    t <- Data$geometry$basis$timeStr
    Lt <- rep(list(t$grid),N)
    Ly <- split(t(expand(Data)), rep(1:ncol(t(expand(Data))), 
                                     each = nrow(t(expand(Data)))))
    interval <- c(t$t0, t$tP)
    
    # I don't actually need this...
    basis <- fda::create.bspline.basis(interval, nbasis = 13, norder = 4)
    
    phi <- KFPCA::KFPCA(Lt, Ly, interval, nK = L, bw = 1, 
                          bwK = 0.045, bwmean = 0.03, dataType = "Dense",
                   nRegGrid = length(t$grid), fdParobj = basis, more = FALSE)$FPC_dis
    
  } else{
    
    Eigen.Cov = eigen.Cov_operator( Cov, nelements = L  )
    phi <- Eigen.Cov$efunctions
    
  }
  
  p <- matrix(NA, nrow = N, ncol = L)
  q <- rep(0, L)
  rho <- rep(0, L)
  # TODO: check that doing in the other way works fine
  #phi <-  matrix(NA, nrow = L, ncol = L)
  
  # for each eigenfunction:
  for( i in 1:L ){
    
    #phi[,i] <- Eigen.Cov$efunctions[,i]
    
    if( is.null(Cov) ){
      # Kendall's Tau function
      values = expand(Data)
    }
    else{
      values = Data$coefficients
    }
    
    for( j in 1:N ){
      # compute projection over eigen-function i
      # TODO: CHECK THIS!!!!
      
      p[j,i] <- geometry::dot(values[j,], phi[,i])
    }
    
    # compute q_i as robust estimate of the variance of projected data
    q[i] <- robustbase::s_Qn(p[,i])
    
    if( i==1 ){
      rho[1] <- 1
    }
    else{
      rho[i] <- q[i] / q[1]
    }
    
  }
  
  return(list(evalues = rho, efunctions = phi))
}



KL_data_generation = function(Data, ecouples, L, N, trial_size){
  
  rho = ecouples$evalues
  phi = ecouples$efunctions
  
  P = dim(expand(Data))[2]
  
  # TODO: check if the dimension is P or L (reduced dimensionality)
  Data_gauss <-  matrix(0, nrow = trial_size, ncol = dim(phi)[1])
  
  # sample trial_size observations 
  for( k in 1: trial_size ){
    for( i in 1:L ){
      Data_gauss[k,] <- Data_gauss[k,] + rho[i]*rnorm(1,0,1)*phi[,i]
    }
  }
  
  # add geometry to generated data
  # TODO: check that I can do this!!
  Data_gauss = as.fData(Data_gauss, Data$geometry)
  # moreover check that expanding the data also for the case of Kendall's Tau
  # is not a problem (I generate data directly expanded)
  
  return(Data_gauss)
}






