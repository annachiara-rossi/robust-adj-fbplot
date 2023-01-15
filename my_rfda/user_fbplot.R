####________________________ User version of fbplot ________________________####
# in this version, an object of class fData2 from roahd can be passed
# which should have a $values field. 
# NB: THIS IS THE LOCAL VERSION

library(here)

#-------------------------------- dependencies --------------------------------#
# import Tarabelloni's code
source( 'basisOp/functional_basis.R') #my_rfda
source( 'geometry/geometry.R')
source( 'fdata/fData.R')
source( 'stats/stats.R')
source( 'stats/generateFD.R')
source( 'stats/expCovariance.R')
source( 'stats/robustEstimators.R' )

# import roahd depths computation
source(here('roahd/band_depths.R'))

# needed for KFPCA
require(interp)

# needed for kMRCD
source(here('kMRCD/kMRCD.R'))
#------------------------------------------------------------------------------#


#----------------------------- fences computation -----------------------------#
# Input: values of the fData2 object
# Output: envelope, fences and observations outside fences (outliers)

.my_fbplot_fData2 = function(Data, Depths = 'MBD', Fvalue = 1.5 )
{
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
#------------------------------------------------------------------------------#



#---------------------- functional boxplot implementation ---------------------#
# Input: 
#   - Data, an fData2 object
#   - adjust, can be a list e.g. list(VERBOSE = FALSE)

# added [mrcd_target, h, alpha, L] to be set by user when asking for adjustment
# and a specific covariance estimator.
# TODO: leave default version with the estimator we think is the best overall
# TODO: add an option to try all estimators
# TODO: give a list of the possible estimators

my_fbplot.fData2 = function( Data,
                         Depths = 'MBD',
                         Fvalue = 1.5,
                         adjust = FALSE,
                         display = TRUE,
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         ... )
{
  # Checking if depths have already been provided or must be computed
  if( is.character( Depths ) )
  {
    # Nice trick to encapsulate the information on the desired definition of
    # depth inside the vector that supposedly should contain depth values
    Depths_spec = Depths
    
    if( Depths_spec == 'MBD' )
    {
      Depths = MBD( Data$values, manage_ties = TRUE )
    } else {
      Depths = eval( parse( text = paste( Depths, '( Data$values )',
                                          sep = '' ) ) )
    }
  } else {
    stopifnot( length( Depths ) == Data$N )
  }
  
  if( ! is.list( adjust ) )
  {
    # Plain functional boxplot with default F value: F = 1.5
    out = .my_fbplot_fData2( Data$values, Depths, Fvalue )
    
  } else {
    
    nodenames = c( 'Cov_estimator', #'mrcd_target', 'alpha', 'L',
                   'N_trials', 'trial_size', 'TPR', 'F_min', 'F_max',
                   'tol', 'maxiter', 'VERBOSE' )
    unused = setdiff( names( adjust ), nodenames )
    
    # Checking for unused parameters
    if( length( unused ) > 0 )
    {
      for( i in unused )
        warning( 'Warning: unused parameter ', i, ' in adjust argument of fbplot' )
    }
    
    #--------------------------------------------------------------------------#
    Cov_estimator = ifelse( is.null( adjust$Cov_estimator ),
                            'OGK_Qn', # the originally used one
                            adjust$Cov_estimator )
    
    #--------------------------------------------------------------------------#
    
    
    # lists of possible estimators
    multivariate_est <- list("Ledoit_Wolf", "OGK_mad", "OGK_Qn", "MRCD", "kMRCD")
    functional_est <- list("Spherical", "Median", "Kendall")
    
    is_functional = Cov_estimator %in% functional_est
    
    if( !is_functional ){
      # function for multivariate case
      Fvalues = multivariate_case(Data, Depths_spec, Depths, Cov_estimator, adjust, ...)
      
    } else {
      # function for functional case
      Fvalues = functional_case(Data, Depths_spec, Cov_estimator, adjust, ...)
    }
    
    Fvalue = mean( Fvalues )
    
    out = .my_fbplot_fData2( Data$values, Depths, Fvalue = Fvalue  )
  }
  
  ID_out = out$ID_out
  
  #------------------------------- Plotting part ------------------------------#
  
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
    
    # TODO: this depends on the object
    time_grid = seq( Data$t0, Data$tP, length.out = Data$P )
    
    xlab = ifelse( is.null( xlab ), '', xlab )
    ylab = ifelse( is.null( ylab ), '', ylab )
    main = ifelse( is.null( main ), '', main )
    
    if( length( ID_out ) > 0 )
    {
      # Plotting non-outlying data
      graphics::matplot( time_grid,
                         t( Data$values[ - ID_out, ] ), lty = 1, type = 'l',
                         col = col_non_outlying,
                         ylim = range( Data$values ),
                         xlab = xlab, ylab = ylab, main = main, cex.main = 2, ... )
      
      # Computing maximum and minimum envelope
      max_envelope_limit = apply( Data$values[ - ID_out, ], 2, max )
      min_envelope_limit = apply( Data$values[ - ID_out, ], 2, min )
    } else {
      # Plotting all data
      graphics::matplot( time_grid,
                         t( Data$values ), lty = 1, type = 'l',
                         col = col_non_outlying,
                         ylim = range( Data$values ),
                         xlab = xlab, ylab = ylab, main = main, cex.main = 2, ... )
      
      # Computing maximum and minimum envelope
      max_envelope_limit = apply( Data$values, 2, max )
      min_envelope_limit = apply( Data$values, 2, min )
    }
    
    
    # Filling in the central envelope
    
    graphics::polygon( c(time_grid, rev( time_grid) ),
                       c( out$min_envelope_central, rev( out$max_envelope_central ) ),
                       col = col_envelope, border = NA)
    graphics::lines( time_grid, out$max_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    graphics::lines( time_grid, out$min_envelope_central, lty = 1, col = col_envelope, lwd = 3 )
    
    # Plotting the sample median
    graphics::lines( time_grid, Data$values[ which.max( Depths ), ], lty = 1, type = 'l',
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
      graphics::matplot( time_grid, t( roahd::toRowMatrixForm( Data$values[ ID_out, ] ) ),
                         lty = 1, type = 'l', col = col_outlying, lwd = 3, add = T )
    }
  }
  #------------------------------- end plotting -------------------------------#
  
  
  return( list( Depth = Depths,
                Fvalue = Fvalue,
                ID_outliers = ID_out) )   # had to remove the estimator in return because of the no-adj case
}


multivariate_case <- function(Data, Depths_spec, Depths, Cov_estimator, adjust, ...)
{
  
  # Estimation of robust covariance matrix
  kwargs = list(...)
  
  switch(Cov_estimator,
         # multivariate estimators
         
         Ledoit_Wolf = {
           Cov = cvCovEst::linearShrinkLWEst(Data$values)
         },
         OGK_mad = {
           Cov = robustbase::covOGK(Data$values, sigmamu = robustbase::s_mad )$cov
         },
         OGK_Qn = {
           Cov = robustbase::covOGK(Data$values, sigmamu = robustbase::s_Qn )$cov
         },
         MRCD = {
           mrcd_target = ifelse( is.null( kwargs$mrcd_target ),
                                 'equicorrelation', 
                                 kwargs$mrcd_target )
           alpha = ifelse( is.null( kwargs$alpha ), 0.75, kwargs$alpha )
           
           Cov = rrcov::CovMrcd(x = Data$values, target=mrcd_target, alpha = alpha)$cov
         },
         kMRCD = {
           alpha = ifelse( is.null( kwargs$alpha ), 0.75, kwargs$alpha )
           
           Cov = CovKmrcd(x = Data$values, alpha = alpha)$cov
         }
  )
  
  # Cholesky factor
  CholCov <- chol( Cov )
  
  # Centerline of the dataset
  centerline = Data$values[ which.max( Depths ), ]
  
  # CALL tuning function
  # N_trials, trial_size, TPR, F_min, F_max, tol, maxiter, VERBOSE
  return( tuning(is_functional = FALSE, Data, Depths_spec, adjust, 
                 centerline = centerline, CholCov = CholCov) )
}



functional_case <- function(Data, Depths_spec, Cov_estimator, adjust, ...)
{
  kwargs = list(...)

  # number of components in the basis, needed for functional operators
  L = ifelse( is.null( kwargs$L ), 10, kwargs$L )
  
  # define the geometry
  timeStr = EvenTimeStr( seq( Data$t0, Data$tP, length.out = Data$P ) )
  fourier.basis = generateBasis( 'Fourier', timeStr, N = L )
  geometry.L2 = Geometry( fourier.basis, space = 'L2' )
  
  # project on the basis
  Data_reduced = fData2(Data$values, geometry.L2)
  
  # Estimation of robust covariance matrix
  switch(Cov_estimator,
         
         # TODO: at this point I need data transformed wrt to a functional basis
         
         # functional operators
         Spherical = {
           Cov = covSph( Data_reduced, ASG = TRUE )   
         },
         Median = {
           Cov = covMed( Data_reduced, ASG = TRUE )   
         },
         Kendall = {
           # in this case I don't have an estimate of the covariance matrix,
           # but directly of the eigenfunctions
           Cov = NULL
         } 
  )
  
  N = nrow( Data_reduced$coefficients )
  
  ecouples = eigencouples_estimation(Data_reduced, Cov, L, N)
  # TODO: che tipo di dati mi servono in: eigencouples_estimation e l'altra mia funzione??
  
  # CALL tuning function
  return( tuning(is_functional = TRUE, Data, Depths_spec, adjust, 
                 geometry = geometry.L2, ecouples = ecouples, L = L, N = N) )
  
}


tuning <- function(is_functional, Data, Depths_spec, adjust, centerline = NULL, CholCov = NULL, 
                    geometry = NULL, ecouples = NULL, L = NULL, N = NULL, ...)
{

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
  
  Fvalues = rep( 0, N_trials )
  
  # The quantity to be minimized is the same
  cost_functional = function( F_curr )( length(
    .my_fbplot_fData2( Data_gauss,
                      Depths = Depths_spec,
                      Fvalue = F_curr )$ID_out ) /
      trial_size - TPR )
  
  for( iTrial in 1 : N_trials )
  {
    if( VERBOSE > 0 )
    {
      cat( ' * * * Iteration ', iTrial, ' / ', N_trials, '\n' )
    }
    
    if(!is_functional){
      Data_gauss = roahd::generate_gauss_fdata( trial_size, centerline, 
                                                CholCov = CholCov )
    } else {
      Data_gauss = expand(KL_data_generation(geometry, ecouples, L, N, trial_size))
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
  
  return(Fvalues)
  
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



KL_data_generation = function(geometry, ecouples, L, N, trial_size){
  
  # only need the geometry, not the original Data
  
  rho = ecouples$evalues
  phi = ecouples$efunctions
  
  P = geometry$basis$timeStr$P
  
  Data_gauss <-  matrix(0, nrow = trial_size, ncol = dim(phi)[1])
  
  # sample trial_size observations 
  for( k in 1: trial_size ){
    for( i in 1:L ){
      Data_gauss[k,] <- Data_gauss[k,] + rho[i]*rnorm(1,0,1)*phi[,i]
    }
  }
  
  # add geometry to generated data
  Data_gauss = as.fData2(Data_gauss, geometry)

  return(Data_gauss)
}










