
covSph = function( fData2, median = NULL, ... )
{
  if( is.null( median ) )
  {
    median = median.fData2( fData2, ... )
  }

  temp = normalise( (fData2 - median)$coefficients, geometry = fData2$geometry )

  C = Matrix( cov( temp ) )

  return( structure( list( coefficients = C,
                          geometry = fData2$geometry ),
                    class = c( 'Spherical_cov', 'Cov_operator' ) ) )
}

covMed = function( fData2, alpha = 3/4, c = 2, median = NULL, ... )
{
  
  # first check the values of hyperparameters
  testit::assert('Hyperparameter c should be positive', {c>0})
  testit::assert('Hyperparameter alpha should be in (0.5,1)', {alpha > 0.5 & alpha <1})
  
  # compute median if not given, using ASG algorithm
  if( is.null( median ) )
  {
    median = median.fData2( fData2, ... )
  }
  
  # start ASG for Median Covariation matrix estimation
  
  L = dim(fData2$coefficients)[2]
  M = matrix(0, nrow = L, ncol = L)
  M_bar = matrix(0, nrow = L, ncol = L)
  
  # for each samples (row)
  for( i in 1:fData2$N){
    
    gamma = c / (max( c(i-1, 1) ))^alpha
    
    X_hat = (fData2$coefficients[i,] - median$coefficients)
    Y = t(X_hat) %*% X_hat
    
    # Stochastic Gradient step
    
    #M = M + gamma * (Y-M)/frobenius_norm(Y-M, geometry = fData2$geometry)
    
    # using non negative modification of the algorithm
    M = M + (Y-M) * min(1.0, gamma /frobenius_norm(Y-M, geometry = fData2$geometry))
    
    # averaging past values
    M_bar = M_bar - (M_bar -M)/i
    
  }
  
  return( structure( list( coefficients = M_bar,
                           geometry = fData2$geometry ),
                     class = c( 'Median_covariation', 'Cov_operator' ) ) )
}


median.fData2 = function( fD, ASG = FALSE, c = 2, alpha = 3/4 )
{
  if( ASG )
  {
    m = rep( 0, fD$geometry$basis$nelements )

    m_ave = m

    for( iObs in 1 : nrow( fD$coefficients ) )
    {
      gamma = c / max( c( iObs - 1, 1) )^alpha

      xmm = fD$coefficients[ iObs, ] - m

      m = m + gamma * xmm / norm( xmm, geometry = fD$geometry )

      m_ave = m_ave - 1 / ( iObs ) * ( m_ave - m )
    }
    return ( as.fData2( m_ave, fD$geometry ) )
  } else {
    stop( ' * * * Finish me! ')
  }
}
