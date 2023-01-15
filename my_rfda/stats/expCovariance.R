######################################################################################
######################################################################################
######################################################################################
######################################################################################
########                                                                      ########
########        Utilities to compute the eigen-decomposition of a             ########
########        generic Ornstein-Uhlenbeck Gaussian process, by solving       ########
########        the related Fredholm integral equation of 2nd kind            ########
########        ( details in "Stochastic Finite Elements" by Ghanem and       ########
########          Spanos, Dover and "Spectral Methods for Uncertainty         ########
########          Quantification", by LeMaitre and Knio, Springer.)           ########
########                                                                      ########
######################################################################################
######################################################################################

ExpCovfunction = function( alpha, beta )
{
  message( ' * * * Think about inserting me into a MaternCovarianceFunction method')

  return( function( s, t )( alpha * exp( - beta * abs( s - t ) ) ) )
}

ExpCov_eigencouples = function( timeStr, K, expCovFunction, eps = pi/1e3, tol = .Machine$double.eps^0.5 )
{
  message( ' * * * allow K = 1 ')

  if( any( c( timeStr$t0, timeStr$tP ) != c( 0, 1 ) ) ){

    message( ' * * * This is still experimental: check me!')

    # Rescaling the time interval to [0,1], by absorbing (b-a) into beta
    beta_old = get( 'beta', environment( expCovFunction ) )
    assign( 'beta', beta_old * ( timeStr$tP - timeStr$t0 ), environment( expCovFunction ) )

    timeStr = EvenTimeStr( ( timeStr$grid - timeStr$t0 ) / ( timeStr$tP - timeStr$t0 )  )
  }

  alpha = get( 'alpha', environment( expCovFunction ) )

  beta = get( 'beta', environment( expCovFunction ) )

  K.even = floor( K / 2 ) + K %% 2

  K.odd = K - K.even

  # Finding omegas
  omega = rep( 0, K )

  omega[ 2 * seq( 1, K.even ) - 1 ] = find_omega_even( K.even, beta, eps, tol )

  omega[ 2 * seq( 1, K.odd ) ] = find_omega_odd( K.odd, beta, eps, tol )

  # Turning omegas into (sorted) eigenvalues
  lambda = compute_eigen_values( omega, alpha, beta )

  # Computing eigenfunctions
  eigen_functions = array( 0, dim = c( K, timeStr$P ) )

  eigen_functions[ 2 * seq( 1, K.even ) - 1, ] = eigen_functions_even( timeStr,
                                                                       omega[ 2 * seq( 1, K.even ) - 1 ] )

  eigen_functions[ 2 * seq( 1, K.odd ), ] = eigen_functions_odd( timeStr,
                                                                 omega[ 2 * seq( 1, K.odd ) ] )

  return( list( eigen_functions = eigen_functions,
                eigen_values = lambda ) )

}

# Find even omega solutions so that odd eigenfunctions can be computed
find_omega_even = function( N, beta, eps = pi/1e3, tol = .Machine$double.eps^0.5 )
{
  # Uniroot employs a bisection root finding algorithm
  f = function(w)( 1 - w * tan( w / 2 ) / beta )

  omega_even = rep(0,N);

  if( N > 1 )
  {
    # eps, pi + eps, 3 * pi + eps, 5 * pi + eps, ...
    lh_extrema = c( 0, seq(1, 2 * (N - 1) - 1, 2) ) * pi + eps

    # pi - eps, 3 * pi - eps, 5 * pi - eps, ...
    rh_extrema = seq(1, 2 * N - 1, 2) * pi - eps
  } else
  {
    lh_extrema = c( 0 + eps )
    rh_extrema = pi - eps
  }

  for( i in 1 : N )
  {
    omega_even[ i ] = uniroot( f, interval = c( lh_extrema[i], rh_extrema[i] ),
                               extendInt = 'downX', tol = tol )[["root"]]
  }

  return( omega_even )
}

# Find odd omega solutions so that odd eigenfunctions can be computed
find_omega_odd = function( N, beta, eps = pi/1e3, tol = .Machine$double.eps^0.5 )
{
  # Here the eps must be much less than before, since the zeros tend to move closer
  # and closer to the singularities as N grows

  # Uniroot employs a bisection root finding algorithm
  f = function(w)( w / beta + tan( w / 2 ) )

  omega_odd = rep(0,N);

  if( N > 1 )
  {
    # pi + eps, 3 * pi + eps, 5 * pi + eps, ...
    lh_extrema = seq(1, 2 * N - 1, 2) * pi + eps

    # 3 * pi - eps, 5 * pi - eps, 7 * pi - eps, ...
    rh_extrema = seq(3, 2 * N + 1, 2) * pi - eps
  } else
  {
    lh_extrema = pi + eps
    rh_extrema = 3 * pi - eps
  }

  for( i in 1 : N )
  {
    omega_odd[ i ] = uniroot( f, interval = c( lh_extrema[i], rh_extrema[i] ),
                              extendInt = 'upX', tol = tol )[["root"]]
  }

  return( omega_odd )
}

# Compute even eigen-functions
eigen_functions_even = function( timeStr, omega_even )
{
  N = length( omega_even )

  eigen_f_even = array( 0, dim = c( N, timeStr$P ) )

  for( i in 1 : N )
  {
    eigen_f_even[ i, ] = cos( omega_even[ i ] * ( timeStr$grid - 0.5 ) ) /
                         sqrt( 0.5 + sin( omega_even[ i ]) / ( 2 * omega_even[ i ] ) )
  }

  return( eigen_f_even )

}

# Compute odd eigen-functions
eigen_functions_odd = function( timeStr, omega_odd )
{
  N = length( omega_odd )

  eigen_f_odd = array( 0, dim = c( N, timeStr$P ) )

  for( i in 1 : N )
  {
    eigen_f_odd[ i, ] = sin( omega_odd[ i ] * ( timeStr$grid - 0.5 ) ) /
                        sqrt( 0.5 - sin( omega_odd[ i ] ) / ( 2 * omega_odd[ i ] ) )
  }

  return( eigen_f_odd )
}

# Compute eigen-values starting from omega, alpha and beta
compute_eigen_values = function( omega, alpha, beta )
{
  return( alpha * ( 2 / beta ) / ( 1 + ( omega / beta )^2 ) )
}


