
source('geometry/geometry.R')


generateFD = function( N, distroRule, geometry )
{
  coefficients = t( sapply( 1 : N , function(...)( distroRule() ) ) )

  if( any( dim( coefficients ) != c( N, geometry$basis$nelements ) ) )
  {
    stop( 'Error in generateFD: arguments provided are mismatching.')
  }

  return( as.fData2( coefficients, geometry ) )
}



distroRule.Gaussian = function( mu, Sigma )
{
  library(mvtnorm)

  if( any( rep( length( mu ), 2 ) != dim( Sigma ) ) )
  {
    stop( 'Error in distroRule.Gaussian: arguments provided are mismatching.')
  }

  return( function() ( rmvnorm( 1,  mean = mu, sigma = Sigma ) ) )
}


