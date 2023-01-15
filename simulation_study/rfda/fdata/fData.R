
source('rfda/geometry/geometry.R')

fData = function( data, geometry = NULL, basis = NULL, space = NULL )
{
  if( is.null( geometry ) )
  {
    if( is.null( basis ) || is.null( space ) )
      stop( 'You have to provide at least geometry or basis and space desired')
    else
      geometry = Geometry( basis, space )
  } else
  {
    if( ! is.null( basis ) && ! identical( basis, geometry$basis ) ||
        ! is.null( space ) && ! identical( space, geometry$space ) )
      stop( 'Basis, space and geometry parameters are not compliant')

  }

  if( is.null( dim( data ) ) )
  {
    data = t( as.matrix( data ) )
  }

  coefficients = matrix( 0, nrow = nrow( data ), ncol = geometry$basis$nelements )

  for( i in 1 : nrow( data ) )
  {
    coefficients[ i, ] = as.numeric( solve( geometry$Mass_matrix,
                                            computeProjectionRHS( geometry, data[ i, ] ) ) )
  }

  return( structure( list( coefficients = coefficients,
                           N = nrow( coefficients ),
                           geometry = geometry ),
                     class = c( 'fData' ) ) )
}

as.fData = function( coefficients, geometry )
{
  if( is.null( dim( coefficients ) ) )
  {
    coefficients = t( as.matrix( coefficients ) )
  }

  return( structure( list( coefficients = coefficients,
                           N = nrow( coefficients ),
                           geometry = geometry ),
                     class = c( 'fData' ) ) )
}


expand = function( ... )
{
  UseMethod( 'expand' )
}

expand.fData = function( fD  )
{
  return( as.matrix( t( apply( fD$coefficients,
                               1,
                               function( x ) ( colSums( fD$geometry$basis$elements * as.vector( x ) ) )
                             )
                      )
                   )
        )
}

computeProjectionRHS = function( geometry, data_vector )
{
  UseMethod( 'computeProjectionRHS' )
}


computeProjectionRHS.L2 = function( geometry, data_vector )
{
    # TRY TO VECTORIZE THIS: IT DEPENDS ON THE POSSIBILITY TO VECTORIZE THE SOLUTION OF THE LINEAR PROBLEM
    as.matrix(
               apply( geometry$basis$elements, 1, function( x ) ( L2_inner_product( geometry$basis$timeStr, x, data_vector ) ) ),
               ncol = 1 )
}


computeProjectionRHS.sH1 = function( geometry, data_vector )
{
  # TRY TO VECTORIZE THIS: IT DEPENDS ON THE POSSIBILITY TO VECTORIZE THE SOLUTION OF THE LINEAR PROBLEM
  apply( geometry$basis$elements, 1, function( x ) ( L2_inner_product( geometry$basis$timeStr,
                                                                       fd_1_II( geometry$basis$timeStr, x ),
                                                                       fd_1_II( geometry$basis$timeStr, data_vector ) ) ) )
}

computeProjectionRHS.H1 = function( geometry, data_vector )
{
  # TRY TO VECTORIZE THIS: IT DEPENDS ON THE POSSIBILITY TO VECTORIZE THE SOLUTION OF THE LINEAR PROBLEM
  apply( geometry$basis$elements, 1, function( x ) ( L2_inner_product( geometry$basis$timeStr,
                                                                       fd_1_II( geometry$basis$timeStr, x ),
                                                                       fd_1_II( geometry$basis$timeStr, data_vector ) ) ) ) +
  computeProjectionRHS.L2( geometry, data_vector )
}

computeProjectionRHS.sH2 = function( geometry, data_vector )
{
  # TRY TO VECTORIZE THIS: IT DEPENDS ON THE POSSIBILITY TO VECTORIZE THE SOLUTION OF THE LINEAR PROBLEM
  apply( geometry$basis$elements, 1, function( x ) ( L2_inner_product( geometry$basis$timeStr,
                                                                       fd_2_II( geometry$basis$timeStr, x ),
                                                                       fd_2_II( geometry$basis$timeStr, data_vector ) ) ) )
}

computeProjectionRHS.H2 = function( geometry, data_vector )
{
  # TRY TO VECTORIZE THIS: IT DEPENDS ON THE POSSIBILITY TO VECTORIZE THE SOLUTION OF THE LINEAR PROBLEM
  apply( geometry$basis$elements, 1, function( x ) ( L2_inner_product( geometry$basis$timeStr,
                                                                          fd_2_II( geometry$basis$timeStr, x ),
                                                                          fd_2_II( geometry$basis$timeStr, data_vector ) ) ) ) +
  computeProjectionRHS.H1( geometry, data_vector )
}

plot.fData = function( fD, lty, ... )
{
    if( missing( lty ) ) lty = 1

      matplot( fD$geometry$basis$timeStr$grid, t( expand( fD ) ), type = 'l', lty = lty, ... )
}

"+.fData" = function( fD, A )
{
  if( class( A ) == 'fData' )
  {

    if( nrow( A$coefficients ) == 1 )
    {
      fD$coefficients = t( t( fD$coefficients ) + as.vector( A$coefficients ) )
    } else {
      fD$coefficients = fD$coefficients + A$coefficients
    }

  } else if( is.null( dim( A ) ) || nrow( A ) == 1 ){

    A = as.vector( A )

    stopifnot( length( A ) %in% c( 1, fD$geometry$basis$timeStr$P ) )

    fD$coefficients = t( t( fD$coefficients ) +
                           as.vector( fData( rep( A,
                                                  fD$geometry$basis$timeStr$P / length( A ) ),
                                             fD$geometry )$coefficients ) )
  } else if( is.matrix( A ) ) {

    stopifnot( ncol( A ) == fD$geometry$basis$timeStr$P )
    stopifnot( nrow( A ) == nrow( fD$coefficients ) )

    fD$coefficients = fD$coefficients + fData( A, fD$geometry )$coefficients

  }

  return( fD )
}

"-.fData" = function( fD, A )
{
  if( class( A ) == 'fData' )
  {
    if( nrow( A$coefficients ) == 1 )
    {
      fD$coefficients = t( t( fD$coefficients ) - as.vector( A$coefficients ) )
    } else {
      fD$coefficients = fD$coefficients - A$coefficients
    }

  } else if( is.null( dim( A ) ) || nrow( A ) == 1 ){

    stopifnot( length( A ) %in% c( 1, fD$geometry$basis$timeStr$P ) )

    fD$coefficients = t( t( fD$coefficients ) -
                           as.vector( fData( rep( A,
                                                  fD$geometry$basis$timeStr$P / length( A ) ),
                                             fD$geometry )$coefficients ) )
  } else if( is.matrix( A ) ) {

    stopifnot( ncol( A ) == fD$geometry$basis$timeStr$P )
    stopifnot( nrow( A ) == nrow( fD$coefficients ) )

    fD$coefficients = fD$coefficients - fData( A, fD$geometry )$coefficients
  }

  return( fD )
}

"*.fData" = function( fD, a )
{
  stopifnot( 1 %in% dim( as.matrix( a ) ) )

  fD$coefficients = fD$coefficients * a

  return( fD )
}

"/.fData" = function( fD, a )
{
  stopifnot( 1 %in% dim( as.matrix( a ) ) )

  fD$coefficients = fD$coefficients / a

  return( fD )
}

