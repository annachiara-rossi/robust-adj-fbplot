# NB: added geigen::geigen because geigen can also be taken from fda
source('fdata/fData.R')

mean.fData2 = function( fD )
{

  return( structure( list( coefficients = matrix( apply( fD$coefficients, 2, mean ), nrow = 1 ),
                           N = 1,
                           geometry = fD$geometry ),
                     class = c( 'fData2' ) ) )
}

cov.operator = function( fD )
{
  coefficients = Matrix( cov( fD$coefficients ),  )

  return( structure( list( coefficients = coefficients,
                          geometry = fD$geometry ),
                          class = 'Cov_operator' ) )
}

as.Cov.function = function( Cov_operator )
{
  return( expand( Cov_operator ) )
}


princomp.fData2 = function( fD, cov_operator = NULL, threshold = 1, components = NULL )
{
  if( is.null( cov_operator ) )
  {
    cov_operator = cov.operator( fD )
  } else if( ! identical( cov.operator$geometry, fD$geometry ) )
      stop( 'You provided fD and cov_operator with mismatching geometries' )

  require( geigen )


  message( ' FIX The use of sparse geigen solver!! ' )
  eigs = geigen::geigen( as.matrix( fD$geometry$Mass_matrix %*%
                              cov_operator$coefficients %*%
                                fD$geometry$Mass_matrix ),
                 as.matrix( fD$geometry$Mass_matrix ), symmetric = TRUE )

  # The eigen-decomposition provided by geigen is in increasing order, therefore since I usually work with decreasing eigen-elements, I have to revert them.
  L = fD$geometry$basis$nelements

  if( is.null( components ) )
  {
    id_subset = L -  which( ( cumsum( rev( eigs$values^2 ) ) / sum( eigs$values^2 ) ) <= threshold ) + 1

    if( length( id_subset ) == 0 ) stop( 'Cumulated variance is always higher than the specified threshold')
  } else {

    id_subset = seq( L, L - components + 1 )
  }

  return( structure( list( nelements = length( id_subset ),
                           evalues = eigs$values[ id_subset ],
                           efunctions = normalise( t( eigs$vectors[ , id_subset ] ), fD$geometry ),
                           geometry = fD$geometry ),
                     class = 'fPC' ) )
}

eigen.Cov_operator = function( cov_operator, nelements = NULL )
{
  require( geigen )

  message( ' FIX The use of sparse geigen solver!! ' )

  eigs = geigen::geigen( as.matrix( cov_operator$geometry$Mass_matrix %*%
                              cov_operator$coefficients %*%
                                cov_operator$geometry$Mass_matrix ),
                 as.matrix( cov_operator$geometry$Mass_matrix ), symmetric = TRUE )

  # The eigen-decomposition provided by geigen is in increasing order, therefore since I usually work with decreasing eigen-elements, I have to revert them.

  if( is.null( nelements ) ){
    nelements = nrow( cov_operator$coefficients )
  }

  indexes = seq( nrow( cov_operator$coefficients ),
                 nrow( cov_operator$coefficients ) - nelements + 1 )

  return( structure( list( nelements = length( indexes ),
                           evalues = eigs$values[ indexes ],
                           efunctions = normalise( t( eigs$vectors[ , indexes ] ),
                                                   cov_operator$geometry ),
                           geometry = cov_operator$geometry ),
                     class = 'fPC' ) )
}

as.fData2.fPC = function( fPC )
{
  return( as.fData2( coefficients = diag( sqrt( fPC$evalues ) ),
                    geometry = Geometry( as.basis( fPC ), fPC$geometry$space ) ) )
}

plot.fPC = function( fPC, ...  )
{
  plot( as.fData2( coefficients = fPC$efunctions,
                                 geometry = fPC$geometry ), ... )
}

expand.Cov_operator = function( Cov_operator )
{
  fValues = Matrix( Cov_operator$coefficients[ 1, 1 ] *
                    outer( Cov_operator$geometry$basis$elements[ 1, ],
                           Cov_operator$geometry$basis$elements[ 1, ] ) )

  for( i in 2 : Cov_operator$geometry$basis$nelements )
  {
    for( j in 1 : ( i - 1 ) )
    {
      fValues = fValues + Cov_operator$coefficients[ i, j ] *
                          outer( Cov_operator$geometry$basis$elements[ i, ],
                                 Cov_operator$geometry$basis$elements[ j, ] )

      fValues = fValues + Cov_operator$coefficients[ j, i ] *
                          outer( Cov_operator$geometry$basis$elements[ j, ],
                                 Cov_operator$geometry$basis$elements[ i, ] )
    }

    fValues = fValues + Cov_operator$coefficients[ i, i ] *
                        outer( Cov_operator$geometry$basis$elements[ i, ],
                               Cov_operator$geometry$basis$elements[ i, ] )

  }

  return( structure( list( values = fValues,
                           geometry = Cov_operator$geometry ),
                     class = 'Cov_function' ) )
}


plot.Cov_operator = function( C, expand = FALSE, ... )
{
  if( expand )
  {
    plot( expand( C ), ...  )
  } else {
    image( C$coefficients,
           colorkey = TRUE,
           useRaster = TRUE,
           col.regions = heat.colors,
           ... )
  }
}

plot.Cov_function = function( C, ... )
{
  image( C$values,
         colorkey = TRUE,
         useRaster = TRUE,
         col.regions = heat.colors,
         ... )
}





