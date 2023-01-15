
# Common plot method for the class "TimeStr"
plot.TimeStr = function( timeStr, ylab, ... )
{
  if( missing( ylab ) ) ylab = 'time'

  plot( 1 : timeStr$P, timeStr$grid, xlab = 'time point ID', ylab, ... )
}

# Constructor of the class "EvenTimeStr"
EvenTimeStr = function( grid, h = NULL, P = NULL )
{
  ! missing( grid ) || stop( ' Error in EvenTimeStr: missing grid in EvenTimeStr' )

  message( ' * * * Handle better the global constant here ' )
  all( abs( diff( unique( diff( grid ) ) ) ) < 1e-14 ) ||
    stop( ' Error in EvenTimeStr: uneven grid passed in EvenTimeStr')

  if( is.null( h ) || is.null( P ) )
  {
    h = grid[ 2 ] - grid[ 1 ]
    P = length( grid )
  } else if ( P != length( grid ) ||
              h < 0 ||
              abs ( h - ( grid[ 2 ] - grid[ 1 ] ) ) > .Machine$double.eps ){

      stop( ' Error in EvenTimeStr: arguments passed to EvenTimeStr not compliant')
  }

  structure( list( grid = grid,
                   h = h,
                   P = P,
                   t0 = grid[ 1 ],
                   tP = grid[ P ]),
             class = c( 'EvenTimeStr', 'TimeStr' ) )
}

# Constructor of the class "UnevenTimeStr"
UnevenTimeStr = function( grid, h = NULL, P = NULL)
{
  ! missing( grid ) || stop( ' Error in UnevenTimeStr: missing grid in UnevenTimeStr' )

  if( is.null( h ) || is.null( P ) )
  {
    h = diff( grid )
    P = length( grid )
  } else if ( P != length( grid ) ||
              any ( h < 0 ) ||
              any ( diff( grid ) - h  > .Machine$double.eps ) ) {

    stop( ' Error in UnevenTimeStr: arguments passed to UnevenTimeStr not compliant')
  }

  structure( list( grid = grid,
                   h = h,
                   P = P,
                   t0 = grid[ 1 ],
                   tP = grid[ P ] ),
             class = c( 'UnevenTimeStr', 'TimeStr' ) )
}

# Function to generate Fourier elements on a specified grid
generateFourierElements = function( timeStr, N, id.subset = NULL )
{
  L =  timeStr$tP - timeStr$t0

  L > 0 || stop( ' Error in generateFourierBasis: a >= b ' )

  basis = array( 0, dim = c( N, timeStr$P ) )

  if( is.null( id.subset ) ) {
    id.subset = ( 1 : N ) - 1
  } else if ( ! is.numeric( id.subset ) ){
    if( id.subset == 'sine' ){
        id.subset = 2 * ( 1 : N ) - 1
    } else if ( id.subset == 'cosine' ){
        id.subset = 2 * ( 1 : N ) - 2
    }
  } else {
      ( length( id.subset ) == N ) ||
       stop( " Error in generateFourierBasis: mismatching arguments provided")

      all( sort( id.subset, decreasing = FALSE ) == id.subset ) ||
       stop( ' Error in generateFourierBasis: you must provide a sorted subset of ids ' )
  }

  for( i in seq_along( id.subset ) )
  {
    if( id.subset[ i ] %% 2 == 0 ){
      if( id.subset[ i ] != 0 )
      {
        basis[ i, ] = sqrt( 2 / L ) * cos( 2 * pi * ceiling( id.subset[ i ] / 2 ) * ( timeStr$grid - timeStr$t0 ) / L )
      } else {
        basis[ i, ] = 1 / sqrt( L )
      }
    } else {
      basis[ i, ] = sqrt( 2 / L ) * sin( 2 * pi * ceiling( id.subset[ i ] / 2 ) * ( timeStr$grid - timeStr$t0 ) / L )
    }
  }

  return( structure( basis, indexes = id.subset ) )
}

# Function to generate a basis of splines
generateBsplineBasis = function( timeStr, N = NULL, degree = NULL, inner.breaks = NULL )
{
  library( splines )
  library( Matrix )

  if( is.null( degree ) ) degree = 3;

  if( is.null( inner.breaks ) ){

    if( is.null( N ) )
      stop( 'Error in generateBsplineBasis: you must provide (at least) either inner.breaks or N' )

    # Building a uniform grid of N - 1 - degree inner knots. For simplicity, I add and remove the boundary knots.
    # Here I assume the time structure can also have an unevenly spaced grid
    inner.breaks = quantile( timeStr$grid, prob = seq( 0, 1, length.out = ( N - 1 - degree ) + 2 ) )[ - c( 1,  N - 1 - degree + 2  ) ]
    names( inner.breaks )  = NULL

    bspline.basis = computeBasisSpline( timeStr, inner.breaks = inner.breaks, degree )
  } else {

    if(  min( inner.breaks ) <= timeStr$t0 || max( inner.breaks ) >= timeStr$tP )
      stop( 'Error in generateBsplineBasis: you have to provide a vector of inner break points, if any.')

    if( ! is.null( N ) && N != degree + 1 + length( inner.breaks ) )
      stop( 'Error in generateBsplineBasis: N must equal degree + 1 + length( inner.breaks )')

      bspline.basis = computeBasisSpline( timeStr, inner.breaks, degree )
    }

  return( bspline.basis )
}

# Generic basis function generator
generateBasis = function( type, timeStr, N, ... )
{

  if( missing( type ) )
    stop( 'Error in generateBasis: you have to specify the basis you want')

  if( missing( N ) )
    stop( 'Error in generateBasis: you have to specify the number of basis elements you want')

  switch ( type,
    'Fourier' = {
      basis = generateFourierElements( timeStr, N, ... )
    },
    'Bspline' = {
      basis = generateBsplineBasis( timeStr, N, ... )
    },
    'Custom' = {
      basis = generateCustomBasis( timeStr, N, ... )
    },
    stop("Error in generateBasis: please, enter a supported basis flag!")
  )

  return( structure( list( timeStr = timeStr,
                           elements = basis,
                           nelements = nrow( basis ),
                           type = type ),
                      class = c( paste( 'fBasis', type, sep = ''), 'fBasis' ) ) )

  ( basis.obj )
}

generateCustomBasis = function( timeStr, N, rule = NULL, elements = NULL )
{
  if( is.null( rule ) && is.null( elements ) )
  {
    stop( ' You must specify either the generating rule or the elements themselves.
            Generating rule is a function that accepts a grid like timeStr$grid and
            returns the matrix of basis elements evaluations ( N elements x timeStr$P ) ')
  }

  if( ! is.null( rule ) && ! is.null( elements ) )
  {
    stop( ' You cannot specify both a rule and the elements' )
  }

  if( ! is.null( rule ) )
  {
    return( rule( timeStr$grid ) )
  } else {
    if( ! is.matrix( elements ) ){
      stop( 'Error in generateCustomBasis: basis elements provided should be in matrix form.')
    } else if( timeStr$P != ncol( elements ) ){
      stop( 'Error in generateCustomBasis: basis and time structure provided mismatch.')
    }
    return( elements )
  }
}

computeBasisSpline = function( timeStr, inner.breaks, degree )
{
  if( any( sapply( inner.breaks, function(x)( x >= timeStr$tP || x <= timeStr$t0 ) ) ) )
  {
    stop( 'Error in computeBasisSpline: knots provided are not inner')
  }

  order = degree + 1

  N.intervals = length( inner.breaks ) + 1

  N = order + length( inner.breaks )

  knots.new = c( timeStr$t0,
                 inner.breaks,
                 timeStr$tP )

  bs.new = array( 0, dim = c( N.intervals, timeStr$P ) )

  for( i in 1 : N.intervals )
  {
    bs.new[ i, ] = 1 * ( sapply( timeStr$grid,
                                 function( x ) ( x >= knots.new[ i ] &&
                                                ifelse( i == nrow( bs.new ), TRUE, x < knots.new[ i + 1 ] ) ) ) )
  }

  if( order == 1 )
  {
    return ( bs.new )
  } else {

    for( K in 2 : order )
    {

      offset.old = 1 - ( 3 - K )
      offset.new = 1 - ( 2 - K )

      knots.old = knots.new
      knots.new = c( timeStr$t0,
                     knots.old,
                     timeStr$tP )

      bs.old = bs.new
      bs.new = array( 0, dim = c( N.intervals + K - 1, timeStr$P ) )

      for( i in ( 2 - K ) : N.intervals  )
      {
        # To understand this, make a "tree" diagram of B1,1 ... B4,1
        # and update up to the third layer
        if( i >= ( 3 - K )  )
        {
          bs.new[ offset.new + i, ] = ( timeStr$grid - knots.new[ offset.new + i ] )  /
                                      ( knots.new[ offset.new + i + K - 1 ] - knots.new[ offset.new + i ] ) *
                                        bs.old[ offset.old + i,  ]
        }

        # To understand this, make a "tree" diagram of B1,1 ... B4,1
        # and update up to the third layer
        if( i + 1 <= N.intervals )
        {

          bs.new[ offset.new + i, ] = bs.new[ offset.new + i, ] + ( knots.new[ offset.new + i + K ] - timeStr$grid ) /
                                                                   ( knots.new[ offset.new + i + K ] - knots.new[ offset.new + i + 1 ] ) *
                                                                     bs.old[ offset.old + i + 1, ]
        }
      }

      K = K + 1;
    }
  }

  # If it needs to be sparse, it is made sparse
  bs.new = Matrix( bs.new )

  attr( bs.new, 'knots' ) = knots.new
  attr( bs.new, 'inner.knots' ) = inner.breaks
  attr( bs.new, 'degree' ) = degree

  return( bs.new )
}


as.basis = function( ... )
{
  UseMethod( 'as.basis' )
}

as.basis.fPC = function( fPC )
{
  elements = expand( as.fData2( fPC$efunctions, fPC$geometry ) )

  return( generateBasis( 'Custom', fPC$geometry$basis$timeStr,
                                   fPC$nelements,
                                   elements = elements ) )
}


# Overload of plot method for fBasis data
plot.fBasis = function( fBasis, lty, xlab, ylab, ... )
{
  if( missing( lty ) ) lty = 1
  if( missing( xlab ) ) xlab = 'time'
  if( missing( ylab ) ) ylab = ''

  matplot( fBasis$timeStr$grid, t( fBasis$elements ), type = 'l', lty = lty, xlab = xlab, ylab = ylab, ... )
}
