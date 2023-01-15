
source('basisOp/functional_basis.R')

# Implementing a method dispatch for inner_product
# based on the class of timeStr
L2_inner_product = function( timeStr, x, y )
{
  UseMethod( 'L2_inner_product' )
}

################################################
# Provide a RCpp implementation of this!!
# ################################################
L2_inner_product.EvenTimeStr = function( timeStr, x, y )
{
  if( length( x ) != length( y ) || timeStr$P != length( x ) )
    stop( ' Error in L2_inner_product.EvenTimeStr: arguments sizes are not compliant.')

    z = x * y;

  return( sum( 0.5 * ( z[ -1 ] + z[ - timeStr$P ] ) * timeStr$h ) )
}

################################################
# Provide also a RCpp implementation of this!!
# ################################################
L2_inner_product.UnevenTimeStr = function( timeStr, x, y )
{
  if( length( x ) != length( y ) || timeStr$P != length( x )  )
    stop( ' Error in L2_inner_product.UnevenTimeStr: arguments sizes are not compliant.')

    z = x * y;

  return( as.numeric( 0.5 * ( z[ -1 ] + z[ - timeStr$P ] ) %*% timeStr$h ) )
}

L2_norm = function( timeStr, x )
{
  UseMethod( 'L2_norm' )
}

################################################
# Provide a RCpp implementation of this!!
# ################################################
L2_norm.EvenTimeStr = function( timeStr, x )
{
  if( timeStr$P != length( x ) )
    stop( ' Error in L2_norm.EvenTimeStr: arguments sizes are not compliant.')

    z = x * x;

  return( sqrt( sum( 0.5 * ( z[ -1 ] + z[ - timeStr$P ] ) * timeStr$h ) ) )
}

################################################
# Provide also RCpp implementation of this!!
# ################################################
L2_norm.UnevenTimeStr = function( timeStr, x )
{
  if( timeStr$P != length( x )  )
    stop( ' Error in L2_norm.UnevenTimeStr: arguments sizes are not compliant.')

  z = x * x;

  return( sqrt( as.numeric( 0.5 * ( z[ -1 ] + z[ - timeStr$P ] ) %*% timeStr$h ) ) )
}


# Implementing a method dispatch for computeMassMatrix
# based on the class of basis
computeMassMatrix = function( basis )
{
  UseMethod( 'computeMassMatrix' )
}

# Specialised method for Fourier class
computeMassMatrix.fBasisFourier = function( fourier.basis )
{
  return( sparseMatrix( i = seq( 1, fourier.basis$nelements ),
                        j = seq( 1, fourier.basis$nelements ),
                        x = rep( 1, fourier.basis$nelements ),
                        symmetric = TRUE ) )
}

# Specialised method for Bspline class
computeMassMatrix.fBasisBspline = function( bspline.basis )
{
  order = attr( bspline.basis$elements, 'degree' ) + 1

  i_list = j_list = v_list = NULL;

  for( iElement in 1 : bspline.basis$nelements )
  {
    # Since B-spline are compactly supported, I can compute only a subset of the
    # required inner products.
    v_temp = sapply( iElement : min( c( iElement + order, bspline.basis$nelements ) ),
                     function( jElement )( L2_inner_product( bspline.basis$timeStr,
                                                             bspline.basis$elements[ iElement, ],
                                                             bspline.basis$elements[ jElement, ] ) ) )

    for( jElement in seq_along( v_temp ) )
    {
      if( abs( v_temp[ jElement ] ) > .Machine$double.eps )
      {
        v_list = c( v_list, v_temp[ jElement ] )
        i_list = c( i_list, iElement )
        j_list = c( j_list, iElement + jElement - 1 )
      }
    }
  }

  return( sparseMatrix( i = i_list,
                        j = j_list,
                        x = v_list,
                        dims = rep( bspline.basis$nelements, 2 ),
                        symmetric = TRUE ) )
}

# Specialised method for custom basis
computeMassMatrix.fBasisCustom = function( custom.basis )
{
#   if( orthonormal )
#   {
#     return( sparseMatrix( i = rep( 1, custom.basis$nelements ),
#                           j = rep( 1, custom.basis$nelements ),
#                           x = rep( 1, custom.basis$nelements ),
#                           symmetric = TRUE ) )
#   } else {

    i_list = j_list = v_list = NULL

    for( iElement in 1 : custom.basis$nelements )
    {
      if( iElement > 1 )
      {
        v_temp = sapply( 1 : ( iElement - 1 ),
                         function( jElement )( L2_inner_product( custom.basis$timeStr,
                                                                 custom.basis$elements[ iElement, ],
                                                                 custom.basis$elements[ jElement, ] ) ) )

        for( jElement in seq_along( v_temp ) )
        {
          if( abs( v_temp[ jElement ] ) > .Machine$double.eps )
          {
            v_list = c( v_list, v_temp[ jElement ] )
            i_list = c( i_list, iElement )
            j_list = c( j_list, jElement )
          }
        }
      }

      v_list = c( v_list, L2_norm( custom.basis$timeStr,
                                   custom.basis$elements[ iElement, ] )^2 )
      i_list = c( i_list, iElement )
      j_list = c( j_list, iElement )
    }

    return( sparseMatrix( i = i_list,
                          j = j_list,
                          x = v_list,
                          dims = rep( custom.basis$nelements, 2 ),
                          symmetric = TRUE ) )
  # }
}

Geometry = function( basis, space )
{
  UseMethod( 'Geometry' )
}

Geometry.fBasisFourier = function( fourier.basis, space )
{
  M = NULL;

  switch( space,
          'L2' = {
            M = computeMassMatrix( fourier.basis )
          },
          'sH1' = {
            M = sparseMatrix( i = seq( 1, fourier.basis$nelements ),
                              j = seq( 1, fourier.basis$nelements ),
                              x = as.vector( computeFourierDerivativeCoefficients( fourier.basis, 1)^2 ),
                              symmetric = TRUE )

          },
          'H1' = {
            M = computeMassMatrix( fourier.basis )
            M = M + sparseMatrix( i = seq( 1, fourier.basis$nelements ),
                                  j = seq( 1, fourier.basis$nelements ),
                                  x = as.vector( computeFourierDerivativeCoefficients( fourier.basis, 1)^2 ),
                                  symmetric = TRUE )

          },
          'sH2' = {
            M = sparseMatrix( i = seq( 1, fourier.basis$nelements ),
                              j = seq( 1, fourier.basis$nelements ),
                              x = as.vector( computeFourierDerivativeCoefficients( fourier.basis, order = 2 )^2 ),
                              symmetric = TRUE )
          },
          'H2' = {
            M = computeMassMatrix( fourier.basis )
            M = M + sparseMatrix( i = seq( 1, fourier.basis$nelements ),
                                  j = seq( 1, fourier.basis$nelements ),
                                  x = as.vector( computeFourierDerivativeCoefficients( fourier.basis, 1) ),
                                  symmetric = TRUE )
            M = M + sparseMatrix( i = seq( 1, fourier.basis$nelements ),
                                  j = seq( 1, fourier.basis$nelements ),
                                  x = as.vector( computeFourierDerivativeCoefficients( fourier.basis, order = 2 )^2 ),
                                  symmetric = TRUE )
          },
          stop("Error in Geometry.fBasisFourier: at present, I support only L2, (s)H1, (s)H2!")
          )



  geometry = structure( list( Mass_matrix = M,
                              basis = fourier.basis,
                              space = space ),
                        class = c( paste( 'Fourier_', space, sep = '' ), space )
                      )
 return( geometry )
}

Geometry.fBasisBspline = function( bspline.basis, space )
{
  M = NULL;

  switch( space,
          'L2' = {
            M = computeMassMatrix( bspline.basis )
          },
          'sH1' = {

            if( attr( bspline.basis$elements, 'degree' ) == 0 )
              stop( 'Erorr: incompatible bspline degree and number of derivatives')

            prev.bspline.basis = generateBasis( type = 'Bspline', timeStr = bspline.basis$timeStr,
                                                N = bspline.basis$nelements - 1,
                                                degree = attr( bspline.basis$elements, 'degree') - 1 )
            M = computeMassMatrix( prev.bspline.basis )
            derCoeffs = computeBsplineDerivativeCoefficients( bspline.basis$timeStr,
                                                              attr( bspline.basis$elements, 'inner.knots' ),
                                                              attr( bspline.basis$elements, 'degree' ) )

            M = derCoeffs %*% M %*% t( derCoeffs )
          },
          'H1' = {

            if( attr( bspline.basis$elements, 'degree' ) == 0 )
              stop( 'Erorr: incompatible bspline degree and number of derivatives.')

            M = computeMassMatrix( bspline.basis )
            prev.bspline.basis = generateBasis( type = 'Bspline', timeStr = bspline.basis$timeStr,
                                                N = bspline.basis$nelements - 1,
                                                degree = attr( bspline.basis$elements, 'degree') - 1 )
            derCoeffs = computeBsplineDerivativeCoefficients( bspline.basis$timeStr,
                                                              attr( bspline.basis$elements, 'inner.knots' ),
                                                              attr( bspline.basis$elements, 'degree' ) )
            M = M + derCoeffs %*% computeMassMatrix( prev.bspline.basis ) %*% t( derCoeffs )
          },
          'sH2' = {

            if( attr( bspline.basis$elements, 'degree' ) <= 1 )
              stop( 'Erorr: incompatible bspline degree and number of derivatives')

            prev2.bspline.basis = generateBasis( type = 'Bspline', timeStr = bspline.basis$timeStr,
                                                 N = bspline.basis$nelements - 2,
                                                 degree = attr( bspline.basis$elements, 'degree') - 2 )
            M = computeMassMatrix( prev2.bspline.basis )
            derCoeffs.1 = computeBsplineDerivativeCoefficients( bspline.basis$timeStr,
                                                                attr( bspline.basis$elements, 'inner.knots' ),
                                                                attr( bspline.basis$elements, 'degree' ) )
            derCoeffs.2 = computeBsplineDerivativeCoefficients( bspline.basis$timeStr,
                                                                attr( bspline.basis$elements, 'inner.knots' ),
                                                                attr( bspline.basis$elements, 'degree' ) - 1 )
            M = derCoeffs.1 %*% derCoeffs.2 %*% M %*% t( derCoeffs.2 ) %*% t ( derCoeffs.1 )

          },
          'H2' = {

            if( attr( bspline.basis$elements, 'degree' ) <= 1 )
              stop( 'Erorr: incompatible bspline degree and number of derivatives')

            M = computeMassMatrix( bspline.basis )

            prev.bspline.basis = generateBasis( type = 'Bspline', timeStr = bspline.basis$timeStr,
                                                N = bspline.basis$nelements - 1,
                                                degree = attr( bspline.basis$elements, 'degree') - 1 )
            derCoeffs.1 = computeBsplineDerivativeCoefficients( bspline.basis$timeStr,
                                                                attr( bspline.basis$elements, 'inner.knots' ),
                                                                attr( bspline.basis$elements, 'degree' ) )

            M = M + derCoeffs.1 %*% computeMassMatrix( prev.bspline.basis ) %*% t( derCoeffs.1 )

            rm( list =  c( 'prev.bspline.basis' ) )

            prev2.bspline.basis = generateBasis( type = 'Bspline', timeStr = bspline.basis$timeStr,
                                                N = bspline.basis$nelements - 2,
                                                degree = attr( bspline.basis$elements, 'degree') - 2 )

            derCoeffs.2 = computeBsplineDerivativeCoefficients( bspline.basis$timeStr,
                                                                attr( bspline.basis$elements, 'inner.knots' ),
                                                                attr( bspline.basis$elements, 'degree' ) - 1)

            M = M + derCoeffs.1 %*% derCoeffs.2 %*% computeMassMatrix( prev2.bspline.basis ) %*% t( derCoeffs.2 ) %*% t( derCoeffs.1 )

          },
          stop("Error in Geometry.fBasisBspline: at present, I support only L2, (s)H1, (s)H2!")
  )



  geometry = structure( list( Mass_matrix = M,
                              basis = bspline.basis,
                              space = space ),
                        class = c( paste( 'Bspline_', space, sep = '' ), space )
  )
  return( geometry )
}

Geometry.fBasisCustom = function( custom_basis, space )
{
  M = NULL;

  switch( space,
          'L2' = {
            M = computeMassMatrix( custom_basis )
          },
          'sH1' = {
            custom_basis_der = generateBasis( 'Custom', timeStr = custom_basis$timeStr,
                                              N = custom_basis$nelements,
                                              elements = t( apply( custom_basis$elements, 1,
                                                                   function( x ) ( fd_1_II( custom_basis$timeStr,
                                                                                            custom_basis$elements ) ) ) ) )
             M = computeMassMatrix( custom_basis_der )
          },
          'H1' = {
            M = computeMassMatrix( custom_basis )

            custom_basis_der = generateBasis( 'Custom',
                                              timeStr = custom_basis$timeStr,
                                              N = custom_basis$nelements,
                                              elements = t( apply( custom_basis$elements, 1,
                                                                   function( x ) ( fd_1_II( custom_basis$timeStr, x ) ) ) ) )
            M = M + computeMassMatrix( custom_basis_der )
          },
          'sH2' = {
            M = computeMassMatrix( custom_basis )

            custom_basis_der = generateBasis( 'Custom',
                                              timeStr = custom_basis$timeStr,
                                              N = custom_basis$nelements,
                                              elements = t( apply( custom_basis$elements, 1,
                                                                   function( x ) ( fd_2_II( custom_basis$timeStr, x ) ) ) ) )
            M = M + computeMassMatrix( custom_basis_der )

          },
          'H2' = {
            M = computeMassMatrix( custom_basis )

            custom_basis_der_1 = generateBasis( 'Custom',
                                                timeStr = custom_basis$timeStr,
                                                N = custom_basis$nelements,
                                                elements = t( apply( custom_basis$elements, 1,
                                                                     function( x ) ( fd_1_II( custom_basis$timeStr, x ) ) ) ) )
            M = M + computeMassMatrix( custom_basis_der_1 )

            rm( list = c( 'custom_basis_der_1' ) )

            custom_basis_der_2 = generateBasis( 'Custom',
                                                timeStr = custom_basis$timeStr,
                                                N = custom_basis$nelements,
                                                elements = t( apply( custom_basis$elements, 1,
                                                                     function( x ) ( fd_2_II( custom_basis$timeStr, x ) ) ) ) )
            M = M + computeMassMatrix( custom_basis_der_2 )

          },
          stop("Error in Geometry.fBasisCustom: at present, I support only L2, (s)H1, (s)H2!")
        )

  geometry = structure( list( Mass_matrix = M,
                              basis = custom_basis,
                              space = space ),
                        class = c( paste( 'Custom_', space, sep = '' ), space )
  )
  return( geometry )
}

computeFourierDerivativeCoefficients = function( basis, order = 1  )
{
  all( order > 0 ) || stop ( ' Error in computeFourierDerivativeCoefficients: orders must be positive')

  L =  timeStr$tP - timeStr$t0

  indexes = attr( basis$elements, 'indexes' )

  coeffs = array( 0, dim = c( length( order ), basis$nelements ) )

  for( i in seq_along( order ) )
  {
    for( j in seq_along( indexes ) )
    {
      if( indexes[ j ] %% 2 == 0 ){
        if( indexes[ j ] != 0 )
        {
          coeffs[ i, j ] = ( - 1 )^floor( ( order[ i ] + 1 ) / 2 ) *
                           ( 2 * pi * ceiling( indexes[ j ] / 2 ) /  L )^( order[ i ] )
        } else {
          coeffs[ i, j ] = 0
        }
      } else {
        coeffs[ i, j ] = ( - 1 )^floor( order[ i ] / 2 ) *
                         ( 2 * pi * ceiling( indexes[ j ] / 2 ) / L )^( order[ i ] )
      }
    }
  }

  return( coeffs )
}


computeBsplineDerivativeCoefficients = function( timeStr, inner.breaks, degree )
{
  require( Matrix )

  if( any( sapply( inner.breaks, function(x)( x >= timeStr$tP || x <= timeStr$t0 ) ) ) )
  {
    stop( 'Error in computeBsplineDerivativeCoefficients: knots provided are not inner')
  }

  if( degree == 0 )
  {
    stop( 'Warning in computeBsplineDerivativeCoefficients: you are differentiating a constant bspline basis')
    return(  Matrix( 0, ncol = length( inner.breaks ) + 1,
                        nrow = length( inner.breaks ) + 1,
                        sparse = TRUE ) )
  }

  order = degree + 1 # order
  N = length( inner.breaks ) + degree + 1
  N.intervals = length( inner.breaks ) + 1

  knots.primitive = c( rep( timeStr$t0, order ),
                       inner.breaks,
                       rep( timeStr$tP, order ) )
  offset.primitive =  ( 1 - ( 2 - order ) )
  offset.derivative = ( 1 - ( 3 - order ) )


  coeffs = Matrix( 0, nrow = length( inner.breaks ) + degree + 1,
                      ncol = length( inner.breaks ) + degree,
                      sparse = TRUE )

  for( i in ( 2 - order ) : N.intervals )
  {
    if( i >= ( 3 - order ) )
    {
      coeffs[ offset.primitive + i, offset.derivative + i ] = ( order - 1 ) /
                                                              ( knots.primitive[ offset.primitive + i + order - 1 ] -
                                                                knots.primitive[ offset.primitive + i ] )
    }

    if( i + 1 <= N.intervals )
    {
      coeffs[ offset.primitive + i, offset.derivative + i + 1 ] = - ( order - 1 ) /
                                                                    ( knots.primitive[ offset.primitive + i + order ] -
                                                                      knots.primitive[ offset.primitive + i + 1 ] )
    }
  }

  return( coeffs )
}

inner_product = function( x, y, geometry = NULL, env = NULL )
{
  # From a small simulation study it seems the passage of environment does
  # not entail a greater efficiency in the use of inner_product, either in
  # terms of memory management or computational time.
  # Think about removing it and unifying the call to a passage-by-copy of
  # the geometry object or even a C++ implementation
  if( is.null( env ) )
  {
    if( is.null( geometry ) )
         stop( 'please, provide a geometry or an environment to inner_product' )
    else
        env = as.environment( geometry )
  } else if( ! is.null( geometry ) )

      stop( 'You cannot provide both geometry and environment to inner_product' )

  # Now x is a column matrix of length( x ) rows
  if( is.vector( x ) ) x = t( as.matrix( x ) )

  # Now y is a column matrix of length( x ) rows
  if( is.vector( y ) ) y = t( as.matrix( y ) )


  if( nrow( x ) != nrow( y ) )
  {
    nrow( x ) == 1 || nrow( y ) == 1 || stop('Mismatching arguments in inner_products')
  }

  # Here I exploit cicrular indexing and sapply over a vector of indexes in order to
  # perform the pairwise inner_products. The - 1 in the vector of indexes is for
  # technical reasons having to do with the circular indexing via %% operator.
  # `get' makes a clever use of the environment built upon the geometry object
  # or given as a function parameter
  res = sapply(  0 : ( max( c( nrow( x ), nrow( y ) ) ) - 1 ),
                function( i )( as.numeric( x[ ( i %% nrow( x ) + 1 ), ] %*%
                                                 get( 'Mass_matrix', env ) %*%
                                                 y[ ( i %% nrow( y ) + 1 ), ] ) ) )

  return( res )
}

norm = function( x, geometry = NULL, env = NULL )
{
  # From a small simulation study it seems the passage of environment does
  # not entail a greater efficiency in the use of norm, either in terms of
  # memory management or computational time.
  # Think about removing it and unifying the call to a passage-by-copy of
  # the geometry object or even a C++ implementation
  if( is.null( env ) )
  {
    if( is.null( geometry ) )
      stop( 'please, provide a geometry or an environment to norm' )
    else
      env = as.environment( geometry )
  } else if( ! is.null( geometry ) )

    stop( 'You cannot provide both geometry and environment to norm' )

  # Now x is a column matrix of length( x ) rows
  if( is.vector( x ) ) x = t( as.matrix( x ) )

  # Here I exploit cicrular indexing and sapply over a vector of indexes in order to
  # perform the pairwise norm computations. The - 1 in the vector of indexes is for
  # technical reasons having to do with the circular indexing via %% operator.
  # `get' makes a clever use of the environment built upon the geometry object
  # or given as a function parameter
  #res = sapply(  0 : ( nrow( x ) - 1 ),
  #               function( i )( sqrt( as.numeric( x[ ( i %% nrow( x ) + 1 ), ] %*%
  #                                                  get( 'Mass_matrix', env ) %*%
  #                                                  x[ ( i %% nrow( x ) + 1 ), ] ) ) ) )

  # NB: this implementation is non-sense. It is equivalent to:
  
  # this gives the norm of every sample (row in the matrix)
  res = sapply(  1 : ( nrow( x ) ),
                  function( i )( sqrt( as.numeric( x[ i, ] %*%
                                                     get( 'Mass_matrix', env ) %*%
                                                     x[ i, ] ) ) ) )
  
  
  return( res )
}

frobenius_norm = function( x, geometry = NULL, env = NULL ){
  
  row_norms = norm(x, geometry, env)
  
  res = sqrt(sum( norm(x, geometry = geometry)^2 ))
  
  return(res)
}

normalise = function( x, geometry )
{
  # Now x is a column matrix of length( x ) rows
  if( is.vector( x ) ) x = t( as.matrix( x ) )

  return( x / norm( x, geometry ) )
}

fd_1_II = function( timeStr, data )
{
  UseMethod( 'fd_1_II' )
}

fd_2_II = function( timeStr, data )
{
  UseMethod( 'fd_2_II' )
}

fd_1_II.EvenTimeStr = function( timeStr, data )
{
  if( is.null( dim( data ) ) ) {
    data = t( as.matrix( data ) )
  }

  der = matrix( 0, ncol = ncol( data ), nrow = nrow( data ) )

  # x_{ i } inner
  der[ , 2 : ( timeStr$P - 1 ) ] = ( data[ , 3 : timeStr$P ] - data[ , 1 : ( timeStr$P - 2 ) ] ) / ( 2 * timeStr$h )

  # x_{ 0 }
  der[ , 1 ] = ( - 3 / 2 * data[ , 1 ] + 2 * data[ , 2 ] - 1 / 2 * data[ , 3 ] ) / timeStr$h

  # x_{ P }
  der[ , timeStr$P ] = ( 3 / 2 * data[ , timeStr$P ] - 2 * data[ , timeStr$P - 1 ] + 1 / 2 * data[ , timeStr$P - 2 ] ) / timeStr$h


  return( der )
}

fd_2_II.EvenTimeStr = function( timeStr, data )
{
  if( is.null( dim( data ) ) ) {
    data = t( as.matrix( data ) )
  }

  der = matrix( 0, ncol = ncol( data ), nrow = nrow( data ) )

  # x_{ i } inner
  der[ , 2 : ( timeStr$P - 1 ) ] = ( data[ , 3 : timeStr$P ] - 2 * data[ , 2 : ( timeStr$P - 1 ) ] + data[ , 1 : ( timeStr$P - 2 ) ] ) / timeStr$h^2

  # x_{ 0 }
  der[ , 1 ] = ( 2 * data[ , 1 ] - 5 * data[ , 2 ] + 4 * data[ , 3 ] - data[ , 4 ] ) / timeStr$h^2

  # x_{ P }
  der[ , timeStr$P ] = ( 2 * data[ , timeStr$P ] - 5 * data[ , timeStr$P - 1 ] + 4 * data[ , timeStr$P - 2 ] - data[ , timeStr$P - 3 ] ) / timeStr$h^2

  return( der )
}


fd_1_II.UnevenTimeStr = function( timeStr, data )
{
  stop( ' THIS IS A MESS!! Implement it!!')
}

fd_2_II.UnevenTimeStr = function( timeStr, data )
{
  stop( ' THIS IS A MESS!! Implement it!!')
}
