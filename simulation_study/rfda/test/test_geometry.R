
rm( list = ls() )

message( ' * * * Fix me! ' )
setwd( '~/University/Phd/research/FDA_utilities/rfda/R/test')


source( '../basisOp/functional_basis.R')
source( '../geometry/geometry.R')
source( '../fdata/fData.R')

# INNER PRODUCTS ----------------------------------------------------------

timeStr = EvenTimeStr( seq( 0, 1, length.out = 1e2 ) )
UtimeStr = UnevenTimeStr( seq( 0, 1, length.out = 1e2 ) )

x = rep( 1, 1e2 )
y = rep( 2, 1e2 )

stopifnot( abs( L2_inner_product( timeStr, x, y ) - 2 ) <= .Machine$double.eps )
stopifnot( abs( L2_inner_product( UtimeStr, x, y ) - 2 ) <= .Machine$double.eps )

stopifnot( abs( L2_norm( timeStr, x ) - sqrt( L2_inner_product( timeStr, x, x ) ) ) <= .Machine$double.eps )
stopifnot( abs( L2_norm( timeStr, y ) - sqrt( L2_inner_product( timeStr, y, y ) ) ) <= .Machine$double.eps )


# MASS MATRIX -------------------------------------------------------------

N = 20

timeStr = EvenTimeStr( seq( 0, 1, length.out = 1e3 ) )

# Fourier
fourier.basis = generateBasis( 'Fourier', timeStr, N )

plot( fourier.basis )

tic = proc.time()
invisible( computeMassMatrix( fourier.basis ) )
toc = proc.time()
print( paste( ' Elapsed time ', (toc - tic)[3] ) )

# Bspline
bspline.basis = generateBasis( 'Bspline', timeStr, N )

plot( bspline.basis )

tic = proc.time()
M  = computeMassMatrix( bspline.basis )
toc = proc.time()
print( paste( ' Elapsed time ', (toc - tic)[3] ) )


# DIFFERENTIATE BSPLINE BASIS ---------------------------------------------

timeStr = EvenTimeStr( seq(0, 1, length.out = 1e3 ) )

degree = 3
inner_breaks = seq( 0, 1, length.out = 11 )[ -c( 1, 11 ) ]
N = length( inner_breaks ) + degree + 1

bspl.basis = generateBasis( 'Bspline', timeStr, N, degree = degree, inner.breaks = inner_breaks )
der.bspl.basis = generateBasis( 'Bspline', timeStr, N - 1, degree = degree - 1, inner.breaks = inner_breaks )

quartz()
par( mfrow = c( 1, 2 ) )
plot( bspl.basis )
plot( der.bspl.basis )

# Direct computation (Finite Differences FD1 )
basis.der.direct = array( 0, dim = dim( bspl.basis$elements ) )

for( i in 1 : nrow( basis.der.direct ) )
{
  if( i == 1 || i == 2 )
  {
    basis.der.direct[ i, ] = c( diff( bspl.basis$elements[ i, ] ) / timeStr$h, 0 )
  } else {
    basis.der.direct[ i, ] = c( 0, diff( bspl.basis$elements[ i, ] ) / timeStr$h )
  }
}

# Automatic computation
der.coefficients = computeBsplineDerivativeCoefficients( timeStr, inner_breaks, degree )

basis.der = array( 0, dim = c( bspl.basis$nelements, timeStr$P ) )

for( i in seq( 1, bspl.basis$nelements ) )
{
  basis.der[ i, ] = as.numeric( t( t( der.bspl.basis$elements ) %*% der.coefficients[ i, ] ) )
}

quartz()
par( mfrow = c( 1, 3 ) )
matplot( timeStr$grid, t( basis.der.direct ), lty = 1, type = 'l' )
matplot( timeStr$grid, t( basis.der ), lty = 1, type = 'l' )
matplot( timeStr$grid, t( basis.der ), lty = 1, type = 'l', lwd = 2, col = 'darkblue' )
matplot( timeStr$grid, t( basis.der.direct ), lty = 2, type = 'l', col = 'red', lwd = 2, add = T )


# MASS MATRIX IN HILBERT SPACES -------------------------------------------

# FOURIER

timeStr = EvenTimeStr( seq( 0, 1, length.out = 1e3 ) )

fourier.basis = generateBasis( 'Fourier', timeStr, N = 5 )

geometry.L2 = Geometry( fourier.basis, space = 'L2' )

geometry.sH1 = Geometry( fourier.basis, space = 'sH1' )

geometry.H1 = Geometry( fourier.basis, space = 'H1' )

geometry.sH2 = Geometry( fourier.basis, space = 'sH2' )

geometry.H2 = Geometry( fourier.basis, space = 'H2' )


x = matrix( rep( 1, 5 ), ncol = 5, nrow = 5 )
y =  matrix( seq(1, 5), ncol = 5, nrow = 5 )

stopifnot( all( abs( inner_product( x, y, env = as.environment( geometry.L2 ) ) - 5 * ( 1 : 5 ) ) <= .Machine$double.eps ) )

stopifnot( all( abs( inner_product( x, y, geometry.L2 ) -  5 * ( 1 : 5 ) ) <= .Machine$double.eps ) )

stopifnot( all( abs ( norm( x, as.environment( geometry.L2 ) ) - rep( sqrt( 5 ), 5 ) ) <= .Machine$double.eps ) )

stopifnot( all( abs ( norm( x, geometry.L2 ) - rep( sqrt( 5 ), 5 ) ) <= .Machine$double.eps ) )

M.L2 = matrix( 0, ncol = fourier.basis$nelements, nrow = fourier.basis$nelements )
M.sH1 = matrix( 0, ncol = fourier.basis$nelements, nrow = fourier.basis$nelements )
M.sH2 = matrix( 0, ncol = fourier.basis$nelements, nrow = fourier.basis$nelements )

for( i in 1 : nrow( M.L2 ) )
{

  vec.L2 = sapply( 0 : ( i - 1 ),
                   function( j ) ( L2_inner_product( timeStr, fourier.basis$elements[ i, ], fourier.basis$elements[ j + 1, ] ) ) )

  vec.L2 = sapply( vec.L2, function( x ) ( ifelse( abs( x ) >= 1e-14, x, 0 )  ) )

  vec.sH1 = sapply( 0 : ( i - 1 ),
                    function( j ) ( L2_inner_product( timeStr,
                                                      fd_1_II( timeStr, fourier.basis$elements[ i, ] ),
                                                      fd_1_II( timeStr, fourier.basis$elements[ j + 1, ] ) ) ) )

  vec.sH1 = sapply( vec.sH1, function( x ) ( ifelse( abs( x ) >= 1e-14, x, 0 )  ) )

  vec.sH2 = sapply( 0 : ( i - 1 ),
                    function( j ) ( L2_inner_product( timeStr,
                                                      fd_2_II( timeStr, fourier.basis$elements[ i, ] ),
                                                      fd_2_II( timeStr, fourier.basis$elements[ j + 1, ] ) ) ) )

  vec.sH2 = sapply( vec.sH2, function( x ) ( ifelse( abs( x ) >= 1e-14, x, 0 )  ) )


  M.L2[ i, 1 : i ] = vec.L2;
  M.L2[ 1 : i, i  ] = vec.L2;

  M.sH1[ i, 1 : i ] = vec.sH1;
  M.sH1[ 1 : i, i  ] = vec.sH1;

  M.sH2[ i, 1 : i ] = vec.sH2;
  M.sH2[ 1 : i, i  ] = vec.sH2;

}

Matrix::norm( Matrix( M.L2, sparse = TRUE )  -  geometry.L2$Mass_matrix , type = 'f' ) / Matrix::norm( geometry.L2$Mass_matrix )

Matrix::norm( Matrix( M.sH1, sparse = TRUE ) - geometry.sH1$Mass_matrix, type = 'f' ) / Matrix::norm( geometry.sH1$Mass_matrix )

Matrix::norm( abs( Matrix( M.sH1 + M.L2, sparse = TRUE ) - geometry.H1$Mass_matrix ), type = 'f' ) / Matrix::norm( geometry.H1$Mass_matrix )

Matrix::norm( Matrix( M.sH2, sparse = TRUE ) - geometry.sH2$Mass_matrix ) / Matrix::norm( geometry.sH2$Mass_matrix )

Matrix::norm( Matrix( M.sH2 + M.sH1 + M.L2, sparse = TRUE ) - geometry.H2$Mass_matrix ) / Matrix::norm( geometry.H2$Mass_matrix )

# BSPLINE

timeStr = EvenTimeStr( seq( 0, 1, length.out = 1e3 ) )

bspline.basis = generateBasis( 'Bspline', timeStr, N = 50, degree = 3 )

geometry.L2 = Geometry( bspline.basis, space = 'L2' )

geometry.sH1 = Geometry( bspline.basis, space = 'sH1' )

geometry.H1 = Geometry( bspline.basis, space = 'H1' )

geometry.sH2 = Geometry( bspline.basis, space = 'sH2' )

geometry.H2 = Geometry( bspline.basis, space = 'H2' )


x = matrix( rep( 1, 5 ), ncol = 50, nrow = 30 )
y =  matrix( seq(1, 30), ncol = 50, nrow = 30 )


stopifnot( all( abs( inner_product( x, y, env = as.environment( geometry.L2 ) ) - ( 1 : 30 )  ) <= 1e-14 ) )

stopifnot( all( abs( inner_product( x, y, geometry.L2 ) - ( 1 : 30 )  ) <= 1e-14 ) )

stopifnot( all( abs( norm( x, as.environment( geometry.L2 ) ) - rep( 1, 30 ) ) <= 1e-15 ) )

stopifnot( all( abs( norm( x, geometry.L2 ) - rep( 1, 30 )  ) <= 1e-15 ) )


M.L2 = matrix( 0, ncol = bspline.basis$nelements, nrow = bspline.basis$nelements )
M.sH1 = matrix( 0, ncol = bspline.basis$nelements, nrow = bspline.basis$nelements )
M.sH2 = matrix( 0, ncol = bspline.basis$nelements, nrow = bspline.basis$nelements )

for( i in 1 : nrow( M.L2 ) )
{

  vec.L2 = sapply( 0 : ( i - 1 ),
                   function( j ) ( L2_inner_product( timeStr, bspline.basis$elements[ i, ], bspline.basis$elements[ j + 1, ] ) ) )

  vec.L2 = sapply( vec.L2, function( x ) ( ifelse( abs( x ) >= 1e-14, x, 0 )  ) )

  vec.sH1 = sapply( 0 : ( i - 1 ),
                    function( j ) ( L2_inner_product( timeStr,
                                                      fd_1_II( timeStr, bspline.basis$elements[ i, ] ),
                                                      fd_1_II( timeStr, bspline.basis$elements[ j + 1, ] ) ) ) )

  vec.sH1 = sapply( vec.sH1, function( x ) ( ifelse( abs( x ) >= 1e-14, x, 0 )  ) )

  vec.sH2 = sapply( 0 : ( i - 1 ),
                    function( j ) ( L2_inner_product( timeStr,
                                                      fd_2_II( timeStr, bspline.basis$elements[ i, ] ),
                                                      fd_2_II( timeStr, bspline.basis$elements[ j + 1, ] ) ) ) )

  vec.sH2 = sapply( vec.sH2, function( x ) ( ifelse( abs( x ) >= 1e-14, x, 0 )  ) )


  M.L2[ i, 1 : i ] = vec.L2;
  M.L2[ 1 : i, i  ] = vec.L2;

  M.sH1[ i, 1 : i ] = vec.sH1;
  M.sH1[ 1 : i, i  ] = vec.sH1;

  M.sH2[ i, 1 : i ] = vec.sH2;
  M.sH2[ 1 : i, i  ] = vec.sH2;

}

Matrix::norm( Matrix( M.L2, sparse = TRUE )  -  geometry.L2$Mass_matrix , type = 'f' ) / Matrix::norm( geometry.L2$Mass_matrix )

Matrix::norm( Matrix( M.sH1, sparse = TRUE ) - geometry.sH1$Mass_matrix, type = 'f' ) / Matrix::norm( geometry.sH1$Mass_matrix )

Matrix::norm( abs( Matrix( M.sH1 + M.L2, sparse = TRUE ) - geometry.H1$Mass_matrix ), type = 'f' ) / Matrix::norm( geometry.H1$Mass_matrix )

Matrix::norm( Matrix( M.sH2, sparse = TRUE ) - geometry.sH2$Mass_matrix ) / Matrix::norm( geometry.sH2$Mass_matrix )

Matrix::norm( Matrix( M.sH2 + M.sH1 + M.L2, sparse = TRUE ) - geometry.H2$Mass_matrix ) / Matrix::norm( geometry.H2$Mass_matrix )

# EFFICIENCY STUDY

bspline.basis = generateBasis( 'Bspline', timeStr, N = 1e2, degree = 3 )

geometry.L2 = Geometry( bspline.basis, space = 'L2' )

x = matrix( rep( 1, bspline.basis$nelements ), ncol = bspline.basis$nelements, nrow = 1e4 )
y = matrix( seq(1, bspline.basis$nelements), ncol = bspline.basis$nelements, nrow = 1e4 )

env = as.environment( geometry.L2 )

gc( reset = TRUE )

tic = proc.time()
temp = inner_product( x, y, env = env  )
toc = proc.time()
print( paste( ' Elapsed time ', (toc - tic)[3] ) )

gc()

gc( reset = TRUE )

tic = proc.time()
temp = inner_product( x, y, geometry.L2 )
toc = proc.time()
print( paste( ' Elapsed time ', (toc - tic)[3] ) )

gc()


# FINITE DIFFERENCES OF ORDER 2 (1st and 2nd) -----------------------------

timeStr = EvenTimeStr( seq( 0, 1, length.out =  1e3 ) )

fourier.basis = generateBasis( 'Fourier', timeStr, N = 5 )

quartz()
par( mfrow = c( 1, 3 ) )
plot( fourier.basis )
matplot( timeStr$grid, t( fd_1_II( timeStr, fourier.basis$elements ) ), type = 'l', lty = 1 )
matplot( timeStr$grid[ -1 ], apply( fourier.basis$elements, 1, function( x )( diff( x ) / timeStr$h ) ), lty = 1, type = 'l' )

quartz()
par( mfrow = c( 1, 3 ) )
plot( fourier.basis )
matplot( timeStr$grid, t( fd_2_II( timeStr, fourier.basis$elements ) ), type = 'l', lty = 1 )
matplot( timeStr$grid[ -c( 1, timeStr$P ) ], apply( fourier.basis$elements, 1, function( x )( diff( diff( x ) ) / timeStr$h^2 ) ), lty = 1, type = 'l' )


# FDATA -------------------------------------------------------------------

timeStr = EvenTimeStr( seq( 0, 1, length.out =  1e3 ) )

# fourier
fourier.basis = generateBasis( 'Fourier', timeStr, N = 5, id.subset = c( 2, 3, 4, 5, 6 ) )

geometry.L2 = Geometry( fourier.basis, space = 'L2' )
geometry.sH1 = Geometry( fourier.basis, space = 'sH1' )
geometry.H1 = Geometry( fourier.basis, space = 'H1' )
geometry.sH2 = Geometry( fourier.basis, space = 'sH2' )
geometry.H2 = Geometry( fourier.basis, space = 'H2' )

fD.L2  = fData( fourier.basis$elements, geometry = geometry.L2 )
fD.sH1 = fData( fourier.basis$elements, geometry = geometry.sH1 )
fD.H1  = fData( fourier.basis$elements, geometry = geometry.H1 )
fD.sH2 = fData( fourier.basis$elements, geometry = geometry.sH2 )
fD.H2  = fData( fourier.basis$elements, geometry = geometry.H2 )

Matrix::norm( Matrix( fD.L2$coefficients - diag( 1, fourier.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, fourier.basis$nelements ) )
Matrix::norm( Matrix( fD.sH1$coefficients - diag( 1, fourier.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, fourier.basis$nelements ) )
Matrix::norm( Matrix( fD.H1$coefficients - diag( 1, fourier.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, fourier.basis$nelements ) )
Matrix::norm( Matrix( fD.sH2$coefficients - diag( 1, fourier.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, fourier.basis$nelements ) )
Matrix::norm( Matrix( fD.H2$coefficients - diag( 1, fourier.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, fourier.basis$nelements ) )

quartz()
plot( fourier.basis, col = 'blue', lwd = 2 )
plot( fD.L2, add = TRUE, col = 'red', lwd = 2, lty = 2 )

quartz()
plot( fourier.basis, col = 'blue', lwd = 2 )
plot( fD.sH1, add = TRUE, col = 'red', lwd = 2, lty = 2 )

quartz()
plot( fourier.basis, col = 'blue', lwd = 2 )
plot( fD.H1, add = TRUE, col = 'red', lwd = 2, lty = 2 )

quartz()
plot( fourier.basis, col = 'blue', lwd = 2 )
plot( fD.sH2, add = TRUE, col = 'red', lwd = 2, lty = 2 )

quartz()
plot( fourier.basis, col = 'blue', lwd = 2 )
plot( fD.H2, add = TRUE, col = 'red', lwd = 2, lty = 2 )

# bspline
bspline.basis = generateBasis( 'Bspline', timeStr, N = 10 )

geometry.L2  = Geometry( bspline.basis, space = 'L2' )
geometry.sH1 = Geometry( bspline.basis, space = 'sH1' )
geometry.H1  = Geometry( bspline.basis, space = 'H1' )
geometry.sH2 = Geometry( bspline.basis, space = 'sH2' )
geometry.H2  = Geometry( bspline.basis, space = 'H2' )

fD.L2  = fData( bspline.basis$elements, geometry = geometry.L2 )
# fD.sH1 = fData( bspline.basis$elements, geometry = geometry.sH1 )
fD.H1  = fData( bspline.basis$elements, geometry = geometry.H1 )
# fD.sH2 = fData( bspline.basis$elements, geometry = geometry.sH2 )
fD.H2  = fData( bspline.basis$elements, geometry = geometry.H2 )

Matrix::norm( Matrix( fD.L2$coefficients - diag( 1, bspline.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, bspline.basis$nelements ) )
# Matrix::norm( Matrix( fD.sH1$coefficients - diag( 1, bspline.basis$nelements ) ) ) /
#               Matrix::norm( diag( 1, bspline.basis$nelements ) )
Matrix::norm( Matrix( fD.H1$coefficients - diag( 1, bspline.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, bspline.basis$nelements ) )
# Matrix::norm( Matrix( fD.sH2$coefficients - diag( 1, bspline.basis$nelements ) ) ) /
#               Matrix::norm( diag( 1, bspline.basis$nelements ) )
Matrix::norm( Matrix( fD.H2$coefficients - diag( 1, bspline.basis$nelements ) ) ) /
  Matrix::norm( diag( 1, bspline.basis$nelements ) )

quartz()
plot( bspline.basis, col = 'blue', lwd = 2 )
plot( fD.L2, add = TRUE, col = 'red', lwd = 2, lty = 2 )

# quartz()
# plot( bspline.basis, col = 'blue', lwd = 2 )
# plot( fD.sH1, add = TRUE, col = 'red', lwd = 2, lty = 2 )

quartz()
plot( bspline.basis, col = 'blue', lwd = 2 )
plot( fD.H1, add = TRUE, col = 'red', lwd = 2, lty = 2 )

# quartz()
# plot( bspline.basis, col = 'blue', lwd = 2 )
# plot( fD.sH2, add = TRUE, col = 'red', lwd = 2, lty = 2 )

quartz()
plot( bspline.basis, col = 'blue', lwd = 2 )
plot( fD.H2, add = TRUE, col = 'red', lwd = 2, lty = 2 )
