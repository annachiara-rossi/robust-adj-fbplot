
rm( list = ls() )

message( ' * * * Fix me! ' )
setwd( '~/University/Phd/research/FDA_utilities/rfda/R/test')

source( '../basisOp/functional_basis.R')
source( '../geometry/geometry.R')
source( '../fdata/fData.R')
source( '../stats/stats.R' )

# FOURIER BASIS -----------------------------------------------------------

require( fda )

timeStr = EvenTimeStr( seq( 0, 8, 1 / (1e2 - 1 ) ) )

my.basis = generateFourierElements( timeStr, 5, id.subset = c( 0, 1, 2, 5, 7 ))

quartz()
matplot( timeStr$grid, t( my.basis ), type = 'l', lty = 1  )


fda.basis = create.fourier.basis( range( timeStr$grid ), nbasis = 9 )

quartz()
plot( fda.basis )




# BSPLINE BASIS -----------------------------------------------------------

require(fda)

T0 = 0
T1 = 1

timeStr = EvenTimeStr( seq(0,1, length.out = 1e2 ) )

degree.list = 1 : 3

N.list = seq( 8, 20)

failed = NULL;

id.failed = data.frame( degree.id = NULL, N.id = NULL );

error = array( 0, dim = c( length( degree.list ), length( N.list ) ) )

for( i in seq_along( degree.list ) )
{
  for( j in seq_along( N.list ) )
  {
    breaks = seq( T0, T1, length.out = N.list[[ j ]] - degree.list[[ i ]] + 1 )

    breaks.inner = breaks[ - c( 1, length(breaks) ) ]

    bspline.basis.fd = create.bspline.basis( c( T0, T1), norder = degree.list[[ i ]] + 1,
                                          nbasis = N.list[[ j ]], breaks = breaks )

    bspline.basis = eval.basis( timeStr$grid, bspline.basis.fd )

    basis = t( generateBsplineBasis( timeStr, N = N.list[[ j ]], degree = degree.list[[ i ]],
                                     inner.breaks = breaks.inner ) )

    error[ i, j ] = max( abs ( basis - bspline.basis ) )

  }
}


# GENERATE BASIS OBJECT ---------------------------------------------------

T0 = 0
T1 = 1

timeStr = EvenTimeStr( seq(0,1, length.out = 1e2 ) )

# Fourier
N = 3

fourier.basis = generateBasis( 'Fourier', timeStr, N, id.subset = 'sine' )

quartz()
plot( fourier.basis )

fourier.basis = generateBasis( 'Fourier', timeStr, N, id.subset = 'cosine' )

quartz()
plot( fourier.basis )

fourier.basis = generateBasis( 'Fourier', timeStr, N )

quartz()
plot( fourier.basis )

# Bspline
N = 10

bspline.basis = generateBasis( 'Bspline', timeStr, N, degree = 4)

quartz()
plot( bspline.basis )


# UNEVEN TIME GRID AND BASIS ----------------------------------------------

uTimeStr = UnevenTimeStr( ( exp( seq( 0, 1, length.out = 1e3 ) ) - 1 ) / ( exp(1) - 1 )  )

plot( uTimeStr )

# fourier
N = 3

fourier.basis = generateBasis( 'Fourier', uTimeStr, N, id.subset = 'sine' )

quartz()
plot( fourier.basis )

# bspline

N = 20
degree = 3

bspline.basis = generateBasis( 'Bspline', uTimeStr, N = N, degree = degree )

quartz()
plot( bspline.basis )






