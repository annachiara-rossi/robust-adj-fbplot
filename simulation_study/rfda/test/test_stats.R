
rm( list = ls() )

message( ' * * * Fix me! ' )
setwd( '~/University/Phd/research/FDA_utilities/rfda/R/test')

source( '../basisOp/functional_basis.R')
source( '../geometry/geometry.R')
source( '../fdata/fData.R')
source( '../stats/stats.R' )
source( '../stats/generateFD.R')
source( '../stats/expCovariance.R')

library(scales)

# STATISTICS --------------------------------------------------------------

require('fda')

timeStr = EvenTimeStr( seq(1, 365) )

# Fourier
fourier.basis = generateBasis( 'Fourier', timeStr, N = 25 )

# L2
geometry.L2 = Geometry( fourier.basis, space = 'L2' )
fD.L2 = fData( t( CanadianWeather$dailyAv[,,1] ), geometry.L2 )

quartz()
plot( fD.L2 )
plot( mean( fD.L2 ), add = TRUE, col = 'black', lwd = 2, lty = 1 )

# H1
geometry.H1 = Geometry( fourier.basis, space = 'H1' )
fD.H1 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.L2 )

# H2
geometry.H2 = Geometry( fourier.basis, space = 'H2' )
fD.H2 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.H2  )

quartz()
par( mfrow = c( 2, 2 ) )
plot( fD.L2 )
plot( fD.H1 )
plot( fD.H2 )

plot( fD.L2, col = 'blue' )
plot( mean( fD.L2 ), col = 'black', lty = 2, add = T )

plot( fD.H1, col = 'red', add = T )
plot( mean( fD.L2 ), col = 'black', lty = 2, add = T )

plot( fD.H2, col = 'forestgreen', add = T )
plot( mean( fD.H2 ), col = 'black', lty = 2, add = T )

quartz()
plot( cov.operator( fD.L2 ) )

quartz()
plot( cov.operator( fD.H1 ) )

quartz()
plot( cov.operator( fD.H2 ) )

fPC = princomp( fD.L2, components = 3 )
quartz()
plot( fPC, lty = 2 , col = c( 'darkred', 'darkblue', 'orange') )

# Bspline
bspline.basis = generateBasis( 'Bspline', timeStr, N = 10 )

# Bspline - L2
geometry.L2 = Geometry( bspline.basis, space = 'L2' )
fD.L2 = fData( t( CanadianWeather$dailyAv[,,1] ), geometry.L2 )

quartz()
plot( fD.L2 )
plot( mean( fD.L2 ), add = TRUE, col = 'black', lwd = 2, lty = 1 )

# Bspline - H1
geometry.H1 = Geometry( bspline.basis, space = 'H1' )
fD.H1 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.L2 )

# Bspline - H2
geometry.H2 = Geometry( bspline.basis, space = 'H2' )
fD.H2 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.H2  )

quartz()
par( mfrow = c( 2, 2 ) )
plot( fD.L2 )
plot( fD.H1 )
plot( fD.H2 )

plot( fD.L2, col = 'blue' )
plot( mean( fD.L2 ), col = 'black', lty = 2, add = T )

plot( fD.H1, col = 'red' )
plot( mean( fD.L2 ), col = 'black', lty = 2, add = T )
plot( fD.L2, col = 'yellow', add = T )


plot( fD.H2, col = 'forestgreen' )
plot( mean( fD.H2 ), col = 'black', lty = 2, add = T )

# PCA

fPC = princomp( fD.L2, components = 5 )
basis.fPC = as.basis( fPC )

quartz()
plot( basis.fPC, lty = 2 )

# L2 geometry for PCA
geometry.fPC.L2 = Geometry( basis.fPC, space = 'L2' )

fD.PCA.L2 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.fPC.L2 )

quartz()
par( mfrow = c( 1, 2 ) )
plot( fD.L2 )
plot( fD.PCA.L2 )

# H1 geometry for PCA
geometry.fPC.H1 = Geometry( basis.fPC, space = 'H1' )

fD.PCA.H1 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.fPC.H1 )

quartz()
par( mfrow = c( 1, 2 ) )
plot( fD.L2 )
plot( fD.PCA.H1 )

# H2 geometry for PCA
geometry.fPC.H2 = Geometry( basis.fPC, space = 'H2' )

fD.PCA.H2 = fData( t( CanadianWeather$dailyAv[ , , 1 ] ), geometry.fPC.H2 )

quartz()
par( mfrow = c( 1, 2 ) )
plot( fD.L2 )
plot( fD.PCA.H2 )


# ALGEBRA OF FDATA --------------------------------------------------------

timeStr = EvenTimeStr( seq( 0, 1, length.out = 1e3 ) )

bspl.basis = generateBasis( 'Bspline', timeStr, N = 20, degree = 3 )

geometry.L2 = Geometry( bspl.basis, space = 'L2' )

quartz()
plot( bspl.basis )

distroRule = distroRule.Gaussian( mu = rep( 0, bspl.basis$nelements),
                                  Sigma = diag( 1, bspl.basis$nelements ) )

N = 1e2
gauss_FD = generateFD( N, distroRule, geometry.L2 )

quartz()
plot( gauss_FD, col = hue_pal()( 10 )  )
plot( mean( gauss_FD ), add = TRUE )

quartz()
plot( mean( gauss_FD ) +
        matrix( sin( 2 * pi * timeStr$grid ), ncol = timeStr$P, nrow = 1 ), col = 'darkred' )
plot( mean( gauss_FD ), add = TRUE )

quartz()
plot( gauss_FD  + gauss_FD, col = 'darkred' )
plot( gauss_FD, col = 'darkblue', add = TRUE )

quartz()
plot( fData( sin( 2 * pi * timeStr$grid ), geometry.L2 ) - 1, col = 'black' )
abline( h = 0, lwd = 2, lty = 2, col = 'black' )

# GENERATION OF SYNTHETIC DATA --------------------------------------------


# GAUSSIAN

timeStr = EvenTimeStr( seq( 0, 1, length.out = 1e3 ) )

bspl.basis = generateBasis( 'Bspline', timeStr, N = 20, degree = 0 )

geometry.L2 = Geometry( bspl.basis, space = 'L2' )

quartz()
plot( bspl.basis )

distroRule = distroRule.Gaussian( mu = rep( 0, bspl.basis$nelements),
                                  Sigma = diag( 1, bspl.basis$nelements ) )

N = 1e2
gauss_FD = generateFD( N, distroRule, geometry.L2 )


quartz()
plot( gauss_FD )

# EXP COVARIANCE

timeStr = EvenTimeStr( seq( 0, 1, length.out = 2e2 ) )

# Amplitude factor
alpha = 0.12

# Correlation decay factor
beta = 0.4

expCov = ExpCovfunction( alpha, beta )

spectrum = ExpCov_eigencouples( timeStr, 10, expCov )


eigen.basis = generateBasis( 'Custom', timeStr, N = length( spectrum$eigen_values),
                             elements = spectrum$eigen_functions )

quartz()
par( mfrow = c(1,2) )
plot( eigen.basis )
matplot( timeStr$grid, t( spectrum$eigen_functions * sqrt( spectrum$eigen_values ) ),
         type = 'l', lty = 1, xlab = 'time', ylab = '' )

# Simulation of syntehtic Gaussian data on Exp cov basis

geometry.L2 = Geometry( eigen.basis, space = 'L2' )

expCovRule.Gauss = distroRule.Gaussian( mu = rep( 0, geometry.L2$basis$nelements ),
                                        Sigma = diag( spectrum$eigen_values ) )

fD.L2 = generateFD( 1e3, expCovRule.Gauss, geometry.L2 ) + sin( 4 * pi * timeStr$grid )

quartz()
plot( fD.L2 )

cov_op = cov.operator( fD.L2 )

cov_fun = expand( cov_op )

quartz()
plot( cov_fun, main = 'Cov. function of simulated data', xlab = 'time', ylab = 'time' )

quartz()
image( Matrix( cov( expand( fD.L2 ) ) ),
       useRaster = TRUE,
       colorkey = TRUE,
       col.regions = heat.colors,
       main = 'Sample function',
       xlab = 'time',
       ylab = 'time' )

quartz()
image( Matrix( outer( timeStr$grid, timeStr$grid, expCov ) ),
       useRaster = TRUE,
       colorkey = TRUE,
       col.regions = heat.colors,
       main = 'Exact covariance function',
       xlab = 'time',
       ylab = 'time' )


