rm( list = ls() )

setwd( '~/University/Phd/research/FDA_utilities/my_rfdaR/test/' )

message( ' * * * Fix me! ' )

source( '../basisOp/functional_basis.R')
source( '../geometry/geometry.R')
source( '../fdata/fData.R')
source( '../stats/stats.R' )
source( '../stats/generateFD.R')
source( '../stats/expCovariance.R')
source( '../stats/robustEstimators.R' )



# MEDIAN ------------------------------------------------------------------

timeStr = EvenTimeStr( seq( 0, 1, length.out = 2e2 ) )

bspl.basis = generateBasis( 'Bspline', timeStr, N = 50, degree = 3 )

geometry.L2 = Geometry( bspl.basis, space = 'L2' )

quartz()
plot( bspl.basis )

distroRule = distroRule.Gaussian( mu = rep( 0, bspl.basis$nelements),
                                  Sigma = diag( 0.1, bspl.basis$nelements ) )

N = 1e2
gauss_FD = generateFD( N, distroRule, geometry.L2 ) + sin( 4 * pi * timeStr$grid )


quartz()
plot( gauss_FD, col = brewer.pal( 9, 'Blues' ) )
plot( median.fData2( gauss_FD, ASG = TRUE  ), col = 'yellow', lwd = 2, lty = 1, add = TRUE )
plot( mean( gauss_FD ), col = 'red', lwd = 2, lty = 2, add = TRUE )


# SPHERICAL COVARIANCE ----------------------------------------------------

# Simulating gaussian data with exponential covariance operator

timeStr = EvenTimeStr( seq( 0, 1, length.out = 2e2 ) )

# Amplitude factor
alpha = 0.12

# Correlation decay factor
beta = 0.4

# Number of principal components to retain
L = 10

expCov = ExpCovfunction( alpha, beta )

spectrum = ExpCov_eigencouples( timeStr, L, expCov )

eigen.basis = generateBasis( 'Custom', timeStr, N = length( spectrum$eigen_values),
                             elements = spectrum$eigen_functions )

geometry.L2 = Geometry( eigen.basis, space = 'L2' )

expCovRule.Gauss = distroRule.Gaussian( mu = rep( 0, geometry.L2$basis$nelements ),
                                        Sigma = diag( spectrum$eigen_values ) )

fD = generateFD( 1e3, expCovRule.Gauss, geometry.L2 ) + sin( 4 * pi * timeStr$grid )

quartz()
plot( fD )

# Computing covariance and spherical covariance

C = cov.operator( fD )

CS = covSph( fD, ASG = TRUE )

quartz()
plot( C, expand = TRUE, main = 'Covariance operator' )

quartz()
plot( CS, expand = TRUE, main = 'Spherical covariance' )

# CHECKING EIGENFUNCTIONS -------------------------------------------------

Eigen.C = eigen.Cov_operator( C, nelements = L )
Eigen.CS = eigen.Cov_operator( CS, nelements = L  )

quartz()
par( mfrow = c( 1, 3 ) )
plot( as.fData2.fPC( Eigen.C ), xlab = 'time', ylab = '', main = 'Sample eigenfunctions' )
plot( as.fData2.fPC( Eigen.CS ), xlab = 'time', ylab = '', main = 'Spherical eigenfunctions' )
matplot( timeStr$grid, t( spectrum$eigen_functions * sqrt( spectrum$eigen_values ) ),
         type = 'l', lty = 1, xlab = 'time', ylab = '', main = 'Exact eigenfunctions' )

# Eigenfunctions have same shape but different magnitudes, as a consequence of the
# scalings by eigenvalues. Spherical eigenvalues, in fact, do not coincide with
# sample one.




