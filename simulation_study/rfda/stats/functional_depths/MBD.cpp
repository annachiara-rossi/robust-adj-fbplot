//
//  MBD.cpp
//  functional_depths
//
//  Created by Nicholas Tarabelloni on 20/10/15.
//
//

#include <algorithm>
#include <vector>
#include <iterator>

#include <Rcpp.h>
// #include <RcppArmadillo.h>

#include "MBD.hpp"

using namespace Rcpp;
using namespace MBD;


// [[Rcpp::export]]
RcppVector_Type compute_MBD_2( const RcppMatrix_Type & D )
{
  // Number of observations
  const UInt N( D.rows () );

  // Number of time instants
  const UInt P( D.cols() );

  // Vector of final values of depths
  RcppVector_Type mbd_vector( N, Real( 0 ) );

  for( UInt iCol( 0 ); iCol < P; ++iCol )
  {
    // Extracting the current feature (std::sort works in-place,
    // and views of columns work by references, thus I need to
    // copy the values)
    RcppVector_Type currentColumn( D( _, iCol ) );

    // Sorting the time values
    std::sort( currentColumn.begin(), currentColumn.end() );

    for( UInt iRow( 0 ); iRow < N; ++iRow )
    {
      // If you prefer, declare iter with auto type. I don't do it
      // to avoid warnings on the use of C++11 features
      RcppVectorIterator_Type iter = std::lower_bound( currentColumn.begin(),
                                                       currentColumn.end(),
                                                       D( iRow, iCol )
                                                     );

      const UInt rank( iter - currentColumn.begin() );

      // I add N_up * N_down + N - 1
      mbd_vector[ iRow ] +=  rank * ( N - rank - 1 ) + N - 1;
    }
  }

  mbd_vector =  mbd_vector / ( 0.5 * P * N * ( N - 1 ) );

  return mbd_vector;
}
