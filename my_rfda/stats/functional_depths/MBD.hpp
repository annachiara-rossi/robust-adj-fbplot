//
//  MBD.hpp
//  functional_depths
//
//  Created by Nicholas Tarabelloni on 20/10/15.
//
//

#ifndef MBD_hpp
#define MBD_hpp

#include <stdio.h>
#include <Rcpp.h>

namespace MBD{

  typedef unsigned int UInt;
  typedef double Real;


  typedef Rcpp::NumericVector RcppVector_Type;
  typedef Rcpp::NumericMatrix RcppMatrix_Type;
  typedef std::vector< Real > STDVector_Type;
  typedef Rcpp::NumericVector::iterator RcppVectorIterator_Type;

}


#endif /* MBD_hpp */
