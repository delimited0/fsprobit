#ifndef FSPROBIT_HMCMNP_H
#define FSPROBIT_HMCMNP_H

#include <Rcpp.h>

Rcpp::List hmcmnp_moments(
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& Precision,
    int choice_index,
    int n_mc);

#endif
