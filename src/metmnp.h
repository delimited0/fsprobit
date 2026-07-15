#ifndef FSPROBIT_METMNP_H
#define FSPROBIT_METMNP_H

#include <Rcpp.h>

Rcpp::List metmnp_moments(
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericMatrix& Sigma,
    int choice_index,
    int n_mc);

#endif
