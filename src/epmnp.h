#ifndef FSPROBIT_EPMNP_H
#define FSPROBIT_EPMNP_H

#include <Rcpp.h>

void epmnp_moments_inplace(
    Rcpp::NumericVector mu_ep,
    Rcpp::NumericMatrix Sigma_ep,
    int choice_index);

#endif
