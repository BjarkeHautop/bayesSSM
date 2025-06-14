#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector resample_multinomial_cpp(int n, NumericVector weights) {
  if (is_true(any(weights < 0))) stop("Weights must be non-negative");
  double total_weight = sum(weights);
  if (total_weight == 0) stop("Sum of weights must be greater than 0");

  NumericVector prob = weights / total_weight;
  IntegerVector indices = Rcpp::sample(n, n, true, prob);
  return indices;
}

// [[Rcpp::export]]
Rcpp::IntegerVector resample_stratified_cpp(int n, Rcpp::NumericVector weights) {
  if (Rcpp::is_true(Rcpp::any(weights < 0))) {
    Rcpp::stop("Weights must be non-negative");
  }
  double total_weight = Rcpp::sum(weights);
  if (total_weight == 0) {
    Rcpp::stop("Sum of weights must be greater than 0");
  }
  NumericVector prob = weights / total_weight;
  NumericVector cum_sum = cumsum(prob);
  IntegerVector u = seq(0, n - 1);
  NumericVector u_num = as<NumericVector>(u);
  NumericVector u_scaled = (u_num + runif(n)) / n;

  IntegerVector indices(n);
  int j = 0;
  for (int i = 0; i < n; i++) {
    while (j < (int)cum_sum.size() - 1 && cum_sum[j] < u_scaled[i]) {
      j++;
    }
    indices[i] = j + 1;  // shift by 1 here for 1-based indexing
  }

  return indices;
}

// [[Rcpp::export]]
Rcpp::IntegerVector resample_systematic_cpp(int n, Rcpp::NumericVector weights) {
  if (Rcpp::is_true(Rcpp::any(weights < 0))) {
    Rcpp::stop("Weights must be non-negative");
  }
  double total_weight = Rcpp::sum(weights);
  if (total_weight == 0) {
    Rcpp::stop("Sum of weights must be greater than 0");
  }
  NumericVector prob = weights / total_weight;
  NumericVector cum_sum = cumsum(prob);
  IntegerVector u = seq(0, n - 1);
  NumericVector u_num = as<NumericVector>(u);
  NumericVector u_scaled = (u_num + R::runif(0.0, 1.0)) / n;
  IntegerVector indices(n);
  int j = 0;
  for (int i = 0; i < n; i++) {
    while (j < (int)cum_sum.size() - 1 && cum_sum[j] < u_scaled[i]) {
      j++;
    }
    indices[i] = j + 1;  // shift by 1 here for 1-based indexing
  }

  return indices;
}
