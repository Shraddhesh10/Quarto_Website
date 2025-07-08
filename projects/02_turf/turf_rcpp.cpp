#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame turf_m(NumericMatrix mat, int m, int keep = 10, SEXP w = R_NilValue) {
  int N = mat.nrow(), J = mat.ncol();
  if (m < 1 || m > J) stop("m must be between 1 and ncol(mat)");
  
  // Handle weights
  NumericVector weights = (Rf_isNull(w) || (Rf_length(w) == 1 && as<double>(w) == 1.0))
    ? NumericVector(N, 1.0)
      : as<NumericVector>(w);
  
  if (weights.size() != N || is_true(any(weights < 0)) || sum(weights) == 0 ||
      is_true(any(!is_finite(weights))) || is_true(any(is_na(weights))))
    stop("Invalid weights");
  
  // Generate combinations of m columns
  IntegerMatrix combos = transpose(as<IntegerMatrix>(Function("combn")(J, m)));
  int M = combos.nrow();
  
  NumericVector reach(M), freq(M);
  for (int i = 0; i < M; ++i) {
    NumericVector y(N);
    for (int j = 0; j < m; ++j) {
      int col = combos(i, j) - 1;
      for (int n = 0; n < N; ++n)
        y[n] += mat(n, col);
    }
    
    double rsum = 0, fsum = 0, wsum = sum(weights);
    for (int n = 0; n < N; ++n) {
      if (y[n] > 0) rsum += weights[n];
      fsum += y[n] * weights[n];
    }
    
    reach[i] = rsum / wsum;
    freq[i]  = fsum / wsum;
  }
  
  // Sort by reach and keep top combos
  IntegerVector ord = as<IntegerVector>(Function("order")(reach, _["decreasing"] = true));
  int n = std::min(M, keep);
  ord = ord[Range(0, n - 1)];
  
  // Build output
  List out = List::create(
    _["size"] = rep(m, n),
    _["rank"] = seq_len(n),
    _["reach"] = reach[ord - 1],
                      _["freq"]  = freq[ord - 1]
  );
  
  for (int j = 0; j < m; ++j) {
    IntegerVector items(n);
    for (int i = 0; i < n; ++i)
      items[i] = combos(ord[i] - 1, j);
    out["item" + std::to_string(j + 1)] = items;
  }
  
  return as<DataFrame>(out);
}

// [[Rcpp::export]]
DataFrame turf(NumericMatrix mat, int j = -1, int keep = 10, SEXP w = R_NilValue) {
  if (j == -1) j = std::min(5, mat.ncol() - 1);
  
  List res;
  for (int m = 1; m <= j; ++m)
    res.push_back(turf_m(mat, m, keep, w));
  
  return as<DataFrame>(Function("rbindlist", Environment::namespace_env("data.table"))(res, Named("fill") = true));
}
