#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]


// Generates an integer sequence similar to R ------------------------------

IntegerVector int_seq(int first, int last) {
  
  IntegerVector y(abs(last - first) + 1);
  
  if (first < last){
    std::iota(y.begin(), y.end(), first);
  } 
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}



// [[Rcpp::export]]

NumericVector HurstFitRollingWindow(IntegerVector M, NumericVector Y, NumericVector Y2){
  
  // S = sqrt(Y2[m]/m - (Y[m]/m)^2)
  // Z = Y[1:m] - (1:m) * Y[m]/m
  
  int L = M.size();
  NumericVector RS(L);
  
  for(int i = 0; i < L; ++i){
    
    int m = M[i];
    double Y2m = Y2[m] / m;
    double Ym = Y[m] / m;
    // S = sqrt(Y2[m]/m - (Y[m]/m)^2)
    // Z = Y[1:m] - (1:m) * Y[m]/m
    double S = sqrt(Y2m - pow(Ym, 2));
    IntegerVector I = int_seq(1, m);
    NumericVector Y_i = Y[I];
    NumericVector Z = Y_i - Rcpp::as<Rcpp::NumericVector >(I) * Y[m] / m;
    RS[i] = (max(Z) - min(Z)) / S;
  }
  return RS;
} 



