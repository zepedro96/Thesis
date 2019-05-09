

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]


// Generates an integer sequence similar to R ------------------------------

IntegerVector int_seq(int &first, int &last) {
  
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


// Hampel filter/identifier -----------------------------------------------


// [[Rcpp::export]]

NumericVector HampelFilterCpp(NumericVector x, int k, double t0, int winType){
  
  // Check for NA's first
  if(all(is_na(x)).is_false()){
    
    NumericVector res = clone(x);
    int len = x.size();
    int it_min(0);
    int it_max(0);
    
    if(winType == 1){ // symmetric window
      it_max = it_max + len;
    }
    else if(winType == 2){ // right-side window
      it_min = it_min + k;
      it_max = it_max + (len - k);
    } 
    else{
      return NumericVector::create(0);
    }

    for(int i = it_min; i < it_max; ++i){
      
      int first = i - k;
      int last = i + k;
      
      if(winType==1){
        
        if(first < 0) // Adjust first window element
          first = 0;
        
        if(last > (len - 1)) // Adjust last window element
          last = len - 1;
      }
      
      IntegerVector win = int_seq(first, last); // integer vector with window range
      NumericVector w0 = res[win]; // data in test window
      double x0 = median(w0); // median in window (using the already corrected data)
      NumericVector absDiff = abs(w0 - x0);
      // 1.4826 is a constant scale factor, which depends on the distribution
      // For normally distributed data k is taken to be:  k = [1 / ( Φ^−1 (3 / 4))] ≈ 1.4826
      double s0 = 1.4826 * median(absDiff); // Calculate normalized median absolute differences

      if(abs(res[i] - x0) > (t0 * s0)){
        res[i] = x0;
      }
    }
    
    return res;
    
  }
  else{
    return NumericVector::create(0);
  }
}




  
  