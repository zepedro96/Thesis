
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


double Lerp(double v0, double v1, double t){
  return (1 - t)*v0 + t*v1;
}

std::vector<double> Quantile(const std::vector<double>& x, const std::vector<double>& probs){
  
  if (x.size() <= 2 || probs.empty()){
    
    throw std::runtime_error("Invalid input");
  }
  
  std::vector<double> data = x;
  std::sort(data.begin(), data.end());
  std::vector<double> quantiles;
  
  for (size_t i = 0; i < probs.size(); ++i){
    
    double center = Lerp(-0.5, data.size() - 0.5, probs[i]);
    
    size_t left = std::max(int64_t(std::floor(center)), int64_t(0));
    size_t right = std::min(int64_t(std::ceil(center)), int64_t(data.size() - 1));
    
    double datLeft = data.at(left);
    double datRight = data.at(right);
    double quantile = Lerp(datLeft, datRight, center - left);
    
    quantiles.push_back(quantile);
  }
  
  return quantiles;
}



std::vector<int> removeDuplicates(std::vector<int> x){

  sort( x.begin(), x.end() );
  x.erase( unique( x.begin(), x.end() ), x.end() );

  return x;
}


// Only works for integer weights since it is based on boolean comparisons...
// A better way would be to implement a threshold based comparison

// [[Rcpp::export]]
NumericVector checkWeights(NumericVector x, 
                           int checkWindow){
  
  // Create new weights vector that will be potentially changed
  // Had to use clone, otherwise just copies by reference...
  NumericVector wts = clone(x);
  
  int L = x.size();
  int inc = std::floor(checkWindow / 2); // symetric increment 
  
  for(int i = 0; i < L; ++i){

    IntegerVector win = IntegerVector::create(i-inc, i+inc);
    win = win[(win >= 0) & (win <= L-1)];
    
    // Weights in the subset window
    NumericVector subsetWeights = x[win];
    
    if(any(x[i] == subsetWeights).is_false()){
      
      double mn = mean(subsetWeights);
      
      if(mn < x[i]){ // Don't want negative bias to have more weight!..
        wts[i] = mn;
      }
    }
  }
  
  return wts;
}



// [[Rcpp::export]]
NumericVector setWeightsCpp(NumericVector x, 
                            IntegerVector indices,
                            NumericVector quantProbs, 
                            NumericVector quantWeights,
                            double baseWeight,
                            bool checkWts,
                            int checkWin){

  // Get unique indices (each value corresponds to a given period, e.g., a year)
  std::vector<int> uniqueIndices = removeDuplicates(Rcpp::as<std::vector<int> >(indices));
  
  int len_x = x.size();
  IntegerVector intIndices = seq_len(len_x);
  
  // Set weights vector with initial value defined in baseWeight
  NumericVector wts(len_x, baseWeight);
  
  int I = uniqueIndices.size();
  int J = quantProbs.size();

  for(int i = 0; i < I; ++i){ // Loop through periods (e.g., year) defined in indices

    // Calculate quantiles for the specific year
    std::vector<double> qtData = Rcpp::as<std::vector<double> >(x[indices == uniqueIndices[i]]);
    std::vector<double> qtsYear = Quantile(qtData,Rcpp::as<std::vector<double> >(quantProbs));
    
    // Integer vector indices for values correspoding to a given period
    IntegerVector loopInd = intIndices[indices == uniqueIndices[i]];
    int K = loopInd.size();
    
    for(int j = 0; j < J; ++j){ // Loop through quantiles
      
      for(int k = 0; k < K; ++k){ // Loop through x values in loopInd
        
        if(x[loopInd[k]-1] >= qtsYear[j]){

          // Set weights if greater than the quantile value
          wts[loopInd[k]-1] = quantWeights[j];
          
        }
      }
    }
  }
  
  if(checkWts){
    NumericVector new_wts = checkWeights(wts, checkWin);
    return(new_wts);
  }else{
    return(wts);
  }
}



