#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <iterator>
#include <algorithm>


const double MAX = 1500.00;

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec count_per_interval(NumericVector &t_times, const int Nbasis) {
  // vec taus = logspace<vec>(log(0.0000000001),log(1.0),Nbasis+1);
  vec taus = linspace<vec>(0.0, 1.0, Nbasis+1);
  vec w = diff(taus);
  vec n_ij = zeros<vec>(Nbasis);
  if(t_times.size() > 0) {
    for(int i=0; i<Nbasis;i++) {
      for(int j = 0; j<t_times.size(); j++) {
        if(t_times[j] >= taus[i] && t_times[j] < taus[i+1]) {
          n_ij[i] = n_ij[i] + 1.0;
        }
      }
    } 
  } 
  return n_ij/sqrt(w);
}


double overlap_calculator(double a0, double a1, double b0, double b1){
  if(a1 <= b0) return 0;
  if(a0 >= b1) return 0;
  if(a1 >= b0 && a1 <=b1) {
    if(a0<b0) {
      return (a1 - b0);
    } else {
      return (a1 - a0);
    }
  }
  if(a0 >= b0 && a0 <= b1) {
    if(a1 > b1) {
      return (b1 - a0);
    } else {
      return (a1 - a0);
    }
  }
}

// [[Rcpp::export]]
vec integrate_Y(NumericVector &enter_times, 
  NumericVector &exit_times,
  const int Nbasis) {
  // vec tau = logspace<vec>(log(0.0000000001),log(1.0),Nbasis+1);
  vec tau = linspace<vec>(0.00, 1.00,Nbasis+1);
  vec w = diff(tau);
  vec g_ij = zeros<vec>(Nbasis);
  if(enter_times.size() > 0){
    for(int j = 0; j<enter_times.size(); j++) {
      for(int i=0; i<Nbasis; i++) {
        g_ij[i] = g_ij[i] + overlap_calculator(tau[i],tau[i+1],enter_times[j],exit_times[j]);
      }
    }
  }
  return g_ij/w;
}


// [[Rcpp::export]]
vec n_function(DataFrame df, 
  String transitions, int P = 100, const double MAX = 1500.00) {

  StringVector tran = df["tran"];
  NumericVector Time = df["Time"];

  int N = tran.size();
  vec N_ij(P,fill::zeros);

  NumericVector transition_times;
  for (int j = 0; j < N; ++j)
  {
    if (tran[j] == transitions) {
      transition_times.push_back(Time[j]/MAX);
    }
  }
  N_ij = count_per_interval(transition_times, P);

  return N_ij;
}

// [[Rcpp::export]]
vec y_function(DataFrame df, 
  String s0, int P = 100, const double MAX = 1500.00) {

  StringVector previous_state = df["previous_state"];
  StringVector Event = df["Event"];
  NumericVector Time = df["Time"];
  // StringVector id = df["id"];

  int N = previous_state.size();

  vec G_ij(P,fill::zeros);

  NumericVector enter_times;
  NumericVector exit_times;

  for (int j = 0; j < N; ++j)
  {
    if (previous_state[j] == s0) {
      exit_times.push_back(Time[j]/MAX);
      if(s0 == "arr") {
        enter_times.push_back(0.0);
      } else {    
        enter_times.push_back(Time[j-1]/MAX);
      }
    }
  }
  G_ij = integrate_Y(enter_times, exit_times, P);
  return G_ij;
}

/*
The following function takes a list of data frame as argument and calculates the 
N_ij function for each 
*/
// [[Rcpp::export]]
List get_n_y(List dat, String s0, String transitions, int P = 100, const double MAX = 1500.00) {
  vec n_ij(P,fill::zeros);
  vec g_ij(P,fill::zeros);
  for (int i = 0; i < dat.size(); ++i)
  {
    DataFrame df = dat[i];
    n_ij = n_ij + n_function(df, transitions, P, MAX);
    g_ij = g_ij + y_function(df, s0, P, MAX);
  }
  List result = List::create(_["n_ij"] = n_ij , _["g_ij"] = g_ij);
  return result;
} 

// [[Rcpp::export]]
double generalized_test(DataFrame df, vec beta0, 
  String s0, String transitions, 
  double penalty, int P, const double MAX = 1500.00) {
  vec n_ij = n_function(df, transitions, P, MAX);
  vec g_ij = y_function(df, s0, P, MAX);
  vec gm = penalty - abs(n_ij - beta0 % g_ij);
  vec temp1 = gm % gm;
  return sum(temp1.elem(arma::find(g_ij>0.0))/g_ij.elem(arma::find(g_ij>0.0)));
}

// [[Rcpp::export]]
vec generalized_test_time(DataFrame df, vec beta0, 
  String s0, String transitions, 
  double penalty, int P, const double MAX = 1500.00) {
  vec n_ij = n_function(df, transitions, P, MAX);
  vec g_ij = y_function(df, s0, P, MAX);
  vec gm = penalty - abs(n_ij - beta0 % g_ij);
  vec temp1 = gm % gm;
  vec result(temp1.n_elem, fill::zeros);
  for (int i = 0; i < result.n_elem; ++i)
  {
    if(g_ij[i]>0.0) {
      result[i] = temp1[i]/g_ij[i];
    }
    else {
      result[i] = 0;
    }
  }
  return result;
}




/***R

*/

