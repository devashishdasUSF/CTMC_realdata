#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>


Rcpp::IntegerVector whichC(Rcpp::LogicalVector x) {
  int count = 0;
  for(int i = 0; i<x.size(); i++) {
    if(x[i]) {
      count++;
    }
  }
  Rcpp::IntegerVector idx(count);

  count = 0;
  for(int i = 0; i<x.size(); i++) {
    if(x[i]) {
      idx[count] = i;
      count++;
    }
  }
  return idx;
}


// [[Rcpp::export]]
double intensityFunc(double t)
{
	return 5*(2 - sin(2*M_PI*t));// +  sin(M_PI*(t-0.75)/2));
}

// [[Rcpp::export]]
double service_rate(double t, double shift = 0.0, double slow = 1.0)
{
  double x = 1.1*intensityFunc(t -shift)/slow;
  return -log((Rcpp::runif(1))[0])/x;
}



// [[Rcpp::export]]
Rcpp::NumericVector homPoisson(double maxrate=50) {
	Rcpp::NumericVector times(100);
	double t = 0, s =0;
	int counter = 0;
	while (t <= 1) {
		t = t - log((Rcpp::runif(1))[0])/maxrate;
		times[counter] = t; counter++;
	}
	Rcpp::NumericVector empt(0);
	if (counter-2 > 0) {
		return times[Rcpp::Range(0, counter-2)];
		// return Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(times.subvec(0, counter-2)));
	} else {
		return empt;
	}
}


// [[Rcpp::export]]
Rcpp::NumericVector inhomPoisson(double maxrate=50) {
	Rcpp::NumericVector times = homPoisson(maxrate);
	Rcpp::LogicalVector index(times.size(), FALSE);
	for(int i = 0; i<times.size(); i++) {
		if((Rcpp::runif(1))[0] <= intensityFunc(times[i])/maxrate) {
			index[i] = TRUE;
		}
	}
	return times[index];
}


// [[Rcpp::export]]
Rcpp::List queue_new(int n_servers = 2, double shift = 0.0, double slow = 1.0) {

  arma::vec arrival = Rcpp::as<arma::vec>(inhomPoisson());
  int n = arrival.size();
  arma::vec departure(n);
  int queue = 0;
  

  arma::vec queue_times(n_servers);
  queue_times.fill(0);
  arma::vec start_time(n);

  for( int i=0; i < n; ++i)
  {
    queue = index_min(queue_times);
    start_time[i] = std::max(arrival[i], queue_times[queue]);
    queue_times[queue] =  start_time[i] + service_rate(start_time[i], shift, slow);
    departure[i] = queue_times[queue];
  }


  // Rcpp::DataFrame time_df = Rcpp::DataFrame::create(Rcpp::Named("arrive_times") = arrival,Rcpp::Named("proc_start_times") = start_time, Rcpp::Named("wait_times") = start_time - arrival, Rcpp::Named("depart_times") = departure);
  
  Rcpp::List output_obj;
  // output_obj["time_df"] = time_df;

  arma::ivec x(n); x.fill(1);
  arma::ivec y(n); y.fill(-1);
  arma::ivec z = join_cols(x,y);
  

  arma::vec event_times = join_cols(arrival, departure);
  z = z.elem(arma::sort_index(event_times));
  event_times = event_times.elem(arma::sort_index(event_times));
  arma::ivec number_in_q = cumsum(z);

  arma::uvec ind;
  ind = arma::find(event_times <= 1.0);

  number_in_q = number_in_q.elem(ind);
  z = z.elem(ind);
  event_times = event_times.elem(ind);

  output_obj["number_in_q"] = number_in_q;
  output_obj["arrival_or_departure"] = z;
  output_obj["times"] = event_times;

  return(output_obj);

}


// [[Rcpp::export]]
arma::vec count_per_interval(Rcpp::NumericVector &t_times, const int Nbasis) {
  // vec taus = logspace<vec>(log(0.0000000001),log(1.0),Nbasis+1);
  arma::vec taus = arma::linspace<arma::vec>(0.0, 1.0, Nbasis+1);
  arma::vec w = arma::diff(taus);
  arma::vec n_ij = arma::zeros<arma::vec>(Nbasis);
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
arma::vec integrate_Y(Rcpp::NumericVector &enter_times, Rcpp::NumericVector &exit_times,
  const int Nbasis) {
  // vec tau = logspace<vec>(log(0.0000000001),log(1.0),Nbasis+1);
  arma::vec tau = arma::linspace<arma::vec>(0.00, 1.00,Nbasis+1);
  arma::vec w = arma::diff(tau);
  arma::vec g_ij = arma::zeros<arma::vec>(Nbasis);
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
Rcpp::List counting_model(Rcpp::List queue_obj) {


  Rcpp::NumericVector number_in_q;
  Rcpp::NumericVector z;
  Rcpp::NumericVector event_times;

  number_in_q = queue_obj["number_in_q"];
  z = queue_obj["arrival_or_departure"];
  event_times = queue_obj["times"];

  Rcpp::NumericVector empty_times;
  empty_times = event_times[number_in_q == 0];


  Rcpp::NumericVector enter_times;
  Rcpp::NumericVector exit_times;
  Rcpp::NumericVector dep_times;

  int num_events = event_times.size();

  if(num_events == 0) { // if nothing happened;

    enter_times = Rcpp::NumericVector(0);
    exit_times = Rcpp::NumericVector(0);
    dep_times = Rcpp::NumericVector(0);

  } else {
    dep_times = event_times[z == -1];
    if(empty_times.size() == 0) {   // no empty times

      enter_times = Rcpp::NumericVector::create(event_times[0]);
      exit_times = Rcpp::NumericVector::create(1.0);

    } else {
      int sz = empty_times.size();
      Rcpp::IntegerVector idx(sz);
      idx = whichC(number_in_q == 0);
      if(number_in_q[num_events-1] == 0) { // last event is a empty queue

        enter_times = Rcpp::NumericVector(sz);
        exit_times = Rcpp::NumericVector(sz);
        enter_times[0] = event_times[0];
        for(int i = 1; i<sz; i++) enter_times[i] = event_times[idx[i-1]+1];
        for(int i = 0; i<sz; i++) exit_times[i] = event_times[idx[i]];

      } else {

        enter_times = Rcpp::NumericVector(sz + 1);
        exit_times = Rcpp::NumericVector(sz + 1);
        enter_times[0] = event_times[0];
        for(int i = 1; i<sz; i++) enter_times[i] = event_times[idx[i-1]+1];
        for(int i = 0; i<sz; i++) exit_times[i] = event_times[idx[i]];

        enter_times[sz] = event_times[idx[sz - 1] + 1];
        exit_times[sz] = 1.0;

      }
    }
  }

  Rcpp::List EnterExitDep;
  EnterExitDep["enter"] = enter_times;
  EnterExitDep["exit"] = exit_times;
  EnterExitDep["departures"] = dep_times;
  return EnterExitDep;



  // last entry is empty;

  
}

// [[Rcpp::export]]
Rcpp::List get_n_and_y(Rcpp::List EnterExitDep_obj, int Nbasis = 1000) {

  Rcpp::NumericVector enter_times;
  Rcpp::NumericVector exit_times;
  Rcpp::NumericVector dep_times;

  enter_times = EnterExitDep_obj["enter"];
  exit_times = EnterExitDep_obj["exit"];
  dep_times = EnterExitDep_obj["departures"];


  arma::vec N;
  arma::vec Y;

  N = count_per_interval(dep_times, Nbasis);
  Y = integrate_Y(enter_times, exit_times, Nbasis);

  Rcpp::List NY_object;
  NY_object["N"] = N;
  NY_object["Y"] = Y;

  return NY_object;

}

// [[Rcpp::export]]
Rcpp::List simulate() {

  Rcpp::List x;
  Rcpp::List y;
  Rcpp::List z;

  int Nbasis = 100;
  int ITER = 100000;

  arma::vec N(Nbasis);
  arma::vec Y(Nbasis);

  for(int i = 0; i<ITER; i++) {
    x = queue_new();
    y = counting_model(x);
    z = get_n_and_y(y, Nbasis);
    arma::vec N0 = z[0];
    arma::vec Y0 = z[1];
    N = N + N0;
    Y = Y + Y0;
  }

  Rcpp::List NY_object;
  NY_object["N"] = N/ITER;
  NY_object["Y"] = Y/ITER;

  return NY_object;

}

arma::vec KS_test_cdf(Rcpp::NumericVector dep_times, int nlin = 1000) {
  int n_dep;
  arma::vec a;
  a = arma::linspace<arma::vec>(0.,1.,nlin);
  arma::vec cdf(nlin);
  cdf.zeros();

  n_dep = dep_times.size();
  if(n_dep > 0) {
    for(int i = 0; i<nlin; i++) {
      cdf[i] = (double)sum(dep_times  <= a[i])/n_dep;
    }  
  } else {
    for(int i = 0; i<nlin; i++) {
      cdf[i] = 0.0;
    }
  }
  return cdf;
}

//[[Rcpp::export]]
arma::vec KS_train() {
  Rcpp::List x;
  Rcpp::List y;
  Rcpp::NumericVector dep_times;
  double nlin = 1000;
  arma::vec cdf(nlin);
  cdf.zeros();
  
  int ITER = 100000;

  for(int I = 0; I<ITER; I++) {
    x = queue_new();
    y = counting_model(x);
    dep_times = y["departures"];
    cdf = cdf +  KS_test_cdf(dep_times,nlin);
  } 
  return cdf/ITER;
}

//[[Rcpp::export]]
arma::mat KS_test_stat(arma::vec KS0_intensity, int n_servers = 2, 
  double shift = 0.0, double slow = 1.0) {

  Rcpp::List x;
  Rcpp::List y;
  Rcpp::NumericVector dep_times;

  arma::vec cdf;

  int Nsample = 100000;
  arma::vec stat0(Nsample);
  for(int i = 0; i < Nsample; i++) {
    x = queue_new(n_servers, shift, slow);
    y = counting_model(x);
    dep_times = y["departures"];
    cdf = KS_test_cdf(dep_times, KS0_intensity.size());
    stat0[i] = max(abs(cdf - KS0_intensity));
  }
  
  return stat0;
} 


//[[Rcpp::export]]
arma::mat test_statistic(arma::vec Beta0, int n_servers = 2, double shift = 0.0, double slow = 1.0) {

  Rcpp::List x;
  Rcpp::List y;
  Rcpp::List z;

  int Nbasis = 100;
  int ITER = 100000;


  // try to test for slowing down of queue network
  arma::vec Beta1 = Beta0/1.1;

  arma::mat G(Nbasis,Nbasis);
  G.fill(0.00);
  

  x = queue_new(n_servers, shift, slow);
  y = counting_model(x);
  z = get_n_and_y(y, Nbasis);


  arma::vec N = z[0];
  arma::vec Y = z[1];



  G.diag() = Y;


  return -(-2*Beta0.t()*N + Beta0.t()*G*Beta0 + 2*Beta1.t()*N - Beta1.t()*G*Beta1);

}


//[[Rcpp::export]]
arma::mat test_stat(arma::vec Beta0, int n_servers = 2, double shift = 0.0, double slow = 1.0) {
  int Nsample = 100000;
  arma::vec stat0(Nsample);
  for(int i = 0; i < Nsample; i++) {
    stat0[i] = (test_statistic(Beta0, n_servers, shift, slow))[0,0];
  }
  return stat0;
}




/*** R
library(microbenchmark)
library(dplyr)



*/