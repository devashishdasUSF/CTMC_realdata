#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>



/*
******************************************************************
START  Codes from queuecomputer package
******************************************************************
*/

// [[Rcpp::export]]
Rcpp::List qloop_numeric(Rcpp::NumericVector times, Rcpp::NumericVector service, int n_servers) {
  int n = times.size();
  arma::vec output = arma::vec(n);
  arma::Col<int> server_output = arma::Col<int>(n);
  int queue = 0;

  arma::vec queue_times = arma::vec(n_servers);
  queue_times.fill(0);

  for( int i=0; i < n; ++i)
  {
    queue = index_min(queue_times);
    queue_times[queue] = std::max(times[i], queue_times[queue]) + service[i];
    output[i] = queue_times[queue];
    server_output[i] = queue + 1;
    if( i % 512 == 0 )
    {
      Rcpp::checkUserInterrupt();
    }
  }

  // List output_obj;
  // output_obj["times"] = Rcpp::wrap(output);
  // output_obj["servers"] = Rcpp::wrap(server_output);

  std::vector<double> output_stl;
  std::vector<int> server_stl;

  output_stl = arma::conv_to<std::vector<double> >::from(output);
  server_stl = arma::conv_to<std::vector<int> >::from(server_output);

  // Rcpp::Named("servers") = server_stl)

  Rcpp::List output_obj;
  output_obj["times"] = Rcpp::wrap(output);
  output_obj["server"] = Rcpp::wrap(server_output);
  output_obj["state"] = Rcpp::wrap(queue_times);

  return(output_obj);

}

// [[Rcpp::export]]
Rcpp::List qloop_qq(Rcpp::NumericVector times, Rcpp::NumericVector service, 
	Rcpp::NumericVector x, Rcpp::IntegerVector y) {

  int n_servers = max(y);

  arma::vec queue_times = arma::vec(n_servers);
  queue_times.fill(arma::datum::inf);

  for(int i = 0; i < y[0]; i++)
  {
    queue_times[i] = 0;
  }

  int n = times.size();
  arma::vec output = arma::vec(n);
  arma::Col<int> server_output = arma::Col<int>(n);
  output.fill(arma::datum::inf);
  int queue = 0;
  double next_time = x[0];

  int current_size = y[0];
  int next_size = y[1];
  int diff_size = 0;
  int iter = 0;

  for( int i=0; i < n; ++i)
  {

    if( all(queue_times >= next_time) | (times[i] >= next_time))
    {
      diff_size = next_size - current_size;

      if(diff_size > 0)
      {
        for(int j = current_size; j < next_size; j++)
        {
          queue_times[j] = next_time;
        }
      }

      if(diff_size < 0)
      {
        for(int j = next_size; j < current_size; j++)
        {
          queue_times[j] = arma::datum::inf;
        }
      }

      current_size = next_size;
      iter += 1;
      next_size = y[iter+1];
      next_time = x[iter];

    }

    queue = index_min(queue_times);
    queue_times[queue] = std::max(times[i], queue_times[queue]) + service[i];

    output[i] = queue_times[queue];
    server_output[i] = queue + 1;

    // in case user presses stop.
    if( i % 512 == 0 )
    {
      Rcpp::checkUserInterrupt();
    }

    // in case number of servers is zero.
    if( current_size == 0 )
    {
      i = i - 1;
      if( next_time == arma::datum::inf )
      {
        break;
      }
    }

  }

  // List output_obj;
  // output_obj["times"] = Rcpp::wrap(output);
  // output_obj["servers"] = Rcpp::wrap(server_output);

  std::vector<double> output_stl;
  std::vector<int> server_stl;

  output_stl = arma::conv_to<std::vector<double> >::from(output);
  server_stl = arma::conv_to<std::vector<int> >::from(server_output);

  // Rcpp::Named("servers") = server_stl)

  Rcpp::List output_obj;
  output_obj["times"] = Rcpp::wrap(output);
  output_obj["server"] = Rcpp::wrap(server_output);
  output_obj["state"] = Rcpp::wrap(queue_times);

  return(output_obj);

}



/*
******************************************************************
END    Codes from queuecomputer package
******************************************************************
*/