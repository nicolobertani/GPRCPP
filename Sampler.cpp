#include <Rcpp.h>
#include <random>
#include <math.h>
using namespace std;
using namespace Rcpp;
static const double pi = 3.141592653589793238462643383280; 


bool valid_parameters(double delta, double gamma, double rho, double N) {
  bool check_bool = ( ((N*gamma*gamma/(N-1) - gamma/(N-1)) < rho) && (rho < gamma) && (gamma <= delta));
  return check_bool;
}

double logposterior(NumericVector P,  
                    double rho, double gamma, double delta, double mu0, 
                    double sum_P, double sum_P2, double sum_P_outer, double sum_P_outer_plus,
                    double N){
  double ret_log = -std::numeric_limits<double>::infinity();
  //double A, B, C;
  
  if(valid_parameters(delta, gamma, rho, N)){    	    
    double mu = mu0/sqrt(1-gamma);
    double S3 = sum_P2 - 2*mu*sum_P + N*mu*mu;
    double S4 = S3 + 2*sum_P_outer - 2*mu*sum_P_outer_plus + N*(N-1)*mu*mu;
    double eta1 = (delta + (N-1)*rho)/(1-gamma);
    double eta2 = (delta-rho)/(1-gamma);
    ret_log = -N/2*log(2*pi) - 0.5*log(eta1) - (N-1)/2*log(eta2) - (N*S3-S4)/(2*N*eta2) - S4/(2*N*eta1);
    //Rcout << "log_lik before exp value is " << ret_log << std::endl;
    ret_log = ret_log - 1.5*log(eta1)- log(eta2) - 3.5*log(1-gamma); // Jeffrey's Prior
    //ret_log = ret_log + 0.5*log(pow(eta1,2) + pow(eta2,2)*(N-1)) + 0.5*log(pow(eta2,2) + pow(eta1,2)*(N-1)) - 2*log(eta1)- 2*log(eta2) - 3.5*log(1-gamma); // Independence Jeffrey's Prior


    
    //A = ( pow(eta1,2)*pow(eta2,2)*(1-gamma) -  (1-gamma) +  eta1*pow(eta2,2)*0.5*pow(mu0,2));
    //B = (pow(eta1,4)*pow(eta2,4)*pow(1-gamma,7));
    //Rcout << "delta " << delta <<  "rho" << rho << "gamma" << gamma  << std::endl;
    //Rcout << "A " << A << std::endl;
    //Rcout << "B " << B << std::endl;      
    //ret_log = ret_log + 0.5*log(A/B);
    //Rcout << "log_lik before exp value is " << ret_log << std::endl;
    //ret_log = ret_log + 0.5*log( ( pow(eta1,2)*pow(eta2,2)*(1-gamma) -  (1-gamma) +  eta1*pow(eta2,2)*0.5*pow(mu0,2))) - 2*log(eta1) - 2*log(eta1) - 3.5*log(1-gamma);
  }
  return ret_log;
}

// [[Rcpp::export]]
NumericMatrix sample_parameters_cpp(NumericVector p, double p0, int num_sample = 100000, int seed = 1){
    
  int N = p.size();
  double mu0 =  R::qnorm5(p0, 0.0, 1.0, 1, 0);
  NumericVector P = qnorm(p, 0.0, 1.0);
  //NumericVector P = p;
     
  double sum_P = sum(P);
  double sum_P2 = sum(P*P);
  double sum_P_outer = 0;
  double sum_P_outer_plus = 0;
  for(int iter1 = 0; iter1 < N; iter1++){
    for(int iter2 = 0; iter2 < iter1; iter2++){
      sum_P_outer += P[iter1] * P[iter2];
      sum_P_outer_plus += P[iter1] + P[iter2];
    } 
  }
  
    
  default_random_engine generator;
  generator.seed(seed);
  exponential_distribution<double> exp_RV(1.0);
  uniform_real_distribution<double> uni_RV(0.0,1.0);
  normal_distribution<double> normal_RV(0.0,1.0);
  
  NumericMatrix samples(num_sample, 4);
  NumericMatrix::Column p_post = samples(_,0);
  NumericMatrix::Column rhos = samples(_,1);
  NumericMatrix::Column gammas = samples(_,2);
  NumericMatrix::Column deltas = samples(_,3);
	
  // Pick starting points:
  rhos[0] = 0; 
  gammas[0] = 1/((double)N+1); 
  deltas[0] = 2;
   
  double width = 10; // Should not be too small
  double L_delta, H_delta, L_gamma, H_gamma, L_rho, H_rho, u_log_lik, param_log_lik, y;
  double rho_cand, delta_cand, gamma_cand;
  double mu_cond, delta_cond;

  //Rcout << "rho value is " << rhos[0] << std::endl;
  //Rcout << "gamma value is " << gammas[0] << std::endl;
  //Rcout << "delta value is " << deltas[0] << std::endl;
  
  for(int iter = 1;  iter < num_sample; iter++){
	    
    // Sample rho, gamma, delta
    u_log_lik = logposterior(P, rhos[iter-1], gammas[iter-1], deltas[iter-1],  mu0, sum_P, sum_P2, sum_P_outer, sum_P_outer_plus, N);
    //Rcout << "u_log_lik before exp value is " << u_log_lik << std::endl;
    y = u_log_lik - exp_RV(generator);
    //Rcout << "u_log_lik after exp value is " << y << std::endl;
      
    L_delta = deltas[iter-1] - uni_RV(generator)*width;
    H_delta = L_delta + width;      
    L_gamma = 0.0; 
    H_gamma = 1.0;
    L_rho = 0.0; 
    H_rho = 1.0;

    //Rcout << "L_delta " << L_delta << std::endl;
    //Rcout << "H_delta " << H_delta << std::endl;
    //Rcout << "L_gamma " << L_gamma << std::endl;
    //Rcout << "H_gamma " << H_gamma << std::endl;
    //Rcout << "L_rho " << L_rho << std::endl;
    //Rcout << "H_rho " << H_rho << std::endl;
      
    param_log_lik = -std::numeric_limits<double>::infinity();
    while(param_log_lik < y){
      rho_cand = L_rho + uni_RV(generator)*(H_rho-L_rho); //runif(1, L_rho, H_rho)
      gamma_cand = L_gamma + uni_RV(generator)*(H_gamma-L_gamma); //runif(1, L_gamma, H_gamma)
      delta_cand =  L_delta + uni_RV(generator)*(H_delta-L_delta); //runif(1, L_delta, H_delta)
      //Rcout << "rho_cand " << rho_cand <<  "(" << L_rho << "," << H_rho << ")" << std::endl;
      //Rcout << "gamma_cand " << gamma_cand <<  "(" << L_gamma << "," << H_gamma << ")" << std::endl;
      //Rcout << "delta_cand " << delta_cand <<  "(" << L_delta << "," << H_delta << ")" << std::endl;
		
      param_log_lik = logposterior(P, rho_cand, gamma_cand, delta_cand, mu0, sum_P, sum_P2, sum_P_outer, sum_P_outer_plus, N);
      //Rcout << "The value is " << param_log_lik << std::endl;
  
      // Fix dimension if needed:
      if(rho_cand < rhos[iter-1]){
	L_rho = rho_cand;
      } else {
	H_rho = rho_cand;
      }
        
      if(delta_cand < deltas[iter-1]){
	L_delta = delta_cand;
      } else {
	H_delta = delta_cand;
      }

      if(gamma_cand < gammas[iter-1]){
	L_gamma = gamma_cand;
      } else {
	H_gamma = gamma_cand;
      }
    }
	  
    rhos[iter] = rho_cand;
    deltas[iter] = delta_cand;
    gammas[iter] = gamma_cand;
	  
    // Compute posterior of Y:
    //if(bias && (p0 != -1)){
    //	 mu_cond = mu0 +  sqrt(1-gamma_cand)*gamma_cand / (delta_cand + (N-1)*rho_cand) * (sum_P - N*mus[iter] / sqrt(1-gamma_cand));
    //} else {
    mu_cond = mu0 +  sqrt(1-gamma_cand)*gamma_cand / (delta_cand + (N-1)*rho_cand) * (sum_P - N*mu0 / sqrt(1-gamma_cand));
    //}
    delta_cond = 1 - gamma_cand*gamma_cand*N / (delta_cand + (N-1)*rho_cand);
    p_post[iter] = 1-R::pnorm(0, mu_cond, sqrt(delta_cond), 1, 0);	  
  }
	
	
  return samples;
}
