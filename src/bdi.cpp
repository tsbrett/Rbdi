#include <Rcpp.h>
#include <math.h>   //for log, exp etc
#include <algorithm> //for min functionality
#include <iostream> //unneeded?

using namespace Rcpp;


//Implement choose function from R:
// [[Rcpp::export]]
double c_choose(double N, double K){
	double result = 1;
	for (int i = 1; i <= K; i++){
	    result *= N - (K - i);
	    result /= i;
	}
	return result;
}


// Can I generalise this for the release version to arbitrary functions
//Infectiousness function. This function could be moved to the R script
double Cbeta(double beta0, double beta1, double t ){
	  return beta0 + beta1 *t;
	}

//	#################################################
//	### Calculating the probability distribution: ###
//	#################################################


// [[Rcpp::export]]
double P_BDC(double x, double x0,
              double lambda, double mu, double t){

  double prob_BD = 0.0;
  double A_BD = mu*(exp((lambda-mu)*t) -1)/(lambda*exp((lambda-mu)*t) - mu); //tick
  double B_BD = lambda*(exp((lambda-mu)*t) -1)/(lambda*exp((lambda-mu)*t) - mu); //tick

  // Rcout << "A_BD " << A_BD << "\n";
  // Rcout << "B_BD " << B_BD << "\n";
  // Rcout << "A_BD +B_BD " << A_BD+B_BD << "\n";



  if((x0 < 0) || (x < 0)){
    Rcerr << "Either x or x0 is less than 0, returning 0" << "\n";
    prob_BD = 0.0;
    }

  else if(x0 == 0){
    // Rcout << "x0 == 0" << "\n";
    if(x == 0){
      prob_BD = 1.0;
    }
    else{
      prob_BD = 0.0;
    }
  }

  else if(x == 0){
    // Rcout << "x == 0" << "\n";
    prob_BD = pow(A_BD,x0);
  }

  else{
    double min_x = std::min(x,x0);
    double prob_j = 0;
    // Rcout << "min_x:" << min_x <<"\n";

    for (int j = 0; j <= min_x; j++ )
    {
      prob_j = c_choose(x0,j)*c_choose(x0 + x -j - 1, x0 -1)* //tick
        pow(A_BD,x0-j)* //tick
        pow(B_BD,x-j)* //tick
        pow((1.0 - A_BD-B_BD),j); //tick
      // Rcout << "j: " << j << " prob_j: " << prob_j <<"\n";

      prob_BD += prob_j;
    //  if(abs(prob_j) > eps){
    //    Rcerr << "BD calculation greater than tolerance (eps =" << eps <<"), returning NaN";
    //    prob_BD = NAN;
    //    break;
    //  }
    }
  }

  // Rcout << "prob_BD: " << prob_BD << "\n";
  return prob_BD;
}


// [[Rcpp::export]]
double P_BDIC(double x, double x0, double v,
              double lambda, double mu, double t){

  //Todo: handling of lambda, v, mu <= 0.
  //Todo: handling of lambda == mu

  double m_BDI = (v/(lambda-mu))*(exp((lambda-mu)*t) -1.0); //tick
  double r_BDI = v/lambda;
  double p_BDI = m_BDI/(m_BDI +r_BDI);

  double prob_BD = 0.0;
  double prob_BDI0 = 0.0;
  double prob_BDI = 0.0;
  double min_x = 0;
  if(x0 == 0){
    // Rcout << "x0==0; no BD" << "\n";

    // No contribution from BD process, just immigration:
    prob_BDI = c_choose(x+r_BDI-1,x)*pow(p_BDI,x)*pow((1-p_BDI),r_BDI);
  }
  else{
    for (int i = 0; i <= x; i++){
      // Rcout << "i: " << i << "\n";

      prob_BD = P_BDC(x-i, x0, lambda, mu, t);

      //Calculate prob_BDI0 -- contribution to prob_BDI from new importations,
      // given by a negative binomial distribution
      prob_BDI0 = c_choose(i+r_BDI-1,i)*pow(p_BDI,i)*pow((1-p_BDI),r_BDI); // tick
      // Rcout << "prob_BDI0: " << prob_BDI0 << "\n";

      prob_BDI += prob_BD*prob_BDI0; //tick
    }
  }

  return prob_BDI;
}

//	#################################################
//	### Calculating the log-likelihood: ###
//	#################################################

// [[Rcpp::export]]
double C_bdi_ll_linear_R0(NumericVector x, double delta_t,
                          double eta, double gamm, double dR0, double R00){
	int n = x.size();
	//double eps = std::sqrt(5)*pow(10, -10); // todo: this is not robust!
	double LogLikelihood = 0.0;
	double lambda, p_m;
	for(int m =1; m < n; m++){ // seems correct
	  lambda = Cbeta(R00, dR0, (m-1)*delta_t)*gamm;
	  p_m = P_BDIC(x[m], x[m-1], eta, lambda, gamm, delta_t);
	  LogLikelihood += std::log(p_m);

		// LogLikelihood += std::log(P_BDIC(x[m], x[m-1], eta, lambda, gamm, delta_t));
    //LogLikelihood += std::log(P_BDIC(x[m],x[m-1], eta/gamm, Cbeta(R00 +eps , dR0+eps,m*delta_t), 1, delta_t) );
      }
      return LogLikelihood;
  }











