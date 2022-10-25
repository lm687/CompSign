
// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // observations (count matrix)
  DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  DATA_MATRIX(z); // matrix for random effects
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  PARAMETER_MATRIX(u_large); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
  PARAMETER_VECTOR(cov_par_RE); // RE. vector with non-redundant entries of the covariance matrix for the random effects variances and covariances
  PARAMETER_VECTOR(logs_sd_RE); // RE. vector with scaling factors for the matrix of covariances of RE
  int d = Y.cols();
  int n = Y.rows();
  int d_min1 = d - 1;

  // covariance matrix for random effects
  using namespace density;

  matrix<Type> theta_prime(n,d_min1); // The probabilities of each event (in ALR)
  vector<Type> Q(n); // The probabilities of each event (marginal of exp of ALR)
  UNSTRUCTURED_CORR_t<Type> nll_mvn(cov_par_RE);

  for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
    nll += VECSCALE(nll_mvn, exp(logs_sd_RE))(u_large.row(i));
  }

  theta_prime = x * beta + z * u_large;


  for(int l=0;l<n;l++){
    Q(l) = 0;
    for(int j=0;j<d_min1;j++){
      Q(l) += exp(theta_prime(l,j));
    }
  }

  matrix<Type> theta(n,d); // The probabilities of each event
  for(int l=0;l<n;l++){
    for(int j=0;j<d_min1;j++){
      theta(l,j) = exp(theta_prime(l,j))/(Type(1.0)+Q(l));
    }
    theta(l,d_min1) = Type(1.0)/(Type(1.0)+Q(l));
  }

  REPORT(theta)

  for(int l=0; l<n; l++){ // Multinomial draws
    vector<Type> lth_row = Y.row(l);
    vector<Type> theta_l = theta.row(l);
    nll -= (dmultinom(lth_row, theta_l, true));
  }
  return nll;

}
