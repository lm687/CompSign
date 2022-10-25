
// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y); // observations (count matrix)
  int d = Y.cols(); // number of features
  int n = Y.rows(); // 2*n samples
  DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  DATA_MATRIX(z); // matrix for random effects
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  PARAMETER_MATRIX(u_large); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
  PARAMETER_VECTOR(logs_sd_RE); // RE
  int d_min1 = d - 1;

  // for covariance matrix for random effects
  using namespace density;

  matrix<Type> log_lambda_vec(n, 1);
  matrix<Type> theta_prime(n,d_min1); // The probabilities of each event (in ALR)
  matrix<Type> theta(n,d); // The probabilities of each event
  matrix<Type> Sigma(d_min1,d_min1);
  Sigma.fill(0);                   // Fill the whole matrix
  Sigma.diagonal() = exp(logs_sd_RE);          // Multiply diagonal by 10 to positive definite Sigma
  MVNORM_t<Type> nll_mvn(Sigma);
  
  for(int i=0;i<num_individuals;i++){ // likelihood for the random effects (multivariate normal)
    nll += nll_mvn(u_large.row(i));
  }

  theta_prime = x * beta + z * u_large;

  vector<Type> Q(n); // The probabilities of each event (marginal of exp of ALR)
  for(int l=0;l<n;l++){
    Q(l) = 0;
    for(int j=0;j<d_min1;j++){
      Q(l) += exp(theta_prime(l,j));
    }
  }

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
//     nll -= multinomial(lth_row, theta_l, d); // this was what was written on 20201025 and ofc doesn't compile (no multinomila() function). I also don't know what this d refers to
    nll -= dmultinom(lth_row, theta_l, true);
  }
  return nll;

}
