
// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>
#include "functions.hpp"

template<class Type>
Type dirichlet_multinomial( vector<Type> Y_row, vector<Type> alpha, int d)
{
  // vector of observed values Y of length d and vector of alpha parameters for a dirichlet distribution of the same length
  Type ll_sum;
  ll_sum = 0;

  Type alpha_sum;
  alpha_sum = 0;
  Type y_sum;
  y_sum = 0;
  Type lgamma_alpha_sum;
  lgamma_alpha_sum = 0;
  Type vector_sum;
  vector_sum = 0;

  for(int j=0; j<d; j++){
    alpha_sum += alpha(j);
    vector_sum += lgamma(alpha(j) + Y_row(j));
    y_sum += Y_row(j);
    lgamma_alpha_sum += lgamma(alpha(j));
  }

  //  lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
  //  - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  ll_sum = lgamma(alpha_sum) + vector_sum - lgamma(alpha_sum + y_sum) - lgamma_alpha_sum;

  return ll_sum;

}

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
  DATA_MATRIX(lambda_accessory_mat); // matrix to get a length 2*n vector lambda from a length 2 vector lambda
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  PARAMETER_MATRIX(u_large); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
  PARAMETER_VECTOR(logs_sd_RE); // RE
  PARAMETER_VECTOR(log_lambda); // log of the parameter for overdispersion in Dirichlet-Multinomial model (2 values, one for each group)
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

  log_lambda_vec = lambda_accessory_mat * log_lambda; // get a vector of lambdas of length 2*n

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
    vector<Type> alpha_l = theta.row(l)*exp(log_lambda_vec(l))*1000;
    nll -= dirichlet_multinomial(lth_row, alpha_l, d);
  }
  return nll;

}
