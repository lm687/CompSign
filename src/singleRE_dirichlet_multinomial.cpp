// singleRE_dirichlet_multinomial
// copied from ME_dirichletmultinomial_single_lambda.cpp but adding a second lambda

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
  int d = Y.cols();
  int n = Y.rows();
  // DATA_INTEGER(d) // number of features;
  // DATA_INTEGER(n) // number of samples;
  DATA_INTEGER(num_individuals); // number of different indivdiduals (each with a random effect)
  DATA_MATRIX(lambda_accessory_mat); // matrix to get a length 2*n vector lambda from a length 2 vector lambda
  DATA_MATRIX(x); // matrix of covariates for fixed effects
  DATA_MATRIX(z); // matrix for random effects
  PARAMETER_MATRIX(beta); // coefficients for the fixed effects
  PARAMETER_MATRIX(u_random_effects); // coefficients for the random effects. Even though it is defined as matrix (for TMB matrix multiplication), it is a vector
  PARAMETER(logSigma_RE); // log of the standard deviation of the random effects coefficients
  PARAMETER_VECTOR(log_lambda); // log of the parameter for overdispersion in Dirichlet-Multinomial model (2 values, one for each group)
  int d_min1 = d - 1;

  matrix<Type> log_lambda_vec(n, 1);
  log_lambda_vec = lambda_accessory_mat * log_lambda; // get a vector of lambdas of length 2*n

  for(int i=0;i<num_individuals;i++){
    nll -= dnorm(u_random_effects(i), Type(0.0), exp(logSigma_RE), true);
  }


  matrix<Type> u_large(num_individuals, d_min1); // matrix with columns equal to u, repeated d-1 times
  for(int j=0;j<d_min1;j++){
    u_large.col(j) = u_random_effects.col(0);
  }

  matrix<Type> theta_prime(n,d_min1); // The probabilities of each event (in ALR)
  theta_prime = x * beta + z * u_large;


  vector<Type> Q(n); // The probabilities of each event (marginal of exp of ALR)
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

  // REPORT(theta)

  for(int l=0; l<n; l++){ // Multinomial draws
    vector<Type> lth_row = Y.row(l);
    vector<Type> alpha_l = theta.row(l)*exp(log_lambda_vec(l));
    nll -= dirichlet_multinomial(lth_row, alpha_l, d);
  }
  return nll;

}
