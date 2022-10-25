
// https://figshare.com/articles/C_code_needed_for_running_R_code/4836503
// Simple multinomial fitting.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  Type nll;           // Define variable that holds the return value (neg. log. lik)
  nll = Type(0.0);    // Assign value 1.2; a cast is needed.

  DATA_MATRIX(Y);
  DATA_INTEGER(d);
  DATA_INTEGER(n);
  DATA_MATRIX(x);
  PARAMETER_MATRIX(beta);
  int d_min1 = d - 1;
  int d_min2 = d - 2;

  matrix<Type> theta_prime(n,d_min1); // The probabilities of each event (in ALR)
  theta_prime = x * beta;

  vector<Type> Q(n); // The probabilities of each event (marginal of exp of ALR)
  for(int i=0;i<n;i++){
    Q(i) = 0;
    for(int j=0;j<d_min1;j++){
      Q(i) += exp(theta_prime(i,j));
    }
  }

  matrix<Type> theta(n,d); // The probabilities of each event
  for(int i=0;i<n;i++){
    for(int j=0;j<d_min1;j++){
      theta(i,j) = exp(theta_prime(i,j))/(Type(1.0)+Q(i));
    }
    theta(i,d_min1) = Type(1.0)/(Type(1.0)+Q(i));
  }

  REPORT(theta)

  for(int i=0; i<n; i++){ // Multinomial draws
    vector<Type> ith_row = Y.row(i);
    vector<Type> theta_i = theta.row(i);
    nll += -(dmultinom(ith_row, theta_i, true));
  }
  return nll;

}
