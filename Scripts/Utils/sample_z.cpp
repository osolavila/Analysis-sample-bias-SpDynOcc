#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]

IntegerMatrix simulate_occurrence_cpp(NumericVector numN, NumericMatrix neighID, int nyear, int nsite,
  NumericVector field, NumericVector beta0, NumericVector beta_omega, NumericVector beta_gamma) {
  //Sample Z

  double psi;
  double autocov;
  double omega;
  double gamma;
  IntegerMatrix z(nsite, nyear);


  for (int t=0; t<nyear; ++t) {
    for (int i=0; i<nsite; ++i) {
      if (t == 0){
        psi = R::plogis(beta0[0] + beta0[1]*field[i], 0,1,true, false);
      } else{
        int sum_neigh = 0;
        for (int neigh_i=0; neigh_i<numN[i]; ++neigh_i){
          sum_neigh += z(neighID(i,neigh_i), t-1);
        } // neigh_i
        autocov = sum_neigh/numN[i];
        omega = R::plogis(beta_omega[0] + beta_omega[1]*autocov, 0, 1, true, false);
        gamma = R::plogis(beta_gamma[0] + beta_gamma[1]*autocov, 0, 1, true, false);
        psi = z(i,t-1) * omega + (1 - z(i,t-1) ) * gamma;
      }
     
      z(i,t) = R::rbinom(1, psi);
      
    } // i
  } // t

  return(z);

}