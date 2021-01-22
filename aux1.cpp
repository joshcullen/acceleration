// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
// [[Rcpp::export]]
IntegerMatrix SummarizeDat(IntegerVector z, IntegerVector dat, int ncateg, 
                           int nbehav, int nobs){
  IntegerMatrix nmat(nbehav,ncateg);
  LogicalVector na1=!is_na(dat);
  
  for (int i = 0; i < nobs; i++){
    if (na1[i]){
      nmat(z[i],dat[i])=nmat(z[i],dat[i])+1;
    }
  }
  return nmat;
}

// This function helps with multinomial draws
// [[Rcpp::export]]
int cat1(double value, NumericVector prob) {
  int res=prob.length()-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector z(prob.nrow());
  
  for(int i=0; i<prob.nrow();i++){
    z[i]=cat1(randu[i],prob(i,_));
  }
  return z;
}

//' This function helps store z from all iterations after burn in
// [[Rcpp::export]]
IntegerMatrix StoreZ(IntegerVector z, IntegerMatrix store_z, int nobs) {
  
  for(int i=0; i<nobs;i++){
    store_z(i,z[i])=store_z(i,z[i])+1;
  }
  return store_z;
}