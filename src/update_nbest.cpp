#include <Rcpp.h>
using namespace Rcpp;

// algorithm 2 page 81
// Input: pbest and nbest and eucdist of pbest as matrix and pbestval and nbestval as vectors
// updated nbest and nbestval


//' @importFrom Rcpp evalCpp
//' @useDynLib emoptim
// [[Rcpp::export]]

List update_nbest(const NumericMatrix eucdist, const  NumericMatrix pbest,  const NumericMatrix nbest, const NumericVector pbestval, const  NumericVector nbestval,  const double sf) {


  Rcpp::NumericMatrix nbest_copy = Rcpp::clone(nbest);
  Rcpp::NumericVector nbestval_copy = Rcpp::clone(nbestval);

  const int n = pbestval.size();
  double  fertmp = NA_REAL, FER = NA_REAL;

  for (int i=0; i < n; ++i){
    for (int j=0; j < n; ++j){


      if (j == 0)
        fertmp = R_NegInf;

       if (eucdist(i, j) != 0) {

         FER = sf*((pbestval(j)-pbestval(i))/eucdist(i, j));

         if (FER > fertmp){
           fertmp = FER;
         nbest_copy(i, _) = pbest(j, _);
         nbestval_copy(i) = pbestval(j);
         /*Rcpp::Rcout << "FER = " << FER << std::endl; */
         }


       }

    }
  }


  List ret;
  ret["nbest"] = nbest_copy;
  ret["nbestval"] = nbestval_copy;
  return ret;
}











