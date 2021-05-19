// Negative Log-likelihood
#include <TMB.hpp>

using namespace Eigen;
using namespace tmbutils;

// https://github.com/kaskr/adcomp/issues/96
/* List of matrices struct*/
template<class Type>
struct LOSM_t : vector<matrix<Type> > {
  LOSM_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRUCT(Xlist, LOSM_t);
  DATA_IVECTOR(censorvec);
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(pi);
  
  Type nll = 0;
  Type jll = 0;
  // Calculate and sum individual log likelihoods
  for(int i=0; i < censorvec.size(); i++){
	  matrix<Type> X = Xlist(i);
	  Type indlik = sum(exp(X*alpha));
	  jll = 0;
	  if(censorvec(i)){
		  for(int j=0; j<theta.size(); j++){
			  jll += pi(j)*exp(-exp(theta(j))*indlik);
		  }
		  nll -= log(jll);
	  }
	  else{
		  int num_cols = X.cols();
	          int num_rows = X.rows() - 1;
		  matrix<Type> Xc = X.block(0,0,num_rows,num_cols);
		  Type indlikc = sum(exp(Xc*alpha));
		  for(int j=0; j<theta.size(); j++){
			  jll += pi(j)*exp(-exp(theta(j)) * indlikc);
			  jll -= pi(j)*exp(-exp(theta(j)) * indlik);
		  }
		  nll -= log(jll);
	 }
  }
  ADREPORT(alpha);
  ADREPORT(theta);
  ADREPORT(pi);
  
  return nll;
}
