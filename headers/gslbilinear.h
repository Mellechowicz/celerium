#ifndef _GSL_BILINEAR_H
#define _GSL_BILINEAR_H

#include <iostream>
#include <vector>
#include <exception>
#include <gsl/gsl_multiroots.h>
#include "gslmatrix.h"

namespace celerium{
namespace gsl{

enum Type{
	Hybrids,
	Hybrid,
	Newton,
	Broyden};

template<size_t N, typename Callable>
int cFunction(const gsl_vector *x, void *params, gsl_vector *f){
  Vector vector(x);

  for(size_t i=0; i<N; ++i)
    gsl_vector_set (f, i, (*static_cast<Callable*>(params))(vector,i));

  return GSL_SUCCESS;
}

template<size_t N>
class Bilinear{
protected:
	const gsl_multiroot_fsolver_type *solverType;
	gsl_multiroot_fsolver *solver;
	gsl_multiroot_function function;
	Vector result;
	double eps;
	int maxeval;
public:
	Bilinear(double _eps = 1e-5, int _maxeval = 1000):eps(_eps),maxeval(_maxeval){

	}

	~Bilinear(){
		gsl_multiroot_fsolver_free(solver);
	}

	int print_state (size_t iter, gsl_multiroot_fsolver* solver){
		std::cout<<"iter: "<<iter<<std::endl;
		std::cout<<"xs: ";
		for(size_t i=0; i<N; ++i) std::cout<<gsl_vector_get(solver->x, i)<<", ";
		std::cout<<"fs: ";
		for(size_t i=0; i<N; ++i) std::cout<<gsl_vector_get(solver->f, i)<<", ";

		std::cout<<std::endl;

		return GSL_SUCCESS;
	}

	void setEps(double& _eps){
	  eps = _eps;
	}

	void setMaxEval(double& _maxeval){
	  maxeval = _maxeval;
	}

	template<typename Callable>
	int Solve(Callable F, Vector& variables, Type t = Type::Hybrids){
	  switch(t){
	  	case Type::Hybrids:
	  		solverType = gsl_multiroot_fsolver_hybrids;
	  		break;
	  	case Type::Hybrid:
	  		solverType = gsl_multiroot_fsolver_hybrid;
	  		break;
	  	case Type::Newton:
	  		solverType = gsl_multiroot_fsolver_dnewton;
	  		break;
	  	case Type::Broyden:
	  		solverType = gsl_multiroot_fsolver_broyden;
	  		break;
	  	default:
	  		throw std::range_error("Wrong solver given!");
	  }

	  solver = gsl_multiroot_fsolver_alloc(solverType,N);

	  function = {cFunction<N,Callable>, N, static_cast<void*>(&F)};
	  gsl_multiroot_fsolver_set(solver,&function,variables());

	  int status, iter=0;
	  print_state(iter, solver);

	  do{
	    iter++;
	    status = gsl_multiroot_fsolver_iterate(solver);
	    
	    print_state(iter, solver);
	    
	    if (status){
	      printf ("status = %s\n", gsl_strerror(status));
	      gsl_vector_memcpy(variables(),solver->x);
	      return GSL_SUCCESS;
	    }

	    if (iter > maxeval){
	      printf ("status = %s\n", gsl_strerror(status));
	      gsl_vector_memcpy(variables(),solver->x);
	      return GSL_FAILURE;
	    }
	    
	    status = gsl_multiroot_test_residual(solver->f, eps);
	    } while (status == GSL_CONTINUE);
	    
	    printf ("status = %s\n", gsl_strerror(status));
	    gsl_vector_memcpy(variables(),solver->x);
	    return GSL_FAILURE;
	}

}; // end of class Bilinear

} // end of namespace gsl
} // end of namespace celerium
#endif
