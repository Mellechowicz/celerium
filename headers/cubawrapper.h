#ifndef _CUBAWRAPPER_H
#define _CUBAWRAPPER_H

#include <array>
#include <algorithm>
#include <functional>
#include <typeindex>
#include <stdexcept>
#include <random>
#include <limits>
#include <cmath>
#include <cuba.h>
#include <iostream>

namespace celerium{
namespace CUBA{


template<auto *function, auto *hyperCube>
static int cFunction(const int *ndim __attribute__((unused)), const double x[] __attribute__((unused)),
	      const int *ncomp __attribute__((unused)), double f[] __attribute__((unused)), void *userdata __attribute__((unused))){

  auto xi = new double(*ndim);
  double jacobian = 1.0;
  for(int i=0; i<*ndim; ++i){
    jacobian *= (hyperCube->operator[](i).second - hyperCube->operator[](i).first);
    xi[i] = hyperCube->operator[](i).first + (hyperCube->operator[](i).second - hyperCube->operator[](i).first)*x[i];
  }

  int result = function(xi,f);

  for(int i=0; i<*ndim; ++i)
	  f[i] *=1.0;// jacobian;

  delete xi;

  return result;
}

template<auto *function>
static int cFunction(const int *ndim __attribute__((unused)), const double x[] __attribute__((unused)),
	      const int *ncomp __attribute__((unused)), double f[] __attribute__((unused)), void *userdata __attribute__((unused))){

  return function(x,f);
}

template<auto *F, size_t NDIM, size_t NCOMP, std::array<std::pair<double,double>,NDIM>* hyperCube = nullptr>
class Cuba{
private:
  int hasFailed;

  std::mt19937 seed_generator;
  std::uniform_int_distribution<> uniform_distribution;

  void* spin;
  double jacobian;

public:
  Cuba(){
    seed_generator = std::mt19937(std::random_device()());
    uniform_distribution = std::uniform_int_distribution<>(0,std::numeric_limits<int>::max());
    spin = NULL;

    jacobian = 1.0;
    for(const auto& limits : *hyperCube)
      jacobian *= limits.second - limits.first;
  }

  Cuba(int _maxeval, int _mineval = 0, double _epsrel = 1.52587890625e-05 /*2^-16*/):
    parameters(_maxeval,_mineval,_epsrel){ 
    seed_generator = std::mt19937(std::random_device()());
    uniform_distribution = std::uniform_int_distribution<>(0,std::numeric_limits<int>::max());
    spin = NULL;

    parameters.epsabs = _epsrel/64;

    jacobian = 1.0;
    for(const auto& limits : *hyperCube)
      jacobian *= limits.second - limits.first;
  }

  struct Parameters{
    const int nDim  = NDIM;
    const int nComp = NCOMP;
    int nvec;
    double epsrel;
    double epsabs;
#ifdef _CUBA_VERBOSE
    const int verbose = 3;
#else
#ifdef _VERBOSE
    const int verbose = 1;
#else
    const int verbose = 0;
#endif
#endif
    const int last = 4;
#ifdef _CUBA_RANDOM_0
    const int level = 1<<8;
#else
#ifdef _CUBA_RANDOM_1
    const int level = 1<<9;
#else
#ifdef _CUBA_RANDOM_2
    const int level = 1<<10;
#else
#ifdef _CUBA_RANDOM_3
    const int level = 1<<11;
#else
#ifdef _CUBA_RANDOM_4
    const int level = 1<<12;
#else
    const int level = 0;
#endif
#endif
#endif
#endif
#endif
    int seed;
    int mineval;
    int maxeval;
    int nstart;
    int nincrease;
    int nbatch;
    int gridno;
    char* statefile; 
    int nnew;
    int nmin;
    double flatness;
    int key1;
    int key2;
    int key3;
    int maxpass;
    double border;
    double maxchisq;
    double mindeviation;
    int ngiven;
    int ldxgiven  = parameters.nDim;
    double* xgiven;
    int nextra;
    int key;

    Parameters(int _maxeval = 1048576*NDIM /*2^20*NDIM*/, int _mineval = 0, double eps = 1.52587890625e-05 /*2^-16*/): 
	nvec(1), epsrel(eps), epsabs(eps*0.015625 /**2^-6*/), seed(0),
       	mineval(_mineval), maxeval(_maxeval), nstart(static_cast<int>(sqrt(maxeval))), nincrease(100), nbatch(100), gridno(0), statefile(NULL),
       	nnew(1000), nmin(2), flatness(5.), key1(47), key2(1), key3(1), maxpass(5), border(0.), maxchisq(10.),
       	mindeviation(.25), ngiven(0), ldxgiven(NDIM), xgiven(NULL), nextra(0), key(0){

 
	}
  } parameters;
  
/*
 *
 * SUAVE
 *
 */
private:
 int suave_explicit(int (*function)(const int*, const double[], const int*, double[], void*), void *userdata,
		    double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
    parameters.seed = uniform_distribution(seed_generator);
    int nregions;

    Suave(parameters.nDim, parameters.nComp, function, userdata, parameters.nvec,
          parameters.epsrel, parameters.epsabs, parameters.verbose | parameters.last | parameters.level, parameters.seed,
          parameters.mineval, parameters.maxeval, parameters.nnew, parameters.nmin, parameters.flatness,
          parameters.statefile, NULL, &nregions, &stepsEvaluated, 
	  &hasFailed, result, errorEstimate, probability);

    for(unsigned i=0; i<NCOMP; ++i){
	    result[i]        *= jacobian;
	    errorEstimate[i] *= jacobian;
    }

    return hasFailed;
  }

public:
 int suave_result(double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
	 if(hyperCube->size()) return suave_explicit(cFunction<F,hyperCube>,NULL,result,errorEstimate,probability,stepsEvaluated);
	 return suave_explicit(cFunction<F>,NULL,result,errorEstimate,probability,stepsEvaluated);

 }

 int suave_result(std::array<double,NCOMP>& result, std::array<double,NCOMP>& errorEstimate, std::array<double,NCOMP>& probability, int& stepsEvaluated){
	 return suave_result(result.data(),errorEstimate.data(),probability.data(),stepsEvaluated);
 }

 /*
  *
  * DIVONNE
  *
  */

private:
 int divonne_explicit(int (*function)(const int*, const double[], const int*, double[], void*), void *userdata,
 		      void (*peakfinder)(const int *, const double [], int *, double [], void*),
		      double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
    parameters.seed = uniform_distribution(seed_generator);
    int nregions=0;

    Divonne(parameters.nDim, parameters.nComp, function, userdata, parameters.nvec,
            parameters.epsrel, parameters.epsabs, parameters.verbose | parameters.last | parameters.level, parameters.seed,
            parameters.mineval, parameters.maxeval, parameters.key1, parameters.key2, parameters.key3,
	    parameters.maxpass, parameters.border, parameters.maxchisq, parameters.mindeviation,
	    parameters.ngiven, parameters.ldxgiven, NULL/*parameters.xgiven*/, parameters.nextra, peakfinder, 
            parameters.statefile, NULL, &nregions, &stepsEvaluated, 
	    &hasFailed, result, errorEstimate, probability);

    for(unsigned i=0; i<NCOMP; ++i){
	    result[i]        *= jacobian;
	    errorEstimate[i] *= jacobian;
    }

    return hasFailed;
  }

public:
 int divonne_result(double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
	 if(hyperCube->size()) return divonne_explicit(cFunction<F,hyperCube>,NULL,NULL,result,errorEstimate,probability,stepsEvaluated);
	 return divonne_explicit(cFunction<F>,NULL,NULL,result,errorEstimate,probability,stepsEvaluated);
 }

 int divonne_result(std::array<double,NCOMP>& result, std::array<double,NCOMP>& errorEstimate, std::array<double,NCOMP>& probability, int& stepsEvaluated){
	 return divonne_result(result.data(),errorEstimate.data(),probability.data(),stepsEvaluated);
 }
}; //end of class Cuba



} //end of namespace CUBA
} //end of namespace celerium

#endif
