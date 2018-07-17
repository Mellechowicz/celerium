#ifndef _NEWCUBAWRAPPER_H
#define _NEWCUBAWRAPPER_H

#include <array>
#include <vector>
#include <algorithm>
#include <functional>
#include <typeindex>
#include <stdexcept>
#include <random>
#include <limits>
#include <memory>
#include <cmath>
#include <cuba.h>
#include <iostream>
#include <tuple>

namespace celerium{
namespace cuba{

template<typename Callable>
int cFunction(const int *ndim __attribute__((unused)), const double x[] __attribute__((unused)),
	      const int *ncomp __attribute__((unused)), double f[] __attribute__((unused)), void *userdata __attribute__((unused))){

  auto function = std::get<0>(*static_cast<std::tuple<Callable&,std::vector<std::pair<double,double>>&,int&>*>(userdata));
  auto hyperCube = std::get<1>(*static_cast<std::tuple<Callable&,std::vector<std::pair<double,double>>&,int&>*>(userdata));
  auto nVec = std::get<2>(*static_cast<std::tuple<Callable&,std::vector<std::pair<double,double>>&,int&>*>(userdata));

  int result = -1;

    double *xi = new double[nVec*hyperCube.size()];
    for(int i=0; i<nVec; ++i)
      for(size_t j=0U; j<hyperCube.size(); ++j){
        xi[hyperCube.size()*i+j] = hyperCube[j].first + (hyperCube[j].second - hyperCube[j].first)*x[i*hyperCube.size()+j];
      }

    result = function(xi,f);

    delete [] xi;
  
  return result;
}

class Cuba{
private:
  int hasFailed;

  std::mt19937 seed_generator;
  std::uniform_int_distribution<> uniform_distribution;

  void* spin;
  double jacobian;

public:
  Cuba(int _maxeval = 100, int _mineval = 0, double _epsrel = 1.52587890625e-05 /*2^-16*/);
  virtual ~Cuba();

  struct Parameters{
    int nDim;
    int nComp;
    const int nvec = 1;
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
    int ldxgiven;
    double* xgiven;
    int nextra;
    int key;

    Parameters(int _nDim = 1, int _maxeval = 1048576 /*2^20*/, int _mineval = 0, double eps = 1.52587890625e-05 /*2^-16*/): 
	nDim(_nDim), nComp(1), epsrel(eps), epsabs(eps*0.015625 /**2^-6*/), seed(0),
       	mineval(_mineval), maxeval(_maxeval), nstart(static_cast<int>(sqrt(maxeval))), nincrease(100), nbatch(100), gridno(0), statefile(NULL),
       	nnew(1000), nmin(2), flatness(5.), key1(47), key2(1), key3(1), maxpass(5), border(0.), maxchisq(10.),
       	mindeviation(.25), ngiven(0), ldxgiven(_nDim), xgiven(NULL), nextra(0), key(0){}
  } parameters;
  
/*
 *
 * SUAVE
 *
 */
private:
 int suave_explicit(int (*function)(const int*, const double[], const int*, double[], void*), void *userdata,
		    double result[], double errorEstimate[], double probability[], int& stepsEvaluated);

public:
 template<typename Callable>
 int suave_result(const Callable& F, std::vector<std::pair<double,double>> _hyperCube, int nComp, double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
	parameters.nDim = _hyperCube.size();
        parameters.nComp = nComp;

	jacobian = 1.0;
	for(const auto& limits : _hyperCube)
	  jacobian *= limits.second - limits.first;

	//auto data = std::make_tuple(F,_hyperCube,parameters.nvec);
        auto data = std::tie(F,_hyperCube,parameters.nvec);

	auto output = suave_explicit(cFunction<Callable>,static_cast<void*>(&data),result,errorEstimate,probability,stepsEvaluated);
	return output;

 }

template<typename Callable, size_t NCOMP>
 int suave_result(const Callable& F, std::vector<std::pair<double,double>> _hyperCube, std::array<double,NCOMP>& result, std::array<double,NCOMP>& errorEstimate, std::array<double,NCOMP>& probability, int& stepsEvaluated){
  return suave_result<Callable>(F,_hyperCube,NCOMP,result.data(),errorEstimate.data(),probability.data(),stepsEvaluated);
 }

template<typename Callable>
 int suave_result(const Callable& F, std::vector<std::pair<double,double>> _hyperCube, std::vector<double>& result, std::vector<double>& errorEstimate, std::vector<double>& probability, int& stepsEvaluated){

  if (result.size() == 0)
    throw std::invalid_argument("celerium::cuba::Cuba::divonne_result: result vector must have non-zero size.");
  
  errorEstimate.resize(result.size());
  probability.resize(result.size());

  return suave_result<Callable>(F,_hyperCube,result.size(),result.data(),errorEstimate.data(),probability.data(),stepsEvaluated);
 }

  
 /*
  *
  * DIVONNE
  *
  */

private:
 int divonne_explicit(int (*function)(const int*, const double[], const int*, double[], void*), void *userdata,
		    double result[], double errorEstimate[], double probability[], int& stepsEvaluated);

public:
 template<typename Callable>
 int divonne_result(const Callable& F, std::vector<std::pair<double,double>> _hyperCube, int nComp, double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
	parameters.nDim = _hyperCube.size();
        parameters.nComp = nComp;

	jacobian = 1.0;
	for(const auto& limits : _hyperCube)
	  jacobian *= limits.second - limits.first;

	//auto data = std::make_tuple(F,_hyperCube,parameters.nvec);
        auto data = std::tie(F,_hyperCube,parameters.nvec);

	auto output = divonne_explicit(cFunction<Callable>,static_cast<void*>(&data),result,errorEstimate,probability,stepsEvaluated);
	return output;

 }

template<typename Callable, size_t NCOMP>
 int divonne_result(const Callable& F, std::vector<std::pair<double,double>> _hyperCube, std::array<double,NCOMP>& result, std::array<double,NCOMP>& errorEstimate, std::array<double,NCOMP>& probability, int& stepsEvaluated){
  return divonne_result<Callable>(F,_hyperCube,NCOMP,result.data(),errorEstimate.data(),probability.data(),stepsEvaluated);
 }

template<typename Callable>
 int divonne_result(const Callable& F, std::vector<std::pair<double,double>> _hyperCube, std::vector<double>& result, std::vector<double>& errorEstimate, std::vector<double>& probability, int& stepsEvaluated){

  if (result.size() == 0)
    throw std::invalid_argument("celerium::cuba::Cuba::divonne_result: result vector must have non-zero size.");
  
  errorEstimate.resize(result.size());
  probability.resize(result.size());

  return divonne_result<Callable>(F,_hyperCube,result.size(),result.data(),errorEstimate.data(),probability.data(),stepsEvaluated);
 }

}; //end of class Cuba



} //end of namespace cuba
} //end of namespace celerium

#endif
