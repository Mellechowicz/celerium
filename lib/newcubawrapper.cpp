#include "../headers/newcubawrapper.h"

celerium::cuba::Cuba::Cuba(int _maxeval, int _mineval, double _epsrel):
	parameters(1, _maxeval,_mineval,_epsrel){ 
		seed_generator = std::mt19937(std::random_device()());
		uniform_distribution = std::uniform_int_distribution<>(0,std::numeric_limits<int>::max());
		spin = NULL;

		parameters.epsabs = _epsrel/64;

	}

celerium::cuba::Cuba::~Cuba(){
}

int celerium::cuba::Cuba::suave_explicit(int (*function)(const int*, const double[], const int*, double[], void*), void *userdata,
		double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
	parameters.seed = uniform_distribution(seed_generator);
	int nregions;

	Suave(parameters.nDim, parameters.nComp, function, userdata, parameters.nvec,
			parameters.epsrel, parameters.epsabs, parameters.verbose | parameters.last | parameters.level, parameters.seed,
			parameters.mineval, parameters.maxeval, parameters.nnew, parameters.nmin, parameters.flatness,
			parameters.statefile, &spin, &nregions, &stepsEvaluated, 
			&hasFailed, result, errorEstimate, probability);

	for(int i=0; i<(parameters.nComp*parameters.nvec); ++i){
		result[i]        *= jacobian;
		errorEstimate[i] *= jacobian;
	}

	return hasFailed;
}

int celerium::cuba::Cuba::divonne_explicit(int (*function)(const int*, const double[], const int*, double[], void*), void *userdata,
		double result[], double errorEstimate[], double probability[], int& stepsEvaluated){
	parameters.seed = uniform_distribution(seed_generator);
	int nregions;

	Divonne(parameters.nDim, parameters.nComp, function, userdata, parameters.nvec,
			parameters.epsrel, parameters.epsabs, parameters.verbose | parameters.last | parameters.level, parameters.seed,
			parameters.mineval, parameters.maxeval, parameters.key1, parameters.key2, parameters.key3,
			parameters.maxpass, parameters.border, parameters.maxchisq, parameters.mindeviation,
			parameters.ngiven, parameters.ldxgiven, parameters.xgiven, parameters.nextra, NULL, 
			parameters.statefile, &spin, &nregions, &stepsEvaluated, 
			&hasFailed, result, errorEstimate, probability);

	for(int i=0; i<parameters.nComp; ++i){
		result[i]        *= jacobian;
		errorEstimate[i] *= jacobian;
	}

	return hasFailed;
}

