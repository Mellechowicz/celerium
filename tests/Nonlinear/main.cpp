#include <iostream>
#include <random>
#include <gsl/gsl_vector.h>
#include "../../headers/gslmatrix.h"
#include "../../headers/gslnonlinear.h"
#define DISK 5.0

double parabola(celerium::gsl::Vector& x, size_t i){
	switch(i){
		case 0:
			return 2.0*(x(0) - 1137.0);
		case 1:
			return 2.0*(x(1) - 1137.0);
	}
	return 0;
}

int main(){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> angle(0.0,2.0*M_PI);
	std::extreme_value_distribution<> radius(DISK,0.5);
	double th,r;

	celerium::gsl::Nonlinear<2> solver;
	celerium::gsl::Vector res(2);
	res(0) = 0.012311;
	res(1) = 0.11234124;

	std::cout<<"Paraboloid:"<<std::endl;
	solver.Solve(parabola,res);

	std::cout<<res;

	r  = radius(gen);
	th =  angle(gen);
	res(0) = r*sin(th);
	res(1) = r*cos(th);
	auto mexicanHat = [](celerium::gsl::Vector& x, size_t i){
		switch(i){
			case 0:
				return 4.0*x(0)*(x(0)*x(0) + x(1)*x(1) - DISK);
			case 1:
				return 4.0*x(1)*(x(0)*x(0) + x(1)*x(1) - DISK);			}
		return 0.0;
	};
	std::cout<<"Mexican hat:"<<std::endl;
	std::cout<<"Hybrids"<<std::endl;
	r  = radius(gen);
	th =  angle(gen);
	res(0) = r*sin(th);
	res(1) = r*cos(th);
	solver.Solve(mexicanHat,res,celerium::gsl::Hybrids);
	std::cout<<"error: "<<fabs(DISK-res(0)*res(0)-res(1)*res(1))<<std::endl;
	std::cout<<"Hybrid"<<std::endl;
	r  = radius(gen);
	th =  angle(gen);
	res(0) = r*sin(th);
	res(1) = r*cos(th);
	solver.Solve(mexicanHat,res,celerium::gsl::Hybrid);
	std::cout<<"error: "<<fabs(DISK-res(0)*res(0)-res(1)*res(1))<<std::endl;
	std::cout<<"Newton"<<std::endl;
	r  = radius(gen);
	th =  angle(gen);
	res(0) = r*sin(th);
	res(1) = r*cos(th);
	solver.Solve(mexicanHat,res,celerium::gsl::Newton);
	std::cout<<"error: "<<fabs(DISK-res(0)*res(0)-res(1)*res(1))<<std::endl;
	std::cout<<"Broyden"<<std::endl;
	r  = radius(gen);
	th =  angle(gen);
	res(0) = r*sin(th);
	res(1) = r*cos(th);
	solver.Solve(mexicanHat,res,celerium::gsl::Broyden);
	std::cout<<"error: "<<fabs(DISK-res(0)*res(0)-res(1)*res(1))<<std::endl;


	return 0;
}
