#include <iostream>
#include <gsl/gsl_vector.h>
#include "../../headers/gslmatrix.h"
#include "../../headers/gslbilinear.h"

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
	celerium::gsl::Bilinear<2> blnr;
	celerium::gsl::Vector res(2);
	res(0) = 0.012311;
	res(1) = 0.11234124;
	blnr.Solve(parabola,res);

	std::cout<<res;

	res(0) = 2.000;
	res(1) = 1.111;
	auto mexicanHat = [](celerium::gsl::Vector x, size_t i){
		switch(i){
			case 0:
				return 4.0*x(0)*(x(0)*x(0) + x(1)*x(1) - 5.0);
			case 1:
				return 4.0*x(1)*(x(0)*x(0) + x(1)*x(1) - 5.0);
			}
		return 0.0;
	};
	blnr.Solve(mexicanHat,res);

	std::cout<<res<<res(0)*res(0)+res(1)*res(1)<<" should be 5"<<std::endl;


	return 0;
}
