#include "../../headers/newcubawrapper.h"
#include <iostream>
#include <cmath>
#include <algorithm>


int S1s(const double x[3], double f[1]){
double zeta = 1.0;
	f[0] = 0.25*M_1_PI*pow(2*zeta,3)*0.5*exp(-2*zeta*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
	return 0;
}

int U1s(const double x[6], double f[1]){
double zeta = 1.0;
	double rr = pow(x[0]-x[3],2) + pow(x[1]-x[4],2) + pow(x[2]-x[5],2);
	if(rr < 1e-64) rr = 2e32;
	else           rr = 2*pow(rr,-0.5);

	f[0]  = pow(0.25*M_1_PI,2)*pow(2*zeta,6)*0.25;
	f[0] *= exp(-2*zeta*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
	f[0] *= rr; 
	f[0] *= exp(-2*zeta*sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]));
	
	return 0;
}	


int main(){
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10.0,10.0));
  std::vector<std::pair<double,double>> b6(6,std::make_pair(-10.0,10.0));

  celerium::cuba::Cuba engine(b6,30000000,10000000,1e-10);

  std::array<double,1> resN, errN, pN;
  int steps = 0;
  engine.suave_result(S1s,b3,resN,errN,pN,steps);
  std::cout<<"S1s= "<<resN[0]<<" +/- "<<errN[0]<<" in "<<steps<<" steps."<<std::endl;


  double resU[1], errU[1], pU[1];
  int stepsU= 0;
  engine.suave_result(U1s,b6,resU,errU,pU,stepsU);
  std::cout<<"U1s= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<stepsU<<" steps."<<std::endl;

  engine.suave_result([](const double x[1], double f[1]){f[0] = sin(x[0]); return 0;},std::vector<std::pair<double,double>>(1,std::make_pair(0.0,M_PI_2)),resU,errU,pU,stepsU);
  return 0;
}
