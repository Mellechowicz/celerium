#include "../../headers/newcubawrapper.h"
#include <iostream>
#include <cmath>
#include <algorithm>


int S1s(const double x[3], double f[1]){
double zeta = 1.0;
	f[0] = 0.25*M_1_PI*pow(2*zeta,3)*0.5*exp(-2*zeta*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
	return 0;
}

struct microscopic{
int operator()(const double x[6], double f[1]){
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
} U1s;


int main(){
  celerium::cuba::Cuba engine(3000000,100000,1e-10);

  // C-style function
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10.0,10.0));
  std::array<double,1> resN, errN, pN;
  int steps = 0;
  // Suave - function
  engine.suave_result(S1s,b3,resN,errN,pN,steps);
  std::cout<<"S1s_/*Suave  */= "<<resN[0]<<" +/- "<<errN[0]<<" in "<<steps<<" steps."<<std::endl;
  // Divonne - function
  engine.divonne_result(S1s,b3,resN,errN,pN,steps);
  std::cout<<"S1s_/*Divonne*/= "<<resN[0]<<" +/- "<<errN[0]<<" in "<<steps<<" steps."<<std::endl;

  // Callable object
  engine.parameters.maxeval = 30000000;
  std::vector<std::pair<double,double>> b6(6,std::make_pair(-10.0,10.0));
  double resU[1], errU[1], pU[1];
  // Suave
  engine.suave_result(U1s,b6,resU,errU,pU,steps);
  std::cout<<"U1s_/*Suave  */= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;
  // Divonne 
  engine.divonne_result(U1s,b6,resU,errU,pU,steps);
  std::cout<<"U1s_/*Divonne*/= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;

  // Lambda
  engine.parameters.maxeval = 3000000;
  // Suave
  engine.suave_result([](const double x[2], double f[1]){f[0] = M_1_PI * exp(-x[0]*x[0] - x[1]*x[1]); return 0;},
		      std::vector<std::pair<double,double>>(2,std::make_pair(-10.0,10.0)),resU,errU,pU,steps);
  std::cout<<"Gaussian2d_/*Suave  */= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;
  // Divonne 
  engine.divonne_result([](const double x[2], double f[1]){f[0] = M_1_PI * exp(-x[0]*x[0] - x[1]*x[1]); return 0;},
		      std::vector<std::pair<double,double>>(2,std::make_pair(-10.0,10.0)),resU,errU,pU,steps);
  std::cout<<"Gaussian2d_/*Divonne*/= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;
  return 0;
}
