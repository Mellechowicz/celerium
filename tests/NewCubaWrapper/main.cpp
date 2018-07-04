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

int SinCosSin2Cos2(const double x[4], double f[4]){
	f[0] = sin(x[0]);
	f[1] = cos(x[1]);
	f[2] = pow(sin(x[2]),2);
	f[3] = pow(cos(x[3]),2);

	return 0;
}


int main(){
  celerium::cuba::Cuba engine(3000000,100000,1e-10);

/*
  // C-style function
  std::vector<std::pair<double,double>> b3(3,std::make_pair(-10.0,10.0));
  std::array<double,1> resN, errN, pN;
  int steps = 0;
  // Suave - function
  engine.suave_result(S1s,b3,resN,errN,pN,steps);
  std::cout<<"S1s_/Suave  /= "<<resN[0]<<" +/- "<<errN[0]<<" in "<<steps<<" steps."<<std::endl;
  // Divonne - function
  engine.divonne_result(S1s,b3,resN,errN,pN,steps);
  std::cout<<"S1s_/Divonne/= "<<resN[0]<<" +/- "<<errN[0]<<" in "<<steps<<" steps."<<std::endl;

  // Callable object
  engine.parameters.maxeval = 30000000;
  std::vector<std::pair<double,double>> b6(6,std::make_pair(-10.0,10.0));
  double resU[1], errU[1], pU[1];
  // Suave
  engine.suave_result(U1s,b6,resU,errU,pU,steps);
  std::cout<<"U1s_/Suave  /= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;
  // Divonne 
  engine.divonne_result(U1s,b6,resU,errU,pU,steps);
  std::cout<<"U1s_/Divonne/= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;

  // Lambda
  engine.parameters.maxeval = 3000000;
  // Suave
  engine.suave_result([](const double x[2], double f[1]){f[0] = M_1_PI * exp(-x[0]*x[0] - x[1]*x[1]); return 0;},
		      std::vector<std::pair<double,double>>(2,std::make_pair(-10.0,10.0)),resU,errU,pU,steps);
  std::cout<<"Gaussian2d_/Suave  /= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;
  // Divonne 
  engine.divonne_result([](const double x[2], double f[1]){f[0] = M_1_PI * exp(-x[0]*x[0] - x[1]*x[1]); return 0;},
		      std::vector<std::pair<double,double>>(2,std::make_pair(-10.0,10.0)),resU,errU,pU,steps);
  std::cout<<"Gaussian2d_/Divonne/= "<<resU[0]<<" +/- "<<errU[0]<<" in "<<steps<<" steps."<<std::endl;


  // Suave nvec = 4
*/
  int stepV;
  double resV[4], errV[4], pV[4];
  engine.parameters.nvec = 4;
  engine.suave_result(SinCosSin2Cos2,std::vector<std::pair<double,double>>(1,std::make_pair(-0.0,2*M_PI)),resV,errV,pV,stepV);
  std::cout<<"Sin(x)= "<<resV[0]<<" +/- "<<errV[0]<<" in "<<stepV<<" steps."<<std::endl;
  std::cout<<"Cos(x)= "<<resV[1]<<" +/- "<<errV[1]<<" in "<<stepV<<" steps."<<std::endl;
  std::cout<<"Sin2(x)= "<<resV[2]<<" +/- "<<errV[2]<<" in "<<stepV<<" steps."<<std::endl;
  std::cout<<"Cos2(x)= "<<resV[3]<<" +/- "<<errV[3]<<" in "<<stepV<<" steps."<<std::endl;
  return 0;
}
