#include "../../headers/cubawrapper.h"
#include <iostream>


int S1s2(const double x[3], double f[1]){
double zeta = 1.0;
	f[0] = 0.25*M_1_PI*pow(2*zeta,3)*0.5*exp(-2*zeta*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
	return 0;
}

int U1s(const double x[6], double f[1]){
double zeta = 1.0;
	double rr = pow(x[0]-x[3],2) + pow(x[1]-x[4],2) + pow(x[2]-x[5],2);
	if(rr < 1e-256) rr = 2e128;
	else           rr = 2*pow(rr,-0.5);

	f[0]  = pow(0.25*M_1_PI,2)*pow(2*zeta,6)*0.25;
	f[0] *= exp(-2*zeta*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));
	f[0] *= rr; 
	f[0] *= exp(-2*zeta*sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]));
	
	return 0;
}	


int main(){
  static std::array<std::pair<double,double>,3> border3d {std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0)};
  static std::array<std::pair<double,double>,6> border6d {std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0),std::make_pair(-10.0,10.0)};
  celerium::CUBA::Cuba<S1s2,3,1,&border3d> normalization(100000,10000,1e-7);
  celerium::CUBA::Cuba<U1s,6,1,&border6d> HubbardU(30000000,1000000,1e-9);

  std::array<double,1> resN, errN, pN;
  int steps = 0;
  normalization.divonne_result(resN,errN,pN,steps);
  std::cout<<"S1s2= "<<resN[0]<<" +/- "<<errN[0]<<" in "<<steps<<" steps."<<std::endl;

  std::array<double,1> resU, errU, pU;
  int stepsU = 0;
  HubbardU.suave_result(resU,errU,pU,stepsU);
  std::cout<<resU[0]<<" +/- "<<errU[0]<<" in "<<stepsU<<" steps."<<std::endl;

  return 0;
}
