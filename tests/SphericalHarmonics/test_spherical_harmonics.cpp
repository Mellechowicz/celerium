#include "../../headers/spherical_harmonics.h"
#include <iostream>
#include <iomanip>

using namespace celerium;


int main(int argc, char *argv[])
{
  std::cout << std::fixed << std::setprecision(12);
  
  double x = 0.236;
  double y = 0.455;
  double z = 0.638;

  double r = sqrt(x*x + y*y + z*z);

  std::cout << "Compare values of spherical harmonics for selected coordinates"
            << "\n\n";
  
  std::cout << std::setw(30) << std::left << "harmonic "
            << std::setw(20) << std::left<< "caluclated"
            << std::setw(20) << std::left<< "exact";
  std::cout << "\n";

  std::cout << std::setw(30) << std::left << "l = 1,  m = 0 (aka p_z): "
            << std::setw(20) << std::left << RealSphericalHarmonic(1,0,x,y,z)
            << std::setw(20) << std::left << sqrt(3.0/4.0/M_PI)*z/r;

  std::cout << "\n";
  
  std::cout << std::setw(30) << std::left << "l = 3,  m = -2 (aka f_xyz): "
            << std::setw(20) << std::left << RealSphericalHarmonic(3,-2,x,y,z)
            << std::setw(20) << std::left << 0.5*sqrt(105.0/M_PI)*x*y*z/r/r/r;

  
  std::cout << "\n\nTEST FINISHED!\n";
  return 0;
}
