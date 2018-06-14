#include "../../headers/interpolator.h"
#include <iostream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <cmath>


using namespace celerium;

int main()
{

  std::cout << "\nTest of the Interpolator class for sine function. \n\n";
  
  // Set displayed number of digits.
  std::cout << std::fixed << std::setprecision(6);

  // Coarsly sample sine function on (0, PI) interval.
  std::vector<sample_struct> samples_sin {
    {0.0, sin(0.0)},
    {0.1, sin(0.1)},
    {0.11, sin(0.11)},
    {1.14, sin(1.14)},
    {1.7, sin(1.7)},
    {2.4, sin(2.4)},
    {M_PI, sin(M_PI)},
  };

  // Coarsly sample cos function on (0, PI) interval.
  std::vector<sample_struct> samples_cos {
    {0.0, cos(0.0)},
    {0.1, cos(0.1)},
    {0.11, cos(0.11)},
    {1.14, cos(1.14)},
    {1.7, cos(1.7)},
    {2.4, cos(2.4)},
    {M_PI, cos(M_PI)},
  };


  // Interpolate between samples using default settings (csplines).
  Interpolator interp(samples_sin);

  std::cout << std::setw(20) << "x";
  std::cout << std::setw(20) << "sin(x)";
  std::cout << std::setw(20) << "interp(x)";
  std::cout << std::setw(20) << "d/dx interp(x)";
  std::cout << std::setw(20) << "d2/dx2 interp(x)";
  std::cout << "\n";
  
  for (size_t i = 0; i <= 20; ++i) {
    std::cout << std::setw(20) << i*M_PI/20;
    std::cout << std::setw(20) << sin(i*M_PI/20);
    std::cout << std::setw(20) << interp(i*M_PI/20);
    std::cout << std::setw(20) << interp.D1(i*M_PI/20);
    std::cout << std::setw(20) << interp.D2(i*M_PI/20);
    std::cout << "\n";
  }

  // Calculate integral of the interpolated function on the interval (0, PI).
  // Exact value for sine function: 2.0.
  std::cout << "\nIntegral of interp(x) over the range (0, PI): "
            << interp.Int(0.0, M_PI) << "\n\n";


  
  std::cout << "\nTest of the Interpolator class for cos function. \n\n";

  // Let us change the data.
  interp.SetSamples(samples_cos);

  std::cout << std::setw(20) << "x";
  std::cout << std::setw(20) << "cos(x)";
  std::cout << std::setw(20) << "interp(x)";
  std::cout << std::setw(20) << "d/dx interp(x)";
  std::cout << std::setw(20) << "d2/dx2 interp(x)";
  std::cout << "\n";
  
  for (size_t i = 0; i <= 20; ++i) {
    std::cout << std::setw(20) << i*M_PI/20;
    std::cout << std::setw(20) << sin(i*M_PI/20);
    std::cout << std::setw(20) << interp(i*M_PI/20);
    std::cout << std::setw(20) << interp.D1(i*M_PI/20);
    std::cout << std::setw(20) << interp.D2(i*M_PI/20);
    std::cout << "\n";
  }
  
  
  // Calculate integral of the interpolated function on the interval (0, PI).
  // Exact value for cos function: 0.0.
  std::cout << "\nIntegral of interp(x) over the range (0, PI): "
            << interp.Int(0.0, M_PI);
  
  
  std::cout << "\n\nTEST FINISHED...\n\n";

  
  return 0;
}
