#include <interpolator3d.h>

#include <iostream>

using namespace celerium;

int main(int argc, char *argv[])
{

  std::vector<double> data;

  for (int i  = 0; i < 2; ++i) {
    for (int j  = 0; j < 2; ++j) {
      for (int k  = 0; k < 2; ++k) {
        data.push_back(i + j);
      }
    }
  }


    for (int i  = 0; i < 2; ++i) {
    for (int j  = 0; j < 2; ++j) {
          std::cout << i+j << " ";
    }
    std::cout << "\n";
  }

  
  
  Interpolator3d interp3d;

  interp3d.Initialize({{{0, 1}, {0, 1}, {0, 1}}}, {{2, 2, 2}}, data);

  std::cout << "\n\n";

  for (int i  = 0; i < 3; ++i) {
    for (int j  = 0; j < 3; ++j) {
          double xx [] = {i*1.0/2, j*1.0/2, 0.1};
          std::cout << interp3d(xx) << " ";
    }
    std::cout << "\n";
  }

  std::cout << "\n\n";
  
  return 0;
}
