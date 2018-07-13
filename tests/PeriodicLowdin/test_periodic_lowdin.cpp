#include <periodic_lowdin.h>
#include <gslmatrixcomplex.h>
#include <iomanip>
#include <algorithm>

using namespace celerium;

int main(int argc, char *argv[])
{

  std::complex<double> I {0, 1};
  
  std::cout << std::setprecision(12);

  std::vector<std::vector<std::vector<std::complex<double>>>> result;

  std::vector<std::pair<std::array<int, 3>,
                        gsl::MatrixComplex>> overlaps;



  overlaps.push_back(
      { {0, 0, 0}, gsl::MatrixComplex(2, {1.0, 0.1, 0.1, 1.0}) });

  overlaps.push_back(
      { {+1, 0, 0}, gsl::MatrixComplex(2, {0.0, 0.5, 0.2, 0.0}) });

  overlaps.push_back(
      { {-1, 0, 0}, gsl::MatrixComplex(2, {0.0, 0.2, 0.5, 0.0}) });
  
  overlaps.push_back(
      { {+2, 0, 0}, gsl::MatrixComplex(2, {0.0, 0.0, 0.0, 0.0}) });

  overlaps.push_back(
      { {-2, 0, 0}, gsl::MatrixComplex(2, {0.0, 0.0, 0.0, 0.0}) });

  overlaps.push_back(
      { {+3, 0, 0}, gsl::MatrixComplex(2, {0.0, 0.0, 0.0, 0.0}) });

  overlaps.push_back(
      { {-3, 0, 0}, gsl::MatrixComplex(2, {0.0, 0.0, 0.0, 0.0}) });


  PeriodicOthogonalization(overlaps,
                           {overlaps.size(), 1, 1},
                           {{0, 0, 0}, {1, 0, 0}},
                           result);


  std::cout << "solution: " << "\n";


  
  for (size_t i = 0; i < result[0].size(); ++i) {
    std::cout << "Wannier #" << i << ":\n";
    for (size_t j = 0; j <  result[0][i].size(); ++j) {
      auto pos = overlaps[j / overlaps[0].second.columnNumber()].first;
      std::cout << j % overlaps[0].second.columnNumber();
      std::cout << ", (" << pos[0] << ", " << pos[1] << ", " << pos[2] << "): ";
      std::cout <<  result[0][i][j] << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "\n\n";


   
  std::cout << "Wannier overlaps: \n";

  for (size_t o1 = 0; o1 < 2; ++o1) {
    for (size_t o2 = 0; o2 < 2; ++o2) {
      std::cout << ScalarProduct(overlaps, result, 0, o1, 0, o2).real() << " ";
    }
    std::cout << "\n";
  }
  


      


  return 0;
}
