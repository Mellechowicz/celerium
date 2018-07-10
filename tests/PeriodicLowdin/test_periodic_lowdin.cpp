#include <periodic_lowdin.h>

using namespace celerium;

int main(int argc, char *argv[])
{

  gsl::Matrix overlap_000(1, {1.0});
  gsl::Matrix overlap_100(1, {0.1});


  std::vector<std::vector<double>> result;

  std::vector<std::pair<std::array<int, 3>, gsl::Matrix>> real_space_overlaps;

  real_space_overlaps.push_back({ {0, 0, 0}, overlap_000 });
  real_space_overlaps.push_back({ {1, 0, 0}, overlap_100 });

  PeriodicLowdin(real_space_overlaps,
                 {4, 1, 1},
                 result);

  for (auto m : result) {
    std::cout << m[0] << "\n";  
  }

  
  
  
  
  return 0;
}
