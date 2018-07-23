#ifndef INTERPOLATOR3D_H
#define INTERPOLATOR3D_H

#include <vector>
#include <array>
#include <cmath>

namespace celerium {

class Interpolator3d {

 public:
  
  Interpolator3d() {};

  void Initialize(
      const std::array<std::pair<double, double>, 3> bounding_box,
      const std::array<int, 3> &cube_dimensions,
      const std::vector<double> &samples) {
    this->cube_dimensions = cube_dimensions;
    this->bounding_box = bounding_box;
    this->samples = samples;

    for (int i = 0; i < 3; ++i)
      this->box_dimensions[i] = bounding_box[i].second - bounding_box[i].first;
    
  }

  double operator()(const double xx[]) {
    
    for (int i = 0; i < 3; ++i)
      if ( xx[i] < bounding_box[i].first ||
           xx[i] > bounding_box[i].second ) return 0.0;

    int index[3];
    double dx[3];

    double x_rel;
    for (int i = 0; i < 3; ++i) {
      x_rel = (xx[i]+bounding_box[i].first)/box_dimensions[i]*cube_dimensions[i];
    index[i] = std::floor(x_rel);
    dx[i] = std::fmod(x_rel, 1.0);
    }
    /*
    for (int  i = 0; i < 3; ++i) 
      if (index[i] == cube_dimensions[i] - 1) {
        index[i]--;
        dx[i]=1;
      }
    */
    const int n0 = cube_dimensions[1]*cube_dimensions[2];
    const int n1 = cube_dimensions[2];
    
    const double c000 = samples[index[0]*n0 + index[1]*n1 + index[2]];
    const double c100 = samples[(index[0]+1)*n0 + index[1]*n1 + index[2]];
    const double c010 = samples[index[0]*n0 + (index[1]+1)*n1 + index[2]];
    const double c001 = samples[index[0]*n0 + index[1]*n1 + index[2] + 1];
    const double c110 = samples[(index[0]+1)*n0 + (index[1]+1)*n1 + index[2]];
    const double c011 = samples[index[0]*n0 + (index[1]+1)*n1 + index[2] + 1];
    const double c101 = samples[(index[0]+1)*n0 + index[1]*n1 + index[2] + 1];
    const double c111 = samples[(index[0]+1)*n0 + (index[1]+1)*n1 + index[2] + 1];

    const double c00 = c000*(1-dx[0]) + c100*dx[0];
    const double c01 = c001*(1-dx[0]) + c101*dx[0];
    const double c10 = c010*(1-dx[0]) + c110*dx[0];
    const double c11 = c011*(1-dx[0]) + c111*dx[0];

    const double c0 = c00*(1-dx[1]) + c10*dx[1];
    const double c1 = c01*(1-dx[1]) + c11*dx[1];

    return c0*(1-dx[2]) + c1*dx[2];
  }
  
  int SaveToFile(const char* file_name);
  int LoadFromFile(const char* file_name);
  

 private:
  std::array<int, 3> cube_dimensions;
  std::array<std::pair<double, double>, 3> bounding_box;
  std::array<double, 3> box_dimensions;
  std::vector<double> samples;
};


}  // celerium

#endif /* INTERPOLATOR3D_H */
