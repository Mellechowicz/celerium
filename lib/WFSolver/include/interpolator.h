#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <vector>
#include <gsl/gsl_spline.h>


struct sample_struct {
  double x;
  double y;
};

class Interpolator {
 public:
  Interpolator(const std::vector<sample_struct> &samples,
               const gsl_interp_type * t);
  ~Interpolator();
  Interpolator(const Interpolator &interp);
  Interpolator(Interpolator &&interp);
  double operator()(double x) const;

  const std::vector<sample_struct> &GetSamples() const;
  
 private:
  gsl_interp_accel *acc; 
  gsl_spline *spline;
  double x_max;
  double x_min;
  const gsl_interp_type *t;
  std::vector<sample_struct> samples;
};

#endif /* INTERPOLATOR_H */


