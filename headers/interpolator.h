// This is a wrapper for the GSL spline interpolation package.
// Discrete data are provided using a vector of sample_struct.
// the Interpolator object performs interpolation and
// provides acces to the interpolated function value through
// operator(), and its two first derivatives: D1 and D2.
// Also, the interpolated function may be integrated using
// the Int method.

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <vector>
#include <gsl/gsl_spline.h>
#include <stdexcept>

namespace celerium {

namespace gsl {


// Stores a single sample:
// x - argument
// y - function value
struct sample_struct {
  double x;
  double y;

  bool operator==(sample_struct rhs) const {
    return (x == rhs.x) && (y == rhs.y);
  }
};

// Performs interpolation and
// provides interface for evaluating
// interpolated function, integration, and
// computing derivatives.
class Interpolator {
 public:

  // Initalizes empty interpolator. Note that and attempt to call
  // operator() without providing samples first
  // will result in calling an exception.
  Interpolator();

  // Initializes the interpolator with the provided samples
  // using defult gsl_interp_cspline interpolation.
  Interpolator(const std::vector<sample_struct> &samples);
  
  // Initializes the interpolator with the provided samples
  // using inteprolation type selected by the user.
  // Valid choices are:
  // gsl_interp_linear,
  // gsl_interp_polynomial,
  // gsl_interp_cspline,
  // gsl_interp_cspline_periodic,
  // gsl_interp_akima,
  // gsl_interp_akima_periodic,
  // gsl_interp_steffen.
  // See
  // https://www.gnu.org/software/gsl/manual/html_node/1D-Interpolation-Types.html#g_t1D-Interpolation-Types
  // for details.
  Interpolator(const std::vector<sample_struct> &samples,
               const gsl_interp_type * t);
  
  ~Interpolator();
  
  Interpolator(const Interpolator &interp);
  
  Interpolator(Interpolator &&interp);

  Interpolator &operator=(const Interpolator &rhs);

  // Resets the Interpolator object and initalizes in with
  // the provided samples and default gsl_interp_cspline interpolation.
  void SetSamples(const std::vector<sample_struct> &samples);

  // Resets the Interpolator object and initalizes in with
  // the provided samples and user-defined inteprolation type.  
  void SetSamples(const std::vector<sample_struct> &samples,
                  const gsl_interp_type * t);

  // Resets the interpolator by removing samples and
  // interpolation data.
  void Clear();

  // Returns the value of interpolated function.
  // The function outside the interplation region
  // (determined by the extreme x-values of the provided samples)
  // are defined to be constant and equal to their
  // respective values at the boundary of the interpolation
  // interval. This implies that the first- and second-
  // derivatives become zero as well in this regime.
  double operator()(double x) const;

  // Returns the first derivative of interpolated function.
  double D1(double x) const;

  // Returns the second derivative of interpolated function.
  double D2(double x) const;
  
  // Returns integral of the interpolated function over the interval (a, b).
  double Int(double a, double b) const;

  // Returns reference to provided samples.
  const std::vector<sample_struct> &GetSamples() const;
  
 private:
  gsl_interp_accel *acc; 
  gsl_spline *spline;
  double x_max;
  double x_min;
  const gsl_interp_type *t;
  std::vector<sample_struct> samples;
};

// Include sources.
#include "../lib/interpolator.cpp"


}  // end namespace gsl

} // end namespace celerium;

#endif /* INTERPOLATOR_H */



