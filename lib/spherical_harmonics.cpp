#include "../headers/spherical_harmonics.h"

std::complex<double> SphericalHarmonic(int l, int m,
                                       double theta, double phi) {
  
  const std::complex<double> I {0, 1.0};
  return exp( I*phi*((double)m) ) *
      gsl_sf_legendre_sphPlm (l, abs(m), cos(theta));
}

double RealSphericalHarmonic(int l, int m,
                             double theta, double phi) {
  if (m > 0) {
    return std::pow(-1, m) * sqrt(2.0) *
        gsl_sf_legendre_sphPlm (l, abs(m), cos(theta)) *
        cos(phi*m);
  }
  else if (m < 0) {
    return std::pow(-1, m) * sqrt(2.0) *
        gsl_sf_legendre_sphPlm (l, abs(m), cos(theta)) *
        sin(-phi*m);
  }
  else return gsl_sf_legendre_sphPlm (l, abs(m), cos(theta));
}

double RealSphericalHarmonic(int l, int m,
                             double x,
                             double y,
                             double z) {

  double cos_theta;
  std::complex<double> exp_i_phi {x, y};

  const double r = sqrt(x*x + y*y + z*z);
  const double r_in_plane = sqrt(x*x + y*y);
  

  if (r == 0.0) {
    return 0.0;
  }
  else if (r_in_plane == 0) {
    cos_theta = z / r;
    exp_i_phi = 1.0;
  }
  else {
    cos_theta = z / r;
    exp_i_phi /= r_in_plane;
  }
  
  exp_i_phi = std::pow(exp_i_phi, abs(m));
  
  if (m > 0) {
    return std::pow(-1, m) * sqrt(2.0) *
        gsl_sf_legendre_sphPlm (l, abs(m), cos_theta) *
        real(exp_i_phi);
  }
  else if (m < 0) {
    return std::pow(-1, m) * sqrt(2.0) *
        gsl_sf_legendre_sphPlm (l, abs(m), cos_theta) *
        imag(exp_i_phi);
  }
  else return gsl_sf_legendre_sphPlm (l, abs(m), cos_theta);
}
