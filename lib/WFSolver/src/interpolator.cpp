#include "../include/interpolator.h"


Interpolator::Interpolator(const std::vector<sample_struct> &samples,
                           const gsl_interp_type * t) {
  
  double x [samples.size()];
  double y [samples.size()];

  for (size_t i = 0; i < samples.size(); ++i) {
    x[i] = samples[i].x;
    y[i] = samples[i].y;
    if (x[i] > x_max || i == 0) x_max = x[i];
    if (x[i] < x_min || i == 0) x_min = x[i];
  }

  this->acc = gsl_interp_accel_alloc ();
  this->spline = gsl_spline_alloc (t, samples.size());
  gsl_spline_init (spline, x, y, samples.size());
  this->t = t;
  this->samples = samples;
}

Interpolator::Interpolator(const Interpolator &interp) {

  this->samples=interp.samples;
  this->t = interp.t;
  this->x_max = interp.x_max;
  this->x_min = interp.x_min;
  
  double x [this->samples.size()];
  double y [this->samples.size()];
  
  for (size_t i = 0; i < this->samples.size(); ++i) {
    x[i] = this->samples[i].x;
    y[i] = this->samples[i].y;
  }
  
  this->acc = gsl_interp_accel_alloc ();
  this->spline = gsl_spline_alloc (this->t, this->samples.size());
  gsl_spline_init (this->spline, x, y, this->samples.size());
}


Interpolator::Interpolator(Interpolator &&interp) {
  
  this->acc = interp.acc;
  this->spline = interp.spline;
  this->x_max = interp.x_max;
  this->x_min = interp.x_min;
  this->t = interp.t;
  this->samples = interp.samples;

  interp.acc = nullptr;
  interp.spline = nullptr;
}

Interpolator::~Interpolator() {
  gsl_spline_free (this->spline);
  gsl_interp_accel_free (this->acc);
}

double Interpolator::operator()(double x) const {
  if (x > this->x_max) {
    return gsl_spline_eval (this->spline, x_max, this->acc);
  }
  else if (x < this->x_min) {
    return gsl_spline_eval (this->spline, x_min, this->acc);
  }  
  return gsl_spline_eval (this->spline, x, this->acc);
}

const std::vector<sample_struct> &Interpolator::GetSamples() const {
  return this->samples;
}
