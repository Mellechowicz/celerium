#include "../headers/interpolator.h"


Interpolator::Interpolator() {
  this->acc = nullptr;
  this->spline = nullptr;
  this->t = nullptr;
  this->x_min = 0;
  this->x_max = 0;
}


Interpolator::Interpolator(const std::vector<sample_struct> &samples) {
  double x [samples.size()];
  double y [samples.size()];

  this->x_max = 0;
  this->x_min = 0;
  
  for (size_t i = 0; i < samples.size(); ++i) {
    x[i] = samples[i].x;
    y[i] = samples[i].y;
    if (x[i] > x_max || i == 0) x_max = x[i];
    if (x[i] < x_min || i == 0) x_min = x[i];
  }
    
  this->acc = gsl_interp_accel_alloc();
  this->spline = gsl_spline_alloc(gsl_interp_cspline, samples.size());
  gsl_spline_init(spline, x, y, samples.size());
  this->t = gsl_interp_cspline;
  this->samples = samples;
}


Interpolator::Interpolator(const std::vector<sample_struct> &samples,
                           const gsl_interp_type *t) {
  
  double x [samples.size()];
  double y [samples.size()];

  this->x_max = 0;
  this->x_min = 0;

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
  this->samples = std::move(interp.samples);

  interp.acc = nullptr;
  interp.spline = nullptr;
}

Interpolator::~Interpolator() {
  gsl_spline_free (this->spline);
  gsl_interp_accel_free (this->acc);
}

Interpolator &Interpolator::operator=(const Interpolator &rhs) {
  double x [rhs.samples.size()];
  double y [rhs.samples.size()];
  
  gsl_spline_free (this->spline);
  gsl_interp_accel_free (this->acc);

  this->x_max = 0;
  this->x_min = 0;

  for (size_t i = 0; i < rhs.samples.size(); ++i) {
    x[i] = rhs.samples[i].x;
    y[i] = rhs.samples[i].y;
    if (x[i] > x_max || i == 0) x_max = x[i];
    if (x[i] < x_min || i == 0) x_min = x[i];
  }

  this->acc = gsl_interp_accel_alloc();
  this->spline = gsl_spline_alloc(gsl_interp_cspline, rhs.samples.size());
  gsl_spline_init(spline, x, y, rhs.samples.size());
  this->t = rhs.t;
  this->samples = rhs.samples;
  return *this;
}



void Interpolator::SetSamples(const std::vector<sample_struct> &samples) {

  gsl_spline_free (this->spline);
  gsl_interp_accel_free (this->acc);

  double x [samples.size()];
  double y [samples.size()];

  this->x_max = 0;
  this->x_min = 0;
  
  for (size_t i = 0; i < samples.size(); ++i) {
    x[i] = samples[i].x;
    y[i] = samples[i].y;
    if (x[i] > x_max || i == 0) x_max = x[i];
    if (x[i] < x_min || i == 0) x_min = x[i];
  }

  this->acc = gsl_interp_accel_alloc ();
  this->spline = gsl_spline_alloc (gsl_interp_cspline, samples.size());
  gsl_spline_init (spline, x, y, samples.size());
  this->t = gsl_interp_cspline;
  this->samples = samples;
}

  void Interpolator::SetSamples(const std::vector<sample_struct> &samples,
                  const gsl_interp_type * t) {

  gsl_spline_free (this->spline);
  gsl_interp_accel_free (this->acc);
    
  double x [samples.size()];
  double y [samples.size()];

  this->x_max = 0;
  this->x_min = 0;

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

void Interpolator::Clear() {
  gsl_spline_free (this->spline);
  gsl_interp_accel_free (this->acc);
  this->acc = nullptr;
  this->spline = nullptr;
  this->t = nullptr;
  this->x_min = 0;
  this->x_max = 0;
  this->samples.clear();
}



double Interpolator::operator()(double x) const {

  if(this->spline == nullptr) {
    throw std::runtime_error("Interpolator is not \
initialized with any data. Interpolator::operator() cannot be executed.");
  }

  
  if (x > this->x_max) {
    return gsl_spline_eval (this->spline, x_max, this->acc);
  }
  else if (x < this->x_min) {
    return gsl_spline_eval (this->spline, x_min, this->acc);
  }  
  return gsl_spline_eval (this->spline, x, this->acc);
}

double Interpolator::D1(double x) const {

  if(this->spline == nullptr) {
    throw std::runtime_error("Interpolator is not \
initialized with any data. Interpolator::D1 cannot be executed.");
  }

  if (x > this->x_max) {
    return 0;
  }
  else if (x < this->x_min) {
    return 0;
  }
  return gsl_spline_eval_deriv(this->spline, x, this->acc);
}

double Interpolator::D2(double x) const {

  if(this->spline == nullptr) {
    throw std::runtime_error("Interpolator is not \
initialized with any data. Interpolator::D2 cannot be executed.");
  }

  
  if (x > this->x_max) {
    return 0;
  }
  else if (x < this->x_min) {
    return 0;
  }
  return gsl_spline_eval_deriv2(this->spline, x, this->acc);
}

double Interpolator::Int(double a, double b) const {

  if(this->spline == nullptr) {
    throw std::runtime_error("Interpolator is not \
initialized with any data. Interpolator::Int cannot be executed.");
  }
  
  double result = 0.0;

  double a_new = a < b ? a : b;
  double b_new = a < b ? b : a;

  result +=  gsl_spline_eval_integ(this->spline,
                                   a_new > x_min ? a_new : x_min,
                                   b_new < x_max ? b_new : x_max,
                                   this->acc);
  if (a_new < x_min) {
    result +=
        gsl_spline_eval(this->spline, this->x_min, this->acc) * (x_min - a_new);
  }

  if (b_new > x_max) {
    result +=
        gsl_spline_eval(this->spline, this->x_max, this->acc) * (b_new - x_max);
  }

  if (b < a) result *= -1;

  return result;
}





const std::vector<sample_struct> &Interpolator::GetSamples() const {
  return this->samples;
}
