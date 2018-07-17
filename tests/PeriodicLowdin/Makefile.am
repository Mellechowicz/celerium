AM_CPPFLAGS = -I$(top_srcdir)/headers
noinst_PROGRAMS = test_periodic_wannier.x
test_periodic_wannier_x_SOURCES = test_periodic_wannier.cpp ../../lib/gslmatrix.cpp ../../lib/gslmatrixcomplex.cpp
test_periodic_wannier_x_LDADD = -lgsl -lgslcblas -lm -lcuba
test_periodic_wannier_x_CXXFLAGS = -pedantic -Wall -std=c++17
test_periodic_wannier_x_LDFLAGS = 

