AM_CPPFLAGS = -I$(top_srcdir)/headers
noinst_PROGRAMS = test_lattice.x
test_lattice_x_SOURCES = test_lattice.cpp ../../lib/gslmatrix.cpp ../../lib/gslmatrixcomplex.cpp ../../lib/newcubawrapper.cpp ../../lib/periodic_wannier.cpp
test_lattice_x_LDADD = -lgsl -lgslcblas -lm -lcuba
test_lattice_x_CXXFLAGS = -pedantic -Wall -std=c++17
test_lattice_x_LDFLAGS = 

