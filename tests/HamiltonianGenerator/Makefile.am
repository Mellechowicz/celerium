AM_CPPFLAGS = -I$(top_srcdir)/headers
noinst_PROGRAMS = test_hamiltonian_generator.x
test_hamiltonian_generator_x_SOURCES = test_hamiltonian_generator.cpp ../../lib/gslmatrix.cpp ../../lib/gslmatrixcomplex.cpp ../../lib/newcubawrapper.cpp ../../lib/periodic_wannier.cpp
test_hamiltonian_generator_x_LDADD = -lgsl -lgslcblas -lm -lcuba
test_hamiltonian_generator_x_CXXFLAGS = -pedantic -Wall -std=c++17
test_hamiltonian_generator_x_LDFLAGS = 

