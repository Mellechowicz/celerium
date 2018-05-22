AM_CPPFLAGS = -I$(top_srcdir)/headers
noinst_PROGRAMS = test_potential_and_wfsolver.x
test_potential_and_wfsolver_x_SOURCES = test_potential_and_wfsolver.cpp
test_potential_and_wfsolver_x_LDADD = $(top_srcdir)/lib/WFSolver/libwfsolver.a -lgsl -lgslcblas -lm
test_potential_and_wfsolver_x_CXXFLAGS = -Wall -fPIC -std=c++17
test_potential_and_wfsolver_x_LDFLAGS = 

