AM_CPPFLAGS = -I$(top_srcdir)/headers
noinst_PROGRAMS = test_wfsolver.x
test_wfsolver_x_SOURCES = test_wfsolver.cpp
test_wfsolver_x_LDADD = $(top_srcdir)/lib/WFSolver/libwfsolver.a -lgsl -lgslcblas -lm
test_wfsolver_x_CXXFLAGS = -Wall -fPIC -std=c++11
test_wfsolver_x_LDFLAGS = 

