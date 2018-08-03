AM_CPPFLAGS = -I$(top_srcdir)/headers
noinst_PROGRAMS = demo_cuo2plane.x
demo_cuo2plane_x_SOURCES = demo_cuo2plane.cpp ../../lib/gslmatrix.cpp ../../lib/gslmatrixcomplex.cpp ../../lib/newcubawrapper.cpp ../../lib/periodic_wannier.cpp
demo_cuo2plane_x_LDADD = -lgsl -lgslcblas -lm -lcuba
demo_cuo2plane_x_CXXFLAGS = -pedantic -Wall -std=c++17
demo_cuo2plane_x_LDFLAGS = 

