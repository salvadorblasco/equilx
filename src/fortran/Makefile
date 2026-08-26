# EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
# by Salvador Blasco <salvador.blasco@uv.es>

MODNAME=fsolver
TESTDIR=tests

FC = gfortran
RM = rm -f
# FFLAGS = -Wall -O2 -fPIC
FFLAGS = -g -O0 -fbacktrace -fbounds-check -ffpe-trap=zero,invalid -Wall
LIBS = -llapack -lblas


LIBEQ_SRCS = cpluss.f90 mb.f90 si.f90 nr.f90 pcf.f90 outer.f90 eq.f90 err.f90
LIBEQ_OBJS = $(LIBEQ_SRCS:.f90=.o)

LIBEMF_SRCS = emf.f90
LIBEMF_OBJS = $(LIBEMF_SRCS:.f90=.o)

.PHONY: all test 

libeq: $(LIBEQ_OBJS)
	$(FC) -shared -o libeq.so $^

pymodule: $(LIBEQ_SRCS)
	uv run f2py -c $(MODNAME).pyf -m $(MODNAME) $^ $(LIBS) 

signature:
	uv run f2py --overwrite-signature -h $(MODNAME).pyf eq.f90 

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS)

test_eq: $(LIBEQ_OBJS) $(TESTDIR)/test_cases.f90 $(TESTDIR)/test_eq.f90 
	$(FC) $(FFLAGS) $^ $(LIBS) -o $(TESTDIR)/$@

test_emf: $(LIBEMF_OBJS)  $(TESTDIR)/test_emf.f90
	$(FC) $(FFLAGS) $^ $(LIBS) -o $(TESTDIR)/$@

debug: FFLAGS = $(DEVFLAGS)
debug: $(OBJS) test_cases.f90 test.f90 
	$(FC) $(DEVFLAGS) $^ $(LIBS) -o test

clean:
	$(RM) *.o test *.mod

deepclean: clean
	$(RM) $(MODNAME)*.so
