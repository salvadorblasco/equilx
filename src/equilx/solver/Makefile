# EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
# by Salvador Blasco <salvador.blasco@uv.es>

MODNAME=fsolver

FC = gfortran
RM = rm -f
FFLAGS = -Wall -O2
LIBS = -llapack -lblas
DEVFLAGS = -g -O0 -fbacktrace -fbounds-check -ffpe-trap=zero,invalid -Wall

SRCS = cpluss.f90 mb.f90 si.f90 nr.f90 pcf.f90 outer.f90 eq.f90
OBJS = $(SRCS:.f90=.o)

.PHONY: all test debug fsolver clean

module: $(SRCS)
	uv run f2py -c $(MODNAME).pyf -m $(MODNAME) $^ $(LIBS) 

signature:
	uv run f2py -h $(MODNAME).pyf eq.f90 

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS)

test: $(OBJS) test_cases.f90 test.f90 
	$(FC) $(FFLAGS) $^ $(LIBS) -o $@

debug: FFLAGS = $(DEVFLAGS)
debug: $(OBJS) test_cases.f90 test.f90 
	$(FC) $(DEVFLAGS) $^ $(LIBS) -o test


clean:
	$(RM) *.o test *.mod

deepclean: clean
	$(RM) $(MODNAME)*.so
