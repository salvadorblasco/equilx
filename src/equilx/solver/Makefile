MODNAME=fsolver
SH=/bin/bash
.PHONY=fortran
FC = gfortran
FFLAGS = -Wall -O2
LIBS = -llapack -lblas
DEBUGFLAGFS = -g -O0 -fbacktrace -fbounds-check -ffpe-trap=zero,invalid,overflow,underflow -Wall
fsolver: pcf.o si.o cpluss.o mb.o nr.f90
	uv run f2py -c -m $@ $^

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS)

test: test.f90 cpluss.o mb.o si.o nr.o eq.o pcf.o
	$(FC) $(FFLAGS) $^ $(LIBS) -o $@

debug_test: test.f90 cpluss.f90 mb.f90 si.f90 nr.f90 eq.f90 pcf.f90
	$(FC) $(DEBUGFLAGS) $^ $(LIBS) -o test

clean:
	rm *.o *.so test

# clibeq: gcc -Wall -shared eq.c
#
# NCP
# 1 0     ^
# 0 1  NEQ+NCP
# 1 1     v
