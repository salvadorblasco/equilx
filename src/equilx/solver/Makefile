SH=/bin/bash
.PHONY=fortran
FC = gfortran
FFLAGS = -Wall -O2

fsolver: pcf.o si.o cpluss.o mb.o nr.o
	uv run f2py -llapack -c -m fsolver pcf.o si.o cpluss.o mb.o nr.o eq.f90

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS)

test: test.f90 cpluss.o mb.o si.o
	$(FC) test.f90 cpluss.o mb.o si.o -o test

clean:
	rm *.o *.so

# clibeq: gcc -Wall -shared eq.c
