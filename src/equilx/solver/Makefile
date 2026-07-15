SH=/bin/bash
.PHONY=fortran


fortran:
	f2py3 -llapack -c -m mfeq mfeq.pyf eq.f90

si:
	uv run f2py -c -m si si.f90

pcf: cpluss.o
	uv run f2py -c -m pcf cpluss.o pcf.f90

mb:
	uv run f2py -c -m mb mb.f90

cpluss.o:
	gfortran -c cpluss.f90

# clibeq: gcc -Wall -shared eq.c
