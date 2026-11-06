# EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
# by Salvador Blasco <salvador.blasco@uv.es>

TESTDIR=tests

FC = gfortran
F2PY = uv run f2py
RM = rm 

DEBUG ?= 0
ifeq ($(BUILD),1)
	FFLAGS = -g -O0 -fPIC -fbacktrace -fbounds-check -ffpe-trap=zero,invalid -Wall
else
	FFLAGS = -Wall -O2 -fPIC
end if

LDFLAGS = -shared
RPATH   = -Wl,-rpath,'$$ORIGIN/../lib'
LLAPACK = -llapack -lblas

SRCDIR  = src
TESTDIR = tests
BUILD   = build
LIBDIR  = $(BUILD)/lib
BINDIR  = $(BUILD)/bin
MODDIR  = $(BUILD)/include
PYMDIR  = $(SRCDIR)/equilx

# cpluss.f90 mb.f90 si.f90 nr.f90 pcf.f90 outer.f90 eq.f90 err.f90
EQ_SRCS = $(wildcard $(SRCDIR)/libeq/*.f90)
EQ_OBJS = $(patsubst $(SRCDIR)/libeq/%.f90,$(BUILD)/libeq/%.o,$(EQ_SRCS))
EQ_LIB   = $(LIBDIR)/libeq.so

EMF_SRCS = $(wildcard $(SRCDIR)/libemf/*.f90)
EMF_OBJS = $(patsubst $(SRCDIR)/libemf/%.f90,$(BUILD)/libemf/%.o,$(EMF_SRCS))
EMF_LIB   = $(LIBDIR)/libemf.so

TEST_EQ   = $(BINDIR)/test_eq
TEST_EMF  = $(BINDIR)/test_emf

.PHONY: all test 


vars:
	@echo $(EQ_SRCS) 
	@echo $(EQ_OBJS) 


#--------------------------------------------------------------------------------
# common rules
#
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ 


#--------------------------------------------------------------------------------
# python module rule
#

PYMDIR = $(SRCDIR)/equilx
MODNAME=fsolver
SIGNATURE=$(PYMDIR)/$(MODNAME).pyf

pymodule: $(EQ_SRCS) $(EMF_SRCS)
	$(F2PY) -c $(SIGNATURE) -m $(BINDIR)/$(MODNAME) $(EQ_SRCS) $^ $(LLAPACK)
	mv $(MODNAME)*.so $(PYMDIR)

signature:
	$(F2PU) --overwrite-signature -h $(SIGNATURE) eq.f90 


#--------------------------------------------------------------------------------
# libeq rules
#

libeq: $(EQ_LIB)

$(EQ_LIB): $(EQ_OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LLAPACK)

$(BUILD)/libeq/%.o: $(SRCDIR)/libeq/%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

#$(EQ_OBJS): %.o: %.f90
#	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@


#--------------------------------------------------------------------------------
# libemf rules
#

libemf: $(EMF_LIB)

$(EMF_LIB): $(EMF_OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LLAPACK)

$(BUILD)/libemf/%.o: $(SRCDIR)/libemf/%.f90
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

#--------------------------------------------------------------------------------
# tests rules

TESTCASES:=$(wildcard $(TESTDIR)/case*.f90)


tests: $(TEST_EQ) $(TEST_EMF)

teq: $(TEST_EQ)

temf: $(TEST_EMF)

$(TEST_EQ):  $(TESTDIR)/test_eq.f90  $(EQ_LIB)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $(TESTDIR)/test_cases.f90 -o $(TESTDIR)/test_cases.o
	$(FC) $(FFLAGS) -J$(MODDIR)  -I$(MODDIR) -o $@ $(TESTDIR)/test_cases.o $< -L$(LIBDIR) -leq $(RPATH)

$(TEST_EMF): $(TESTDIR)/test_emf.f90 $(EMF_LIB)
	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR)  -o $@ $< \
	      -L$(LIBDIR) -lemf -leq $(RPATH)

# teqd: $(TESTDIR)/test_eq.f90 $(TESTDIR)/test_cases.f90 $(EQ_OBJS)
# 	$(FC) $(FFLAGS) -J$(MODDIR) -c $(TESTDIR)/test_cases.f90 -o $(TESTDIR)/test_cases.o
# 	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR)  -o $(TEST_EQ) $(TESTCASES) $(TESTDIR)/test_cases.o $(EQ_OBJS) $< $(LLAPACK)
# 
# temfd: $(TESTDIR)/test_emf.f90 $(EMF_OBJS) $(EQ_OBJS)
# 	$(FC) $(FFLAGS) -I$(MODDIR) -J$(MODDIR)  -o $(TEST_EMF) $(TESTCASES) $<  $(EMF_OBJS) $(EQ_OBJS) $(LLAPACK)

run-tests: tests
	@echo "Running tests..."
	@echo ">>> test_eq"
	@$(TEST_EQ) || exit 1
	@echo ">>> test_emf"
	@$(TEST_EMF) || exit 1
	@echo "All tests passed."


#--------------------------------------------------------------------------------
# other rules

clean:
	$(RM) $(EQ_OBJS) $(EMF_OBJS) $(MODDIR)/*
