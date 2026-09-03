# EQ FORTRAN LIB FOR EQUILIBRIA RESOLUTION
# by Salvador Blasco <salvador.blasco@uv.es>

MODNAME=fsolver
TESTDIR=tests

FC = gfortran
RM = rm -f
FFLAGS = -Wall -O2 -fPIC
LIBFLAGS = -Wall -O2 -fPIC
TESTFLAGS = -g -O0 -fPIC -fbacktrace -fbounds-check -ffpe-trap=zero,invalid -Wall
LDFLAGS = -shared
RPATH   = -Wl,-rpath,'$$ORIGIN/../lib'
LLAPACK = -llapack -lblas

SRCDIR  = src
TESTDIR = tests
BUILD   = build
LIBDIR  = $(BUILD)/lib
BINDIR  = $(BUILD)/bin
MODDIR  = $(BUILD)/include

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

#--------------------------------------------------------------------------------
# libeq rules
#

libeq: $(EQ_LIB)

$(EQ_LIB): $(EQ_OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LLAPACK)

$(BUILD)/libeq/%.o: $(SRCDIR)/libeq/%.f90
	$(FC) $(LIBFLAGS) -J$(MODDIR) -c $< -o $@


#--------------------------------------------------------------------------------
# libemf rules
#

libemf: $(EMF_LIB)

$(EMF_LIB): $(EMF_OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LLAPACK)

$(BUILD)/libemf/%.o: $(SRCDIR)/libemf/%.f90
	$(FC) $(LIBFLAGS) -J$(MODDIR) -c $< -o $@

#--------------------------------------------------------------------------------
# tests rules


tests: $(TEST_EQ) $(TEST_EMF)

teq: $(TEST_EQ)

temf: $(TEST_EMF)

$(TEST_EQ): $(TESTDIR)/test_cases.f90 $(TESTDIR)/test_eq.f90  $(EQ_LIB)
	$(FC) $(TESTFLAGS) -J$(MODDIR)  -I$(MODDIR) -o $@ $< \
	      -L$(LIBDIR) -leq $(RPATH)

$(TEST_EMF): $(TESTDIR)/test_emf.f90 $(EMF_LIB)
	$(FC) $(TESTFLAGS) -I$(MODDIR) -o $@ $< \
	      -L$(LIBDIR) -lemf -leq $(RPATH)

run-tests: tests
	@echo "Running tests..."
	@echo ">>> test_eq"
	@$(TEST_EQ) || exit 1
	@echo ">>> test_emf"
	@$(TEST_EMF) || exit 1
	@echo "All tests passed."


#--------------------------------------------------------------------------------
# other rules

# pymodule: $(LIBEQ_SRCS)
# 	uv run f2py -c $(MODNAME).pyf -m $(MODNAME) $^ $(LIBS) 
# 
# signature:
# 	uv run f2py --overwrite-signature -h $(MODNAME).pyf eq.f90 
# 
# %.o: %.f90
# 	$(FC) -c $< -o $@ $(FFLAGS)
# 
# test_eq: $(LIBEQ_OBJS) $(TESTDIR)/test_cases.f90 $(TESTDIR)/test_eq.f90 
# 	$(FC) $(FFLAGS) $^ $(LIBS) -o $(TESTDIR)/$@
# 
# test_emf: $(LIBEMF_OBJS) $(TESTDIR)/test_emf.f90 
# 	$(FC) $(FFLAGS) $^ $(LIBS) -L. -leq -o $(TESTDIR)/$@
# 
# debug: FFLAGS = $(DEVFLAGS)
# debug: $(OBJS) test_cases.f90 test.f90 
# 	$(FC) $(DEVFLAGS) $^ $(LIBS) -o test
# 
# clean:
# 	$(RM) *.o test *.mod
# 
# deepclean: clean
# 	$(RM) $(MODNAME)*.so
