SRCDIR   := ../src/equilx
BUILDDIR := ../build
TESTDIR  := .

FC = gfortran
FFLAGS = -Wall -O2 -fPIC
FFLAGS += -J$(BUILDDIR)

SRC_FILES := emf.f90
SRC_OBJS  := $(addprefix $(BUILDDIR)/,$(SRC_FILES:.f90=.o))

TEST_SRC  := test_emf.f90
TEST_OBJS := $(addprefix $(BUILDDIR)/,$(TEST_SRC:.f90=.o))

.PHONY: all clean
all: test_emf

test_emf: $(SRC_OBJS) $(TEST_OBJS)
	$(FC) $(FFLAGS) $^ -o $(TESTDIR)/$@

# Static pattern rules – now safe because BUILDDIR has no “..”
$(SRC_OBJS): $(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_OBJS): $(BUILDDIR)/%.o: $(TESTDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDDIR):
	mkdir -p $@

clean:
	rm -rf $(BUILDDIR) test_emf
