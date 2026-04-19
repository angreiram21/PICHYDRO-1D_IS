FC := gfortran
FFLAGS ?= -ffixed-line-length-none -mcmodel=large

BUILDDIR := build
MODDIR := $(BUILDDIR)
TARGET := $(BUILDDIR)/hydro

SRC := src/ic_smash.f90 src/fluidization_smash.f90 src/hydro_1D_IS.f src/cornelius.f90
OBJ := $(BUILDDIR)/ic_smash.o \
       $(BUILDDIR)/fluidization_smash.o \
       $(BUILDDIR)/hydro_1D_IS.o \
       $(BUILDDIR)/cornelius.o

.PHONY: all clean

all: $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/ic_smash.o: src/ic_smash.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(BUILDDIR)/fluidization_smash.o: src/fluidization_smash.f90 $(BUILDDIR)/ic_smash.o
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(BUILDDIR)/hydro_1D_IS.o: src/hydro_1D_IS.f $(BUILDDIR)/ic_smash.o $(BUILDDIR)/fluidization_smash.o
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(BUILDDIR)/cornelius.o: src/cornelius.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) -c $< -o $@

$(TARGET): $(OBJ)
	$(FC) -mcmodel=large $(OBJ) -o $(TARGET)

clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod $(TARGET)
