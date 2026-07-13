# Build FieldsForge.
#
#   make                              # local workstation paths
#   make PLATFORM=hpc                 # cluster paths (see below)
#   make LAMMPS_DIR=... MOL_FORGE=... DFTD3=...   # override any path directly
#   make clean
#   make install PREFIX=/some/prefix  # installs fforge to $(PREFIX)/bin

FC       := mpif90
TARGET   := fforge
PLATFORM ?= local

SRCS := parameters.f90 def_input.f90 def_lammps.f90 def_geom.f90 def_max.f90 def_SNAP.f90 def_force_field.f90 main.f90
OBJS := $(SRCS:.f90=.o)

ifeq ($(PLATFORM),hpc)
LAMMPS_DIR ?= /home/postdoc6/lammps
MOL_FORGE  ?= /home/postdoc6/MolForge
DFTD3      ?= /home/postdoc6/dftd3-lib-0.9

LIB_FIX_DIR   ?= $(HOME)/lib_fix
SCALA_FILE    ?= /opt/cp2k-2024.3/tools/toolchain/install/scalapack-2.2.1/lib/libscalapack.a
OPENBLAS_FILE ?= /lib64/libopenblaso.so.0
MOL_FORGE_LIB := $(shell find $(MOL_FORGE) -name "libMolForge.a" -o -name "libMolForge.so" 2>/dev/null | head -n 1)
LAMMPS_FORTRAN_LIB := $(LAMMPS_DIR)/examples/COUPLE/fortran2/liblammps_fortran.a

FFLAGS  ?= -g -O2
LDFLAGS ?= -Wl,--allow-shlib-undefined
LDLIBS  := -L$(LIB_FIX_DIR) -lstdc++ -latomic -L$(LAMMPS_DIR)/src -llammps_mpi \
           $(MOL_FORGE_LIB) $(LAMMPS_FORTRAN_LIB) -L$(DFTD3)/lib -ldftd3 \
           $(SCALA_FILE) $(OPENBLAS_FILE)
else
LAMMPS_DIR ?= /home/valerio/Desktop/lammps
MOL_FORGE  ?= /home/valerio/Desktop/MolForge
DFTD3      ?= /home/valerio/Downloads/dftd3-lib/dftd3-lib-0.9

FFLAGS  ?= -g
LDFLAGS ?=
LDLIBS  := -L$(DFTD3)/lib -L$(LAMMPS_DIR)/src -L$(LAMMPS_DIR)/examples/COUPLE/fortran2 -L$(MOL_FORGE)/libs \
           -lblas -llammps_fortran -llammps_mpi -lMolForge -ldftd3 -llapack -lmpi_cxx -lstdc++ -lm
endif

INCLUDES := -I$(DFTD3)/lib -I$(MOL_FORGE)/Modules -I$(LAMMPS_DIR)/examples/COUPLE/fortran2

PREFIX ?= /usr/local

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(SRCS)
	$(FC) -c $(FFLAGS) $(INCLUDES) $(SRCS)
	$(FC) $(OBJS) $(FFLAGS) $(LDFLAGS) $(LDLIBS) -o $@
	rm -f $(OBJS)

clean:
	rm -f *.o *.mod $(TARGET)

install: $(TARGET)
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET) $(PREFIX)/bin/

