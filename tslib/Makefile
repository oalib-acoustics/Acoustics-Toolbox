
AR = ar

# Intel compiler needed the following definition for the archiver:
# AR = xiar

ARFLAGS = r

RANLIB = ranlib

OBJECTS = bandpassc.o cans.o cfft.o hilbert.o preenv.o rfft.o sourceMod.o

all: $(OBJECTS) libts.a

install: all

clean:
	-rm -f *.a *.mod *.o *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 

# Declare explicit dependencies (and explicit build rules as needed)

libts.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@
	@echo "Routines needed for time series library built"

bandpassc.o: bandpassc.f90
cans.o: cans.f90
cfft.o: cfft.f90
hilbert.o: hilbert.f90
iso.o: iso.f90
preenv.o: preenv.f90
rfft.o: rfft.f90
#source.o: source.f90

sourceMod.o sourcemod.mod: sourceMod.f90
	$(FC) -c $(FFLAGS) $< 
