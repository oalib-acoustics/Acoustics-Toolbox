
EXECUTABLES = bounce.exe kraken.exe krakenc.exe

BOUNCE_OBJECTS = KrakencMod.o  BCImpedancecMod.o                                 Scattering.o bounce.o ../misc/libmisc.a
BOUNCE_LIBS = -L../misc -lmisc

KRAKEN_OBJECTS = KrakenMod.o   BCImpedanceMod.o  InverseIterationMod.o RootFinderBrent.o Scattering.o kraken.o ../misc/libmisc.a
KRAKEN_LIBS = -L../misc -lmisc

KRAKENC_OBJECTS = KrakencMod.o BCImpedancecMod.o InverseIterationMod.o           Scattering.o krakenc.o ../misc/libmisc.a
KRAKENC_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "KRAKEN built"
	@echo "************"
	@echo " "

install: all
	for f in $(EXECUTABLES) ; do \
		echo "----- Installing $$f"; cp -p $$f ../bin; \
	done

clean:
	-rm -f *.o *.mod *.exe *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 

# Declare explicit dependencies (and explicit build rules as needed)

bounce.exe: $(BOUNCE_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BOUNCE_OBJECTS) $(BOUNCE_LIBS)

kraken.exe: $(KRAKEN_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(KRAKEN_OBJECTS) $(KRAKEN_LIBS)

krakenc.exe: $(KRAKENC_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(KRAKENC_OBJECTS) $(KRAKENC_LIBS)

InverseIterationMod.o inverseiterationmod.mod: InverseIterationMod.f90
	$(FC) -c $(FFLAGS) $<

KrakenMod.o krakenmod.mod: KrakenMod.f90
	$(FC) -c $(FFLAGS) $< 

KrakencMod.o krakencmod.mod: KrakencMod.f90
	$(FC) -c $(FFLAGS) $< 

BCImpedanceMod.o BCImpedanceMod.mod: BCImpedanceMod.f90
	$(FC) -c $(FFLAGS) $< 

BCImpedancecMod.o BCImpedancecMod.mod: BCImpedancecMod.f90
	$(FC) -c $(FFLAGS) $< 

bounce.o: krakencmod.mod bounce.f90

kraken.o:  krakenmod.mod  kraken.f90

krakenc.o: krakencmod.mod krakenc.f90

Scattering.o: Scattering.f90

RootFinderBrent.o rootfinderbrent.mod: RootFinderBrent.f90


