
EXECUTABLES = field.exe field3d.exe

FIELD_OBJECTS = ReadModes.o EvaluateMod.o EvaluateADMod.o EvaluateCMMod.o ../misc/libmisc.a field.o
FIELD_LIBS = -L../misc -lmisc

FIELD3D_OBJECTS = ReadModes.o ElementMod.o Evaluate3DMod.o EvaluateGBMod.o EvaluatepdqMod.o ../misc/libmisc.a field3d.o 
FIELD3D_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "KRAKEN field routines built"
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

EvaluateMod.o evaluatemod: EvaluateMod.f90

EvaluateCMMod.o evaluatecmmod.mod: readmodes.mod EvaluateCMMod.f90

EvaluateADMod.o evaluateadmod.mod: readmodes.mod EvaluateADMod.f90

EvaluatepdqMod.o evaluatepdqmod.mod: elementmod.mod EvaluatepdqMod.f90

Evaluate3DMod.o evaluate3dmod.mod: Evaluate3DMod.f90

EvaluateGBMod.o evaluategbmod.mod: EvaluateGBMod.f90

ElementMod.o elementmod.mod: ElementMod.f90

ReadModes.o readmodes.mod: ReadModes.f90
	$(FC) -c $(FFLAGS) $<    # this rule is necessary, otherwise the make utility thinks *.mod is a modula file and uses m2c

field.o:   ReadModes.mod evaluatemod.mod evaluateadmod.mod evaluatecmmod.mod field.f90

field3d.o: ReadModes.mod elementmod.mod evaluate3dmod.mod evaluatepdqmod.mod evaluategbmod.mod field3d.f90

field.exe: $(FIELD_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELD_OBJECTS) $(FIELD_LIBS)

field3d.exe: $(FIELD3D_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(FIELD3D_OBJECTS) $(FIELD3D_LIBS)
