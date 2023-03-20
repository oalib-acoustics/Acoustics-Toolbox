
EXECUTABLES = bellhop.exe bellhop3d.exe

BELLHOP_OBJECTS   = bellhopMod.o angleMod.o ArrMod.o bdryMod.o   sspMod.o ReflectMod.o                       WriteRay.o              influence.o               Step.o   ReadEnvironmentBell.o bellhop.o ../misc/libmisc.a
BELLHOP_LIBS = -L../misc -lmisc

BELLHOP3D_OBJECTS = bellhopMod.o angleMod.o ArrMod.o bdry3DMod.o sspMod.o ReflectMod.o Reflect3DMod.o Cone.o WriteRay.o RayNormals.o influence.o influence3D.o Step3DMod.o ReadEnvironmentBell.o bellhop3D.o ../misc/libmisc.a
BELLHOP3D_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "Bellhop built"
	@echo "*************"
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

bellhop.exe: $(BELLHOP_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BELLHOP_OBJECTS) $(BELLHOP_LIBS)

bellhop3d.exe: $(BELLHOP3D_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(BELLHOP3D_OBJECTS) $(BELLHOP3D_LIBS)

angleMod.o anglemod.mod: angleMod.f90
	$(FC) -c $(FFLAGS) $< 

ArrMod.o arrmod.mod: ArrMod.f90
	$(FC) -c $(FFLAGS) $< 

bdryMod.o bdrymod.mod: bdryMod.f90
	$(FC) -c $(FFLAGS) $< 

bdry3DMod.o bdry3dmod.mod: bdry3DMod.f90
	$(FC) -c $(FFLAGS) $< 

bellhopMod.o bellhopmod.mod: bellhopMod.f90
	$(FC) -c $(FFLAGS) $< 

Cone.o cone.mod: Cone.f90
	$(FC) -c $(FFLAGS) $<

influence.o influence.mod: influence.f90
	$(FC) -c $(FFLAGS) $< 

influence3D.o influence3d.mod: influence3D.f90
	$(FC) -c $(FFLAGS) $< 

RayNormals.o raynormals.mod: RayNormals.f90
	$(FC) -c $(FFLAGS) $<

ReadEnvironmentBell.o readenvironmentbell.mod: ReadEnvironmentBell.f90
	$(FC) -c $(FFLAGS) $< 

ReflectMod.o   reflectmod.mod:   ReflectMod.f90   Cone.o
	$(FC) -c $(FFLAGS) $<

Reflect3DMod.o reflect3dmod.mod: Reflect3DMod.f90 Cone.o
	$(FC) -c $(FFLAGS) $<

sspMod.o sspmod.mod: sspMod.f90
	$(FC) -c $(FFLAGS) $< 

Step.o step.mod: Step.f90
	$(FC) -c $(FFLAGS) $<

Step3DMod.o step3dmod.mod: Step3DMod.f90
	$(FC) -c $(FFLAGS) $<

WriteRay.o writeray.mod: WriteRay.f90
	$(FC) -c $(FFLAGS) $<

bellhop.o:   arrmod.mod anglemod.mod bdrymod.mod   readenvironmentbell.mod bellhopmod.mod sspmod.mod reflectmod.mod                  writeray.mod raynormals.mod influence.mod   step.mod      bellhop.f90

bellhop3D.o: arrmod.mod anglemod.mod bdry3dmod.mod readenvironmentbell.mod bellhopmod.mod sspmod.mod reflectmod.mod reflect3dmod.mod writeray.mod raynormals.mod influence3d.mod step3dmod.mod bellhop3D.f90

