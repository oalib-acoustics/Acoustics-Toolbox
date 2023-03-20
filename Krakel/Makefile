
%.o: %.f90
	$(FC) -c -I../misc $(FFLAGS) $< 

%.mod: %.f90
	$(FC) -c -I../misc $(FFLAGS) $< 

EXECS=krakel.exe
LIBSMODS=-L../misc -lmisc
MKLPATH = /Library/Frameworks/Intel_MKL.framework/Versions/9.1/lib/universal

all :	$(EXECS)
	@echo "KRAKEL built"
	@echo "*************"
	@echo " "
 
install : $(EXECS)
	for f in $(EXECS) ; do \
		echo "----- Installing $$f.exe"; cp -p $$f ../bin; \
	done

clean:
	-rm -f *.o *.exe

# ****** krakel ******
KRAKELOBJS= krakelmod.o krakel.o dgbdi.o dgbtrf.o dgbtrs.o ../misc/libmisc.a
KRAKELOBJS= krakelmod.o krakel.o dgbdi.o                   ../misc/libmisc.a

krakel.o:  krakelmod.mod krakel.f90

krakel.exe:	$(KRAKELOBJS)
#	$(FC) $(FFLAGS) -o krakel.exe $(KRAKELOBJS) $(LIBSMODS) -L$(MKLPATH) -lmkl_lapack -lmkl -lguide -lpthread
	$(FC) $(FFLAGS) -o krakel.exe $(KRAKELOBJS) $(LIBSMODS) -llapack
