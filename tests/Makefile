
ENVFIL = pekeris twersky  double scholte flused elsed atten norm -fal ice kupert kuperb kuperm -f

install: clean $(ENVFIL)

$(ENVFIL):
	-rm -f field.flp
	cp fieldbat.flp field.flp
	(kraken.exe  < $@; mv $@.prt pcprt\$@r.prt; rm -f $@.mod; mv $@.asc $@r.asc)
	(krakenc.exe < $@; mv $@.prt pcprt\$@c.prt; rm -f $@.mod; mv $@.asc $@c.asc)
	(scooter.exe < $@; mv $@.prt pcprt\$@s.prt; rm -f $@.grn; mv $@.asc $@s.asc)
	#fc pcprt\$@r.prt $@r.prt|more
	#fc pcprt\$@c.prt $@c.prt|more
	#fc pcprt\$@s.prt $@s.prt|more
	#c:\examdiff $@r.prt sgiprt\$@r.prt
	#c:\examdiff $@c.prt sgiprt\$@c.prt
	#c:\examdiff $@s.prt sgiprt\$@s.prt
	-rm -f field.flp

clean:
	-rm -f *.prt
	-rm -f *.mod
	-rm -f *.grn
	-rm -f *.shd
	-rm -f *.rts
	-rm -f *.ray
	-rm -f *.brc
	-rm -f *.irc
	-rm -f *.dat
	-rm -f *.asc
	-rm -f fort.*
	-rm -f *.moA
	-rm -f *.mat
	(cd halfspace; 	make clean)
	(cd arctic;	make clean)
	(cd Munk;	make clean)
	(cd MunkLeaky;	make clean)
	(cd calib;	make clean)
	(cd Gulf;	make clean)
	(cd Dickins;	make clean)
	(cd sduct;	make clean)
	(cd SBCX;       make clean)
	(cd BeamPattern; make clean)
	(cd TabRefCoef;	make clean)
	(cd head;	make clean)
	(cd free;	make clean)
	(cd PointLine;	make clean)
	(cd ParaBot;	make clean)
	(cd Ellipse;	make clean)
	(cd 3DAtlantic;	make clean)
	(cd TLslices;	make clean)
	(cd wedge;	make clean)
	(cd BeamPattern;	make clean)
	(cd Noise;	make clean)


