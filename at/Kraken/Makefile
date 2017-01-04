
%.o: %.f90
	$(FC) $(FFLAGS) -c $< 


%.mod: %.f90
	$(FC) $(FFLAGS) -c $< 

EXECS=kraken.exe krakenc.exe bounce.exe field.exe field3d.exe modasc.exe modbin.exe
MODULES=krakmod.f90 krakcmod.f90 ElementMod.f90
LIBSMODS=-L../misc -lmisc

all :	$(EXECS)
	@echo "KRAKEN built"
	@echo "*************"
	@echo " "
 
install : $(EXECS)
	for f in $(EXECS) ; do \
	    ( if test -e $$f ; then \
		echo "Installing $$f"; cp -p $$f ../bin; \
	    else \
		echo "$$f not found"; \
	    fi ; ) \
	done

clean : clean_execs
	-rm -f *.o *.mod

clean_execs : 
	for f in $(EXECS) ; do \
	  if test -e $$f ; then \
		(echo "Removing $$f" ; rm $$f) ; \
	  fi ;  \
	done

modules: $(MODULES)

# ****** kraken ******
KRAKENOBJS=zsecx.o zbrentx.o bcimpedance.o kuping.o PekerisRoot.o sinvitd.o krakmod.o ../misc/libmisc.a
kraken.exe:	kraken.f90 $(KRAKENOBJS) krakmod.o
	$(FC) $(FFLAGS) -o kraken.exe kraken.f90 $(KRAKENOBJS) $(LIBSMODS)

bcimpedance.o:	bcimpedance.f90 krakmod.o
	$(FC) $(FFLAGS) -c bcimpedance.f90

# ****** krakenc ******
KRAKENCOBJS=zseccx.o sinvitz.o bcimpedancec.o	kuping.o PekerisRoot.o krakcmod.o ../misc/libmisc.a
krakenc.exe: krakenc.f90 $(KRAKENCOBJS) krakcmod.o
	$(FC) $(FFLAGS) -o krakenc.exe krakenc.f90 $(KRAKENCOBJS) $(LIBSMODS)
bcimpedancec.o: bcimpedancec.f90 krakcmod.o ../misc/libmisc.a
	$(FC) $(FFLAGS) -c bcimpedancec.f90

# ****** bounce ******
BOUNCEOBJS=bcimpedancec.o kuping.o PekerisRoot.o krakcmod.o ../misc/libmisc.a
bounce.exe:	bounce.f90 $(BOUNCEOBJS)
	$(FC) $(FFLAGS) -o bounce.exe bounce.f90 $(BOUNCEOBJS) $(LIBSMODS)

# ****** field ******
FIELDOBJS=getmodes.o PekerisRoot.o eval.o evalad.o evalcm.o ../misc/libmisc.a
field.exe:	field.f90 $(FIELDOBJS)
	$(FC) $(FFLAGS) -o field.exe field.f90 $(FIELDOBJS) $(LIBSMODS)

# ****** field3d ******
FIELD3DOBJS=ElementMod.o evalpdq.o eval3d.o evalgb.o getmodes.o PekerisRoot.o ../misc/libmisc.a
field3d.exe: field3d.f90 $(FIELD3DOBJS) ElementMod.o
	$(FC) $(FFLAGS) -o field3d.exe field3d.f90 $(FIELD3DOBJS) $(LIBSMODS)

# ****** modasc/bin

modasc.exe: modasc.f90
	$(FC) $(FFLAGS) -o modasc.exe modasc.f90
modbin.exe: modbin.f90
	$(FC) $(FFLAGS) -o modbin.exe modbin.f90

# ****** covar ******
COVAROBJS=getmod.o PekerisRoot.o ../misc/libmisc.a
covar.exe:	covar.f90 $(COVAROBJS)
	$(FC) $(FFLAGS) -o covar.exe covar.f90 $(COVAROBJS) $(LIBSMODS)

