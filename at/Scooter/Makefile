
%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

%.mod: %.f90
	$(FC) $(FFLAGS) -c $< 

EXECS=scooter.exe sparc.exe fields.exe 

LIBSMODS= -L../tslib  -L../misc  -lts -lmisc

all:	$(EXECS)
	@echo "Scooter and Sparc built"
	@echo "***********************"
	@echo " "

clean : clean_execs
	-rm -f *.o *.mod

clean_execs : 
	for f in $(EXECS) ; do \
	  if test -e $$f ; then \
		(echo "Removing $$f" ; rm $$f) ; \
	  fi ;  \
	done

install: all
	for f in $(EXECS) ; do \
	    ( if test -e $$f ; then \
		echo "Installing $$f"; cp -p $$f ../bin; \
	    else \
		echo "$$f not found"; \
	    fi ; ) \
	done

SCOOTEROBJS = scomod.o facnpZ.o ../misc/libmisc.a
scooter.exe: scooter.f90 $(SCOOTEROBJS)
	$(FC) $(FFLAGS) -o scooter.exe scooter.f90 $(SCOOTEROBJS) $(LIBSMODS)

SPARCOBJS=spamod.o facnp.o ../tslib/libts.a ../misc/libmisc.a
sparc.exe:	sparc.f90 $(SPARCOBJS)
	$(FC) $(FFLAGS) -o sparc.exe sparc.f90 $(SPARCOBJS) $(LIBSMODS)

FIELDSOBJS=hts.o pade.o polyc.o ../misc/libmisc.a
fields.exe:	fields.f90 $(FIELDSOBJS)
	$(FC) $(FFLAGS) -o fields.exe fields.f90 $(FIELDSOBJS) $(LIBSMODS)


