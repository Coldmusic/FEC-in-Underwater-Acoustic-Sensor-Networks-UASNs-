
%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

%.mod: %.f90
	$(FC) $(FFLAGS) -c $< 

EXECS=bellhop.exe
LIBSMODS=-L../misc -lmisc

all: $(EXECS)
	@echo "Bellhop built"
	@echo "*************"
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

# ***  bellhop ***

BELLHOPOBJS =	bellMod.o \
		angleMod.o \
                bdryMod.o  \
		sspMod.o \
		ArrMod.o \
		beampatternMod.o \
		readin.o  \
		splined.o \
		munk.o \
		../misc/libmisc.a

bellhop.exe: bellhop.f90 $(BELLHOPOBJS)
	$(FC) $(FFLAGS) -o bellhop.exe bellhop.f90 $(BELLHOPOBJS) $(LIBSMODS)

readin.o:	readin.f90
	$(FC) $(FFLAGS) -c readin.f90

