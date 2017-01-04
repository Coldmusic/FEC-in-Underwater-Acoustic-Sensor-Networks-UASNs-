
%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

%.mod: %.f90
	$(FC) $(FFLAGS) -c $< 

EXECS=krakel.exe
MODULES=krakelmod.f90
LIBSMODS=-L../misc -lmisc
MKLPATH = /Library/Frameworks/Intel_MKL.framework/Versions/9.1/lib/universal

all :	$(EXECS)
	@echo "KRAKEL built"
	@echo "*************"
	@echo " "
 
install : $(EXECS)
	for f in $(EXECS) ; do \
	    ( if test -e $$f ; then \
		echo "Installing $$f.exe"; cp -p $$f ../bin; \
	    else \
		echo "$$f not found"; \
	    fi ; ) \
	done

clean : clean_execs
	-rm -f *.o

clean_execs : 
	for f in $(EXECS) ; do \
	  if test -e $$f ; then \
		(echo "Removing $$f" ; rm $$f) ; \
	  fi ;  \
	done

modules: $(MODULES)

# ****** krakel ******
KRAKELOBJS= krakelmod.o dgbdi.o dgbtrf.o dgbtrs.o ../Kraken/zsecx.o ../misc/libmisc.a
KRAKELOBJS= krakelmod.o dgbdi.o ../Kraken/zsecx.o ../misc/libmisc.a

krakel.exe:	krakel.f90 $(KRAKELOBJS)
#	$(FC) $(FFLAGS) -o krakel.exe krakel.f90 $(KRAKELOBJS) $(LIBSMODS) -L$(MKLPATH) -lmkl_lapack -lmkl -lguide -lpthread
	$(FC) $(FFLAGS) -o krakel.exe krakel.f90 $(KRAKELOBJS) $(LIBSMODS) -llapack

