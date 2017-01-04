
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

all: libts.a $(EXECS)

install: all
	for f in $(EXECS) ; do \
	    ( if test -e $$f.exe ; then \
		echo "Installing $$f.exe"; cp -p $$f.exe ../bin; \
	    else \
		echo "$$f.exe not found"; \
	    fi ; ) \
	done

clean : clean_execs
	-rm -f *.o
	-rm -f libts.a

clean_execs : 
	for f in $(EXECS) ; do \
	  if test -e $$f.exe ; then \
		(echo "Removing $$f" ; rm $$f.exe) ; \
	  fi ;  \
	done

# *** General purpose subroutines ***

OBJECTS=source.o \
	cans.o \
	preenv.o \
	hilbert.o \
	bndpasc.o \
	cfft.o \
	rfft.o \
	polyr.o

libts.a: libts.a($(OBJECTS))
	echo "$(OBJECTS)"; \
	ar rv libts.a $?; \
	ranlib libts.a
	@echo "Routines needed for time series built"

iso: iso.f90 library
	$(FC) $(FFLAGS) -o iso.exe iso.f90 -L../misc -lmisc


