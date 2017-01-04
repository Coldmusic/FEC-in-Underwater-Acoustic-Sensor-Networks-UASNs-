
%.o: %.f90
	$(FC) -c $(FFLAGS) $< 

all:	munk.o readin.o RefCoMod.o SdRdRMod.o twersky.o \
	WriteHeader.o ReadHeader.o subtab.o subtabd.o sorti.o sortid.o \
	mergev.o weight.o splinec.o polyz.o xerror.o \
	errout.o libmisc.a

clean:
	-rm -f libmisc.a
	-rm -f *.o
	-rm -f *.mod
	-rm -f *.lib
	-rm -f *.lst
	-rm -f *.BAK
	-rm -f *.fwd

install: all

# *** General purpose subroutines ***

SdRdRMod.mod:	SdRdRMod.f90
	$(FC) $(FFLAGS) SdRdRMod.f90

RefCoMod.mod:	RefCoMod.f90
	$(FC) $(FFLAGS) RefCoMod.f90


libmisc.a:	munk.o \
	readin.o \
	RefCoMod.o \
	SdRdRMod.o  \
	twersky.o \
	WriteHeader.o \
	ReadHeader.o \
	subtab.o subtabd.o \
	sorti.o sortid.o \
	mergev.o \
	weight.o \
	splinec.o \
        polyz.o \
	xerror.o \
	errout.o
	ar -r libmisc.a errout.o WriteHeader.o mergev.o munk.o polyz.o \
		ReadHeader.o readin.o RefCoMod.o sorti.o sortid.o splinec.o \
		SdRdRMod.o subtab.o subtabd.o twersky.o weight.o
	ranlib libmisc.a

	@echo "Routines needed by SCOOTER\SPARC\BELLHOP\BOUNCE built"
