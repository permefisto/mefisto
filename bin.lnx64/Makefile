#
# Makefile pour obtenir xvtest1 TEST de Mefisto/xvue
# Alain Perronnet       Fevrier 1996
#
LIBS= /usr/lib/Motif1.2/libXm.sl \
	/usr/lib/X11R5/libXhp11.sl \
	/usr/lib/X11R5/libXt.sl \
	/usr/lib/X11R5/libXext.a \
	/usr/lib/X11R5/libX11.sl \
	/usr/lib/X11R5/libX11.a

F77= f77

FFLAGS= +ppu -g -I../

all: xvtest1
	rm *.o ; echo "" ; echo "xvtest1 pret a etre execute"

xvtest1:xvuelc.o xvouvrir.o xvinit.o palcde.o coutip.o \
	codtei.o tiprvb.o chargefonte.o choixfonte.o
	$(F77) $(FFLAGS) ../prpr/xvtest1.f *.o $(LIBS) -o $@

.f.o:	xvouvrir.f xvinit.f palcde.f coutip.f \
	codtei.f tiprvb.f chargefonte.f choixfonte.f
	$(F77) $(FFLAGS) -c $<

.c.o: xvuelc.c
	ccxvue
