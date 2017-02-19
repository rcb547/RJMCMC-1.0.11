
FCC = gfortran
FCFLAGS = -c -g -Wall

AR = ar
ARFLAGS = rcs

OBJS = dbdsdc.o \
	dbdsqr.o \
	dgebd2.o \
	dgebrd.o \
	dgelq2.o \
	dgelqf.o \
	dgeqr2.o \
	dgeqrf.o \
	dgesdd.o \
	disnan.o \
	dlabrd.o \
	dlacpy.o \
	dlaed6.o \
	dlaisnan.o \
	dlamch.o \
	dlamrg.o \
	dlange.o \
	dlanst.o \
	dlapy2.o \
	dlarfb.o \
	dlarf.o \
	dlarfg.o \
	dlarft.o \
	dlartg.o \
	dlas2.o \
	dlascl.o \
	dlasd0.o \
	dlasd1.o \
	dlasd2.o \
	dlasd3.o \
	dlasd4.o \
	dlasd5.o \
	dlasd6.o \
	dlasd7.o \
	dlasd8.o \
	dlasda.o \
	dlasdq.o \
	dlasdt.o \
	dlaset.o \
	dlasq1.o \
	dlasq2.o \
	dlasq3.o \
	dlasq4.o \
	dlasq5.o \
	dlasq6.o \
	dlasr.o \
	dlasrt.o \
	dlassq.o \
	dlasv2.o \
	dorg2r.o \
	dorgbr.o \
	dorgl2.o \
	dorglq.o \
	dorgqr.o \
	dorm2r.o \
	dormbr.o \
	dorml2.o \
	dormlq.o \
	dormqr.o \
	ieeeck.o \
	iladlc.o \
	iladlr.o \
	ilaenv.o \
	iparmq.o \
	lsame.o \
	xerbla.o

TARGETS = liblapack.a

all : $(TARGETS)

liblapack.a : $(OBJS)	
	$(AR) $(ARFLAGS) liblapack.a $(OBJS)

%.o : %.f
	$(FCC) $(FCFLAGS) -o $*.o $*.f

clean :
	rm -f $(OBJS) $(TARGETS)



