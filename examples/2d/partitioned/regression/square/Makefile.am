
bin_PROGRAMS = square

square_SOURCES = square.c
square_CPPFLAGS = -I../../../../../include
square_LDADD = ../../../../../src/.libs/librjmcmc.a

if WITHOPENMPI
square_LDADD += -lmpi
endif

EXTRA_DIST = mkdata.py \
	data.txt


