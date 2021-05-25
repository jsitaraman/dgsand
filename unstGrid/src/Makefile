# Generated automatically from Makefile.in by configure.
ALL: default
##### User configurable options #####

SHELL       = /bin/sh
ARCH        = LINUX
CC          = gcc
F90         = gfortran
CLINKER     = $(CC)
CP          = cp

MAKE        = make --no-print-directory
#
srcdir= $(PWD)

### End User configurable options ###

.SUFFIXES:	.c .C .f90
CFLAGS	  = $(OPTFLAGS) 
CCFLAGS	  = $(CFLAGS)
LIBS	  = -lm -lgfortran

EXECS = ugrid make_domain

ALL: $(EXECS)

default: $(EXECS)
OBJS =  ugrid_test.o ugrid_sub.o
ugrid:  $(OBJS)
	$(F90) $(CFLAGS) $(OBJS) $(LIBS) -o ugrid
	cp ugrid ../bin

make_domain: make_domain.o
	$(CC) $(CFLAGS) make_domain.o $(LIBS) -o make_domain
	cp make_domain ../bin

%.o:%.c
	$(CC) -c $(CFLAGS) $< -o $*.o

%.o:%.f90
	$(F90) -c $(FFLAGS) $< -o $*.o
	

clean:
	rm *.o
	rm $(EXECS)
