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
CFLAGS	  = -g $(OPTFLAGS) 
CCFLAGS	  = $(CFLAGS)
LIBS	  = -lm 

EXECS = dgsand

ALL: $(EXECS)

default: $(EXECS)
OBJS =  dgsand.o
dgsand:  $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o dgsand

%.o:%.c *.h
	$(CC) -c $(CFLAGS) $< -o $*.o

%.o:%.f90
	$(F90) -c $(FFLAGS) $< -o $*.o
	

clean:
	rm *.o
	rm $(EXECS)
