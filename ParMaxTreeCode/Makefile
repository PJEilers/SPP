# Set the following macros to their appropriate values.

CC       = gcc
# Uncomment for use on opteron-based machines
%CFLAGS   = -ansi -pedantic -Wall -O3 -march=opteron -fexpensive-optimizations -funroll-loops
CFLAGS   = -ansi -pedantic -Wall -O3 -march=i686 -fexpensive-optimizations -funroll-loops
#CFLAGS   = -g -O0
CPPFLAGS =
LDFLAGS  = -O
SETESTOBJS = ReadWriteData.o Volume.o avs_io.o

# --------------------------------------------
# Below this line nothing ought to be changed 

VERSION = 0.0.0
INCLUDES= 
LIBS    = -lpthread -lm
DEFS    = -DVERSION="\"$(VERSION)\""
COMPILE = $(CC) -c $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
LINK    = $(CC) $(LDFLAGS) -o $@
SHELL   = /bin/sh

all: mt3D1

clean:
	rm mt3D2 *~ *.o

mt3D2: mt3D2.o $(SETESTOBJS)
	${LINK} mt3D2.o $(SETESTOBJS) $(LIBS)

mt3D1: mt3D1.o $(SETESTOBJS)
	${LINK} mt3D1.o $(SETESTOBJS) $(LIBS)

raw2avs: raw2avs.o $(SETESTOBJS)
	${LINK} raw2avs.o $(SETESTOBJS) $(LIBS)

12bitto8bit: 12bitto8bit.o $(SETESTOBJS)
	${LINK} 12bitto8bit.o $(SETESTOBJS) $(LIBS)

mt3D1area: mt3D1area.o $(SETESTOBJS)
	${LINK} mt3D1area.o $(SETESTOBJS) $(LIBS)

mt3D2.o: mt3D2.c decision3D.c
	${COMPILE} mt3D2.c 

mt3D1.o: mt3D1.c decision3D.c
	${COMPILE} mt3D1.c 

raw2avs.o: raw2avs.c
	${COMPILE} raw2avs.c 

12bitto8bit.o: 12bitto8bit.c
	${COMPILE} 12bitto8bit.c 

mt3D1area.o: mt3D1area.c decision3D.c
	${COMPILE} mt3D1area.c 

Volume.o: Volume.c 
	${COMPILE} Volume.c 

ReadWriteData.o: ReadWriteData.c 
	${COMPILE} ReadWriteData.c 

avs_io.o: avs_io.c 
	${COMPILE} avs_io.c 
