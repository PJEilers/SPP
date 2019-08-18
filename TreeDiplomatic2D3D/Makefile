default: shared

CC = gcc
TARGET = treediplomatic
CFLAGS = -O3 -Wall -std=gnu99 -pedantic -fexpensive-optimizations -funroll-loops -c
EXTRAFLAGS = -fPIC
STATICLIB = lib$(TARGET).a
SHAREDLIB = lib$(TARGET).so
LIBRARIES = -lpthread -lcfitsio -lfreeimage

shared:
	$(CC) $(CFLAGS) $(EXTRAFLAGS) -L$(HOME)/lib -I$(HOME)/lib Source/*.c $(LIBRARIES)
	$(CC) -shared -o $(SHAREDLIB) *.o
	cp $(SHAREDLIB) $(HOME)/lib
	cp Source/TreeDiplomatic.h $(HOME)/lib
	rm $(SHAREDLIB)
	rm *.o
static:
	$(CC) $(CFLAGS) -L$(HOME)/lib -I$(HOME)/lib main.c Source/*.c $(LIBRARIES)
	ar rcs $(STATICLIB) *.o
	cp $(STATICLIB) $(HOME)/lib
	cp Source/TreeDiplomatic.h $(HOME)/lib
	rm $(STATICLIB)
	rm *.o