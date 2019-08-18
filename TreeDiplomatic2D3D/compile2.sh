#!/bin/bash
#gcc -O2 -std=gnu99 -pedantic -march=opteron -fexpensive-optimizations -funroll-loops -c -L$HOME/lib -I$HOME/lib main.c quanttree.c handleimages.c quantizedimage.c refinetree.c radixsort.c filter.c -lpthread -lcfitsio -lm -lfreeimage

gcc -O3 -Wall -std=gnu99 -pedantic -fexpensive-optimizations -fPIC -funroll-loops -c -L$HOME/lib -I$HOME/lib main.c quanttree.c datafunctions.c handleimages.c quantizedimage.c refinetree.c radixsort.c filter.c -lpthread -lcfitsio -lfreeimage

