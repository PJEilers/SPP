#!/bin/bash
#./compile.sh
g++ -o main -L$HOME/lib -I$HOME/lib  main.o handleimages.o -ltreediplomatic -lpthread -lcfitsio -lfreeimage
rm *.o
