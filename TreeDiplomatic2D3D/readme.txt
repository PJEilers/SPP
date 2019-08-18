README for usage of parallel max-tree library,

Usage:
make or make static
When making a shared library, remember to change the LD_LIBRARY_PATH to ~/lib at runtime: 
LD_LIBRARY_PATH=~/lib
export LD_LIBRARY_PATH
./compile.sh (All library files should be in $HOME/lib)
./build.sh (All library files should be in $HOME/lib)
./main <nthreads> <input image> <lambda> <bits per pixel> <output image> [attrib] 
Where attrib is: 
        0 - Area
        1 - Area of min. enclosing rectangle
        2 - Square of diagonal of min. enclosing rectangle
        3 - Moment of Inertia
        4 - (Moment of Inertia) / (area)^2
        5 - Mean X position
        6 - Mean Y position
        7 - Mean Z position

This main program uses cfitsio to read images: https://heasarc.gsfc.nasa.gov/fitsio/ and treediplomatic library for parallel tree-building and filtering.

