>>>>>>>>./tempexe 4 4 img16bit.tif 550 16 0 l.tif 0
>>>>>>>>diff l.tif truthimg16.tif

RANDMAX=2147483648.000000 ; 4294967296.000000
Command: ./tempexe
Filtering image 'img16bit.tif' using attribute area with lambda=550.000000.
FreeImage version: 3.10.0
Image: Width=2048 Height=1 Depth=1489 Size=3049472 Size2D=2048. nthreads for sorting and quantize the max tree: 4.
Size of int=4. Size of long=8. Size of float=4. Size of double=8.
Min=1122. Max=43524. MULFACTOR=1
Radix Sort (steps=1)
===========================================================
/*** Calculate the quantized image ***/
/*** Build the max tree of the quantized image. ***/
Pilot max-tree built.
/*** Refine ***/ nthreads for the parallel refinement with najman couprie: 4. 
Refined max-tree built.
Init filtering
3) lwb=2285568; upb=3049472.
1) lwb=761856; upb=1523712.
2) lwb=1523712; upb=2285568.
0) lwb=0; upb=761856.
End filtering
Sorting: 0.020000 s.
Create Quantized Image: 0.020000 s.
Quantized Tree: 0.290000 s.
Refinement Tree: 1.230000 s.
Filtering: 0.120000 s.
Wall-Clock time: 1.680000 s.
Image written to 'l.tif'
