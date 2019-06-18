/*************************************************************************************/
/* setest.c  :   Performs erosions, dilations, openings and closings on 3D .fld files*/
/* Usage     :   Used with volume structures - see: Volume.c                         */
/* Author    :   Georgios K. Ouzounis                                                */
/* Dep.      :   IWI - University of Groningen                                       */
/* Vs&Date   :   vs 1.0 - 10th Feb. 2005                                             */  
/*************************************************************************************/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Convert.h"
#include "Volume.h"
#include "se.h"
#include "types.h"
#include "avs_io.h"

int main(int argc, char *argv[])
{
   ushort g[10000], h[10000], h2[10000];
   VolumeStruct VolStruct_In, VolStruct_Out;
   avs_header header;
   int choice, datatype;
      
   char *in_filename, *out_filename="out.fld";
   ulong k = 3;

   if (argc<2)
   {
      printf("Usage: %s <input volume> [k] [output volume]\n", argv[0]);
      return(0);
   }
   in_filename = argv[1];
   if (argc>=3)  k = atol(argv[2]);
   if (argc>=4)  out_filename = argv[3];

   /* Read Input Volume */
   readAVS(in_filename, &header, &VolStruct_In);
   datatype = VolStruct_In.DataType;
   if(datatype==1)
      ConvertByte2Short(&VolStruct_In);
   
   /* Preapare output Volume */
   InitVolume(&VolStruct_Out, VolStruct_In.VolWidth, VolStruct_In.VolHeight, 
               VolStruct_In.VolDepth, VolStruct_In.DataType);

   printf("\nChoose operation : 1. Erosion\n");
   printf("                   2. Dilation\n");
   printf("                   3. Opening\n");
   printf("                   4. Closing\n");
   scanf("%d",&choice);
   
   switch(choice){
      case 1:
         VolGrayErodeHor(&VolStruct_In, k, g, h, h2, &VolStruct_Out);
         VolGrayErodeVer(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayErodeDepth(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         break;
      case 2:  
         VolGrayDilateHor(&VolStruct_In, k, g, h, h2, &VolStruct_Out);
         VolGrayDilateVer(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayDilateDepth(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
	 break;
      case 3:
         VolGrayErodeHor(&VolStruct_In, k, g, h, h2, &VolStruct_Out);
         VolGrayErodeVer(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayErodeDepth(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
	 VolGrayDilateHor(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayDilateVer(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayDilateDepth(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
	 break;
      case 4:
         VolGrayDilateHor(&VolStruct_In, k, g, h, h2, &VolStruct_Out);
         VolGrayDilateVer(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayDilateDepth(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
	 VolGrayErodeHor(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayErodeVer(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
         VolGrayErodeDepth(&VolStruct_Out, k, g, h, h2, &VolStruct_Out);
	 break;
      default:
         printf("Invalid selection of operator\n");
   }
   if(datatype==1)
      ConvertShort2Byte(&VolStruct_Out);
      
   writeAVS(out_filename, &header, &VolStruct_Out); 
   
   DestroyVolume(&VolStruct_In);
   DestroyVolume(&VolStruct_Out);
   return(0);
} /* main */
