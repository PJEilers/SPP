#ifndef SE_H
#define SE_H

void DilateHorLine(ushort *f, ulong width, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r);
void ErodeHorLine( ushort *f, ulong width, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r);
void DilateVerLine(ushort *f, ulong width, ulong height, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r);
void ErodeVerLine( ushort *f, ulong width, ulong height, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r);
void DilateDepthLine(ushort *f, ulong width, ulong height, ulong depth, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r);
void ErodeDepthLine( ushort *f, ulong width, ulong height, ulong depth, ulong k, ushort *g, ushort *h, ushort *h2, ushort *r);
void VolGrayDilateHor(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out);
void VolGrayErodeHor( VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out);
void VolGrayDilateVer(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out);
void VolGrayErodeVer( VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out);
void VolGrayDilateDepth(VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out);
void VolGrayErodeDepth( VolumeStruct *vol_in, ulong k, ushort *g, ushort *h, ushort *h2, VolumeStruct *vol_out);




#endif /* SE_H */
