#ifndef AMR_MACROS_H
#define AMR_MACROS_H

#define ARRAY_2D(P,LO_X,HI_X,LO_Y,HI_Y,I,J) P((J-LO_Y)*(HI_X-LO_X+1)+I-LO_X+1)
#define ARRAY_3D(P,LO_X,HI_X,LO_Y,HI_Y,LO_Z,HI_Z,I,J,K) P((K-LO_Z)*(HI_X-LO_X+1)*(HI_Y-LO_Y+1) + (J-LO_Y)*(HI_X-LO_X+1)+I-LO_X+1)

#ifdef CUDA 
#endif

#endif
