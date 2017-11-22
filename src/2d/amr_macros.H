#ifndef AMR_MACROS_H
#define AMR_MACROS_H

#ifdef CUDA 

#define FFLUXPTR_D  1
#define CFLUXPTR_D  2
#define NODE_DATA_SIZE  2

#define NESTLEVEL_D  1
#define NDILO_D  2
#define NDIHI_D  3
#define NDJLO_D  4
#define NDJHI_D  5
#define NODE_STAT_SIZE  5

#endif

#endif
