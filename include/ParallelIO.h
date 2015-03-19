
#ifndef __PARALLELIO_H__
#define __PARALLELIO_H__

#ifdef USEH5HUT
#  include "H5hut-io.h"
#endif

#ifdef PHDF5
#  include "phdf5.h"
#endif

#include "ipicfwd.h"

void WriteFieldsH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle);
void WritePartclH5hut(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);

void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid);
void ReadFieldsH5hut(int nspec, EMfields3D *EMf,       Collective *col, VCtopology3D *vct, Grid3DCU *grid);

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);


void WriteFieldsVTK(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle);
void WritePartclVTK(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle);
void ByteSwap(unsigned char * b, int n);
#endif
