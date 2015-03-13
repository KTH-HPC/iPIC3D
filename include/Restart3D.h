// Developed by Stefano Markidis, and Giovanni Lapenta
#ifndef _RESTART3D_H_
#define _RESTART3D_H_

#include "ipicfwd.h"
#include "arraysfwd.h"
#include "aligned_vector.h"
//#include "Particles3Dcomm.h"
//#include "PSKhdf5adaptor.h"

#include <iosfwd>

using std::string;

/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part);

/** this restart function writes the last restart with the last cycle */
void writeRESTART(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part, bool fool);

/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART_ES(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part);

/** this restart function writes the last restart with the last cycle */
void writeRESTART_ES(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part, bool fool);

void read_field_restart(
  const Collective* col,
  const VCtopology3D* vct,
  const Grid* grid,
  arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
  arr3_double Ex, arr3_double Ey, arr3_double Ez,
  array4_double* rhons_, int ns);

void read_particles_restart(
    const Collective* col,
    const VCtopology3D* vct,
    int species_number,
    vector_double& u,
    vector_double& v,
    vector_double& w,
    vector_double& q,
    vector_double& x,
    vector_double& y,
    vector_double& z,
    vector_double& t);

#endif // _RESTART_H_
