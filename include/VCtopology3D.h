/***************************************************************************
  VCtopology3D.h  -  a 3D Virtual cartesian topology
  A virtual topology is a mechanism for naming the processes
  in a communicator in a way that fits the communication
  pattern better. Since our processes will communicate mainly
  with the nearest neighbours after the fashion of a two-dimensional
  grid, we create a virtual topology to reflect this fact
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL Luveun
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef VCtopology3D_H
#define VCtopology3D_H

#include "mpi.h"

/**
 *  
 * Virtual cartesian topology
 * A virtual topology is a mechanism for naming the processes
 * in a communicator in a way that fits the communication
 * pattern better. Since our processes will communicate mainly
 * with the nearest neighbours after the fashion of a two-dimensional
 * grid, we create a virtual topology to reflect this fact
 * @version 2.0
 */

class Collective;

class VCtopology3D //:public VirtualTopology3D
{
public:
  /** constructor: Define topology parameters: dimension, domain decomposition,... */
  VCtopology3D(const Collective& col);
  /** destructor */
  ~VCtopology3D();
  /** Find the neighbors in the new communicator  */
  void setup_vctopology(MPI_Comm comm_old);
  /** Print topology info */
  void Print();
  /** Print the mapping of topology */
  void PrintMapping();

  int getXLEN()const{ return (XLEN); }
  int getYLEN()const{ return (YLEN); }
  int getZLEN()const{ return (ZLEN); }
  int getNprocs()const{ return (nprocs); }
  bool getPERIODICX()const{ return (PERIODICX); }
  bool getPERIODICY()const{ return (PERIODICY); }
  bool getPERIODICZ()const{ return (PERIODICZ); }

  // legacy names
  //
  int getCartesian_rank()const{ return (cartesian_rank); }
  int getXleft_neighbor()const{ return (xleft_neighbor); }
  int getXright_neighbor()const{ return (xright_neighbor); }
  int getYleft_neighbor()const{ return (yleft_neighbor); }
  int getYright_neighbor()const{ return (yright_neighbor); }
  int getZleft_neighbor()const{ return (zleft_neighbor); }
  int getZright_neighbor()const{ return (zright_neighbor); }
  int getXleft_neighbor_P()const{ return (xleft_neighbor); }
  int getXright_neighbor_P()const{ return (xright_neighbor); }
  int getYleft_neighbor_P()const{ return (yleft_neighbor); }
  int getYright_neighbor_P()const{ return (yright_neighbor); }
  int getZleft_neighbor_P()const{ return (zleft_neighbor); }
  int getZright_neighbor_P()const{ return (zright_neighbor); }

  // new interface
  //
  int getXleft()const{ return (xleft_neighbor); }
  int getXrght()const{ return (xright_neighbor); }
  int getYleft()const{ return (yleft_neighbor); }
  int getYrght()const{ return (yright_neighbor); }
  int getZleft()const{ return (zleft_neighbor); }
  int getZrght()const{ return (zright_neighbor); }

  bool isPeriodicXlower()const{ return _isPeriodicXlower; }
  bool isPeriodicXupper()const{ return _isPeriodicXupper; }
  bool isPeriodicYlower()const{ return _isPeriodicYlower; }
  bool isPeriodicYupper()const{ return _isPeriodicYupper; }
  bool isPeriodicZlower()const{ return _isPeriodicZlower; }
  bool isPeriodicZupper()const{ return _isPeriodicZupper; }

  bool noXleftNeighbor()const{ return _noXleftNeighbor; }
  bool noXrghtNeighbor()const{ return _noXrghtNeighbor; }
  bool noYleftNeighbor()const{ return _noYleftNeighbor; }
  bool noYrghtNeighbor()const{ return _noYrghtNeighbor; }
  bool noZleftNeighbor()const{ return _noZleftNeighbor; }
  bool noZrghtNeighbor()const{ return _noZrghtNeighbor; }

  bool hasXleftNeighbor()const{ return !_noXleftNeighbor; }
  bool hasXrghtNeighbor()const{ return !_noXrghtNeighbor; }
  bool hasYleftNeighbor()const{ return !_noYleftNeighbor; }
  bool hasYrghtNeighbor()const{ return !_noYrghtNeighbor; }
  bool hasZleftNeighbor()const{ return !_noZleftNeighbor; }
  bool hasZrghtNeighbor()const{ return !_noZrghtNeighbor; }

  bool isBoundaryProcess()const{ return _isBoundaryProcess; }

  bool isXlower()const{ return coordinates[0]==0; }
  bool isYlower()const{ return coordinates[1]==0; }
  bool isZlower()const{ return coordinates[2]==0; }
  bool isXupper()const{ return coordinates[0]==dims[0]-1; }
  bool isYupper()const{ return coordinates[1]==dims[1]-1; }
  bool isZupper()const{ return coordinates[2]==dims[2]-1; }

  bool getcVERBOSE()const{ return (cVERBOSE); }
  int getCoordinates(int dir)const{ return (coordinates[dir]); }
  const int *getCoordinates()const{ return (coordinates); }
  const int *getDims()const{ return dims; }
  //const int *getDivisions()const{ return getDims(); } // old name
  int getPeriods(int dir)const{ return (periods[dir]); }
  MPI_Comm getComm()const{ return (CART_COMM); }


private:
  /** New communicator with virtual cartesian topology */
  MPI_Comm CART_COMM;
  /** New communicator with virtual cartesian topology for Particles*/
  //MPI_Comm CART_COMM_P;
  /** MPI status during sending and receiving communication */
  MPI_Status status;
  /** Direction X for shift MPI_Cart_Shift*/
  int XDIR;
  /** Direction Y for shift MPI_Cart_Shift*/
  int YDIR;
  /** Direction Z for shift MPI_Cart_Shift*/
  int ZDIR;
  /** RIGHT = +    upwards   shift */
  int RIGHT;
  /** LEFT  = -    downwards shift */
  int LEFT;
  /** dimension of virtual topology */
  int PROCDIM;
  /** number of subdomains - Direction X */
  int XLEN;
  /** number of subdomains - Direction Y */
  int YLEN;
  /** number of subdomains - Direction Z */
  int ZLEN;
  /** nprocs = number of processors */
  int nprocs;
  /** periodicity on boundaries - DIRECTION X*/
  bool PERIODICX;
  /** periodicity on boundaries - DIRECTION Y*/
  bool PERIODICY;
  /** periodicity on boundaries - DIRECTION Z*/
  bool PERIODICZ;
  /** periodicity on boundaries - DIRECTION X*/
  //bool PERIODICX_P;
  /** periodicity on boundaries - DIRECTION Y*/
  //bool PERIODICY_P;
  /** periodicity on boundaries - DIRECTION Z*/
  //bool PERIODICZ_P;
  /** rank may be reordered     */
  int reorder;
  /** arrays for Create_Cart_create  */
  int dims[3]; // i.e. divisions
  /** periodicity */
  int periods[3];
  //int periods_P[3];
  /** coordinates on processors grid */
  int coordinates[3];
  /** cartesian rank */
  int cartesian_rank;
  /** cartesian rank of XLEFT neighbor */
  int xleft_neighbor;
  /** cartesian rank of XRIGHT neighbor */
  int xright_neighbor;
  /** cartesian rank of YLEFT neighbor */
  int yleft_neighbor;
  /** cartesian rank of YRIGHT neighbor */
  int yright_neighbor;
  /** cartesian rank of ZRIGHT neighbor */
  int zleft_neighbor;
  /** cartesian rank of ZLEFT neighbor */
  int zright_neighbor;
  /** cartesian rank of XLEFT neighbor */
  //int xleft_neighbor_P;
  /** cartesian rank of XRIGHT neighbor */
  //int xright_neighbor_P;
  /** cartesian rank of YLEFT neighbor */
  //int yleft_neighbor_P;
  /** cartesian rank of YRIGHT neighbor */
  //int yright_neighbor_P;
  /** cartesian rank of ZRIGHT neighbor */
  //int zleft_neighbor_P;
  /** cartesian rank of ZLEFT neighbor */
  //int zright_neighbor_P;
  
  /** indicators of whether this is a periodic boundary */
  bool _isPeriodicXlower;
  bool _isPeriodicXupper;
  bool _isPeriodicYlower;
  bool _isPeriodicYupper;
  bool _isPeriodicZlower;
  bool _isPeriodicZupper;

  /** indicators of whether this lacks a neighbor */
  bool _noXrghtNeighbor;
  bool _noXleftNeighbor;
  bool _noYrghtNeighbor;
  bool _noYleftNeighbor;
  bool _noZrghtNeighbor;
  bool _noZleftNeighbor;

  int _isBoundaryProcess;

  /** if cVERBOSE == true, print to the screen all the comunication */
  bool cVERBOSE;
};

typedef VCtopology3D VirtualTopology3D;

#endif
