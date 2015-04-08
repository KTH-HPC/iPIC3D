#include <mpi.h>
#include <fstream>

#include "ParallelIO.h"
#include "MPIdata.h"
#include "TimeTasks.h"
#include "Collective.h"
#include "Grid3DCU.h"
#include "VCtopology3D.h"
#include "Particles3Dcomm.h"
#include "EMfields3D.h"
#include <algorithm>
#include <iostream>
#include <sstream>

/*! Function used to write the EM fields using the parallel HDF5 library */
void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef PHDF5
  timeTasks_set_task(TimeTasks::WRITE_FIELDS);

  stringstream filenmbr;
  string       filename;

  bool         bp;

  /* ------------------- */
  /* Setup the file name */
  /* ------------------- */

  filenmbr << setfill('0') << setw(5) << cycle;
  filename = col->getSaveDirName() + "/" + col->getSimName() + "_" + filenmbr.str() + ".h5";

  /* ---------------------------------------------------------------------------- */
  /* Define the number of cells in the globa and local mesh and set the mesh size */
  /* ---------------------------------------------------------------------------- */

  int nxc = grid->getNXC();
  int nyc = grid->getNYC();
  int nzc = grid->getNZC();

  int    dglob[3] = { col ->getNxc()  , col ->getNyc()  , col ->getNzc()   };
  int    dlocl[3] = { nxc-2,            nyc-2,            nzc-2 };
  double L    [3] = { col ->getLx ()  , col ->getLy ()  , col ->getLz ()   };

  /* --------------------------------------- */
  /* Declare and open the parallel HDF5 file */
  /* --------------------------------------- */

  PHDF5fileClass outputfile(filename, 3, vct->getCoordinates(), vct->getComm());

  bp = false;
  if (col->getParticlesOutputCycle() > 0) bp = true;

  outputfile.CreatePHDF5file(L, dglob, dlocl, bp);

  // write electromagnetic field
  //
  outputfile.WritePHDF5dataset("Fields", "Ex", EMf->getExc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Ey", EMf->getEyc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Ez", EMf->getEzc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Bx", EMf->getBxc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "By", EMf->getByc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Bz", EMf->getBzc(), nxc-2, nyc-2, nzc-2);

  /* ---------------------------------------- */
  /* Write the charge moments for each species */
  /* ---------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream ss;
    ss << is;
    string s_is = ss.str();

    // charge density
    outputfile.WritePHDF5dataset("Fields", "rho_"+s_is, EMf->getRHOcs(is), nxc-2, nyc-2, nzc-2);
    // current
    //outputfile.WritePHDF5dataset("Fields", "Jx_"+s_is, EMf->getJxsc(is), nxc-2, nyc-2, nzc-2);
    //outputfile.WritePHDF5dataset("Fields", "Jy_"+s_is, EMf->getJysc(is), nxc-2, nyc-2, nzc-2);
    //outputfile.WritePHDF5dataset("Fields", "Jz_"+s_is, EMf->getJzsc(is), nxc-2, nyc-2, nzc-2);
  }

  outputfile.ClosePHDF5file();

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}

/*! Function to write the EM fields using the H5hut library. */
void WriteFieldsH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){
  if(col->field_output_is_off())
    return;
#ifdef USEH5HUT
  timeTasks_set_task(TimeTasks::WRITE_FIELDS);

  H5output file;


  /* ---------------- */
  /* Write the fields */
  /* ---------------- */

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetNameCycle(filename, cycle);

  file.OpenFieldsFile("Node", nspec, col->getNxc()+1, col->getNyc()+1, col->getNzc()+1, vct->getCoordinates(), vct->getDims(), vct->getComm());

  file.WriteFields(EMf->getEx(), "Ex", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getEy(), "Ey", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getEz(), "Ez", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBx(), "Bx", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBy(), "By", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBz(), "Bz", grid->getNXN(), grid->getNYN(), grid->getNZN());

  for (int is=0; is<nspec; is++) {
    stringstream  ss;
    ss << is;
    string s_is = ss.str();
    file.WriteFields(EMf->getRHOns(is), "rho_"+ s_is, grid->getNXN(), grid->getNYN(), grid->getNZN());
  }

  file.CloseFieldsFile();

//--- SAVE FIELDS IN THE CELLS:
//
//  file.OpenFieldsFile("Cell", nspec, col->getNxc(), col->getNyc(), col->getNzc(), vct->getCoordinates(), vct->getDims(), vct->getComm());
//
//  file.WriteFields(EMf->getExc(), "Exc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getEyc(), "Eyc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getEzc(), "Ezc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getBxc(), "Bxc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getByc(), "Byc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getBzc(), "Bzc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//
//  for (int is=0; is<nspec; is++) {
//    stringstream  ss;
//    ss << is;
//    string s_is = ss.str();
//    file.WriteFields(EMf->getRHOcs(is), "rhoc_"+ s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
//  }
//
//  file.CloseFieldsFile();
//
//--- END SAVE FIELDS IN THE CELLS.

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}

/*! Function to write the particles using the H5hut library. */
void WritePartclH5hut(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){
#ifdef USEH5HUT
  timeTasks_set_task(TimeTasks::WRITE_PARTICLES);

  H5output file;

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetNameCycle(filename, cycle);

  /* ------------------- */
  /* Write the particles */
  /* ------------------- */

  file.OpenPartclFile(nspec, vct->getComm());
  for (int i=0; i<nspec; i++){
    // this is a hack
    part[i].convertParticlesToSynched();
    file.WriteParticles(i, part[i].getNOP(),
                           part[i].getQall(),
                           part[i].getXall(),
                           part[i].getYall(),
                           part[i].getZall(),
                           part[i].getUall(),
                           part[i].getVall(),
                           part[i].getWall(),
                           vct->getComm());
  }
  file.ClosePartclFile();

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}

#if 0
void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid){
#ifdef USEH5HUT

  H5input infile;
  double L[3] = {col->getLx(), col->getLy(), col->getLz()};

  infile.SetNameCycle(col->getinitfile(), col->getLast_cycle());
  infile.OpenPartclFile(nspec);

  infile.ReadParticles(vct->getCartesian_rank(), vct->getNproc(), vct->getDims(), L, vct->getComm());

  for (int s = 0; s < nspec; s++){
    part[s].allocate(s, infile.GetNp(s), col, vct, grid);

    infile.DumpPartclX(part[s].getXref(), s);
    infile.DumpPartclY(part[s].getYref(), s);
    infile.DumpPartclZ(part[s].getZref(), s);
    infile.DumpPartclU(part[s].getUref(), s);
    infile.DumpPartclV(part[s].getVref(), s);
    infile.DumpPartclW(part[s].getWref(), s);
    infile.DumpPartclQ(part[s].getQref(), s);
  }
  infile.ClosePartclFile();

//--- TEST PARTICLE LECTURE:
//  for (int s = 0; s < nspec; s++){
//    for (int n = 0; n < part[s].getNOP(); n++){
//      double ix = part[s].getX(n);
//      double iy = part[s].getY(n);
//      double iz = part[s].getZ(n);
//      if (ix<=0 || iy<=0 || iz <=0) {
//        cout << " ERROR: This particle has negative position. " << endl;
//        cout << "        n = " << n << "/" << part[s].getNOP();
//        cout << "       ix = " << ix;
//        cout << "       iy = " << iy;
//        cout << "       iz = " << iz;
//      }
//    }
//  }
//--- END TEST

#endif
}
#endif

#if 0
void ReadFieldsH5hut(int nspec, EMfields3D *EMf, Collective *col, VCtopology3D *vct, Grid3DCU *grid){
#ifdef USEH5HUT

  H5input infile;

  infile.SetNameCycle(col->getinitfile(), col->getLast_cycle());

  infile.OpenFieldsFile("Node", nspec, col->getNxc()+1,
                                       col->getNyc()+1,
                                       col->getNzc()+1,
                                       vct->getCoordinates(),
                                       vct->getDims(),
                                       vct->getComm());

  infile.ReadFields(EMf->getEx(), "Ex", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getEy(), "Ey", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getEz(), "Ez", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBx(), "Bx", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBy(), "By", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBz(), "Bz", grid->getNXN(), grid->getNYN(), grid->getNZN());

  for (int is = 0; is < nspec; is++){
    std::stringstream  ss;
    ss << is;
    std::string s_is = ss.str();
    infile.ReadFields(EMf->getRHOns(is), "rho_"+s_is, grid->getNXN(), grid->getNYN(), grid->getNZN());
  }

  infile.CloseFieldsFile();

  // initialize B on centers
  MPI_Barrier(MPI_COMM_WORLD);

  // Comm ghost nodes for B-field
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  grid->interpN2C(EMf->getBxc(), EMf->getBx());
  grid->interpN2C(EMf->getByc(), EMf->getBy());
  grid->interpN2C(EMf->getBzc(), EMf->getBz());

  // Comm ghost cells for B-field
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  // communicate E
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  for (int is = 0; is < nspec; is++)
    grid->interpN2C(EMf->getRHOcs(), is, EMf->getRHOns());

//---READ FROM THE CELLS:
//
//  infile.OpenFieldsFile("Cell", nspec, col->getNxc(),
//                                       col->getNyc(),
//                                       col->getNzc(),
//                                       vct->getCoordinates(),
//                                       vct->getDims(),
//                                       vct->getComm());
//
//  infile.ReadFields(EMf->getExc(), "Exc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getEyc(), "Eyc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getEzc(), "Ezc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getBxc(), "Bxc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getByc(), "Byc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getBzc(), "Bzc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//
//  for (int is = 0; is < nspec; is++){
//    std::stringstream  ss;
//    ss << is;
//    std::string s_is = ss.str();
//    infile.ReadFields(EMf->getRHOcs(is, 0), "rhoc_"+s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
//  }
//
//  infile.CloseFieldsFile();
//
//  // initialize B on nodes
//  grid->interpC2N(EMf->getBx(), EMf->getBxc());
//  grid->interpC2N(EMf->getBy(), EMf->getByc());
//  grid->interpC2N(EMf->getBz(), EMf->getBzc());
//
//  for (int is = 0; is < nspec; is++)
//    grid->interpC2N(EMf->getRHOns(), is, EMf->getRHOcs());
//
//---END READ FROM THE CELLS

#endif
}
#endif

void WriteFieldsVTK(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, const string & outputTag ,int cycle){

	//All VTK output at grid cells excluding ghost cells
	const int nxn  =grid->getNXN(),nyn  = grid->getNYN(),nzn =grid->getNZN();
	const int dimX =col->getNxc() ,dimY = col->getNyc(), dimZ=col->getNzc();
	const double spaceX = dimX>1 ?col->getLx()/(dimX-1) :col->getLx();
	const double spaceY = dimY>1 ?col->getLy()/(dimY-1) :col->getLy();
	const double spaceZ = dimZ>1 ?col->getLz()/(dimZ-1) :col->getLz();
	const int    nPoints = dimX*dimY*dimZ;
	MPI_File     fh;
	MPI_Status   status;

	if (outputTag.find("B", 0) != string::npos || outputTag.find("E", 0) != string::npos
			 || outputTag.find("Je", 0) != string::npos || outputTag.find("Ji", 0) != string::npos){

		const string tags0[]={"B", "E", "Je", "Ji"};
		float writebuffer[nzn-3][nyn-3][nxn-3][3];

		for(int tagid=0;tagid<4;tagid++){
		 if (outputTag.find(tags0[tagid], 0) == string::npos) continue;

		 char   header[1024];
		 if (tags0[tagid].compare("B") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getBxTot(ix+1, iy+1, iz+1);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getByTot(ix+1, iy+1, iz+1);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getBzTot(ix+1, iy+1, iz+1);
					  }

	      //Write VTK header
		  sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Magnetic Field from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d\n"
						   "VECTORS B float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (tags0[tagid].compare("E") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getEx(ix+1, iy+1, iz+1);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getEy(ix+1, iy+1, iz+1);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getEz(ix+1, iy+1, iz+1);
					  }

		  //Write VTK header
		   sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Electric Field from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS E float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (tags0[tagid].compare("Je") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 0);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 0);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 0);
					  }

		  //Write VTK header
		   sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Electron current from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS Je float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 }else if (tags0[tagid].compare("Ji") == 0){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  writebuffer[iz][iy][ix][0] = (float)EMf->getJxs(ix+1, iy+1, iz+1, 1);
						  writebuffer[iz][iy][ix][1] = (float)EMf->getJys(ix+1, iy+1, iz+1, 1);
						  writebuffer[iz][iy][ix][2] = (float)EMf->getJzs(ix+1, iy+1, iz+1, 1);
					  }

		  //Write VTK header
		   sprintf(header, "# vtk DataFile Version 2.0\n"
						   "Ion current from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "VECTORS Ji float\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);
		 }


		 if(EMf->isLittleEndian()){

			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][0],4);
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][1],4);
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][2],4);
					  }
		 }

		  int nelem = strlen(header);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_"<< tags0[tagid] << "_" << cycle << ".vtk";
		  MPI_File_open(vct->getComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

		  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
		  }

	      int err = MPI_File_set_view(fh, disp, MPI_FLOAT, EMf->getProcviewXYZ(), "native", MPI_INFO_NULL);
	      if(err){
	          dprintf("Error in MPI_File_set_view\n");

	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }

	      err = MPI_File_write_all(fh, writebuffer[0][0][0],3*(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
	      if(err){
		      int tcount=0;
		      MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			  dprintf(" wrote %i",  tcount);
	          dprintf("Error in write1\n");
	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }
	      MPI_File_close(&fh);
		}
	}

	if (outputTag.find("rho", 0) != string::npos){

		float writebufferRhoe[nzn-3][nyn-3][nxn-3];
		float writebufferRhoi[nzn-3][nyn-3][nxn-3];
		char   headerRhoe[1024];
		char   headerRhoi[1024];

		for(int iz=0;iz<nzn-3;iz++)
		  for(int iy=0;iy<nyn-3;iy++)
			  for(int ix= 0;ix<nxn-3;ix++){
				  writebufferRhoe[iz][iy][ix] = (float)EMf->getRHOns(ix+1, iy+1, iz+1, 0)*4*3.1415926535897;
				  writebufferRhoi[iz][iy][ix] = (float)EMf->getRHOns(ix+1, iy+1, iz+1, 1)*4*3.1415926535897;
			  }

		//Write VTK header
		sprintf(headerRhoe, "# vtk DataFile Version 2.0\n"
						   "Electron density from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "SCALARS rhoe float\n"
						   "LOOKUP_TABLE default\n",  dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		sprintf(headerRhoi, "# vtk DataFile Version 2.0\n"
						   "Ion density from iPIC3D\n"
						   "BINARY\n"
						   "DATASET STRUCTURED_POINTS\n"
						   "DIMENSIONS %d %d %d\n"
						   "ORIGIN 0 0 0\n"
						   "SPACING %f %f %f\n"
						   "POINT_DATA %d \n"
						   "SCALARS rhoi float\n"
						   "LOOKUP_TABLE default\n",  dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints);

		 if(EMf->isLittleEndian()){
			 for(int iz=0;iz<nzn-3;iz++)
				  for(int iy=0;iy<nyn-3;iy++)
					  for(int ix= 0;ix<nxn-3;ix++){
						  ByteSwap((unsigned char*) &writebufferRhoe[iz][iy][ix],4);
						  ByteSwap((unsigned char*) &writebufferRhoi[iz][iy][ix],4);
					  }
		 }

		  int nelem = strlen(headerRhoe);
		  int charsize=sizeof(char);
		  MPI_Offset disp = nelem*charsize;

		  ostringstream filename;
		  filename << col->getSaveDirName() << "/" << col->getSimName() << "_rhoe_" << cycle << ".vtk";
		  MPI_File_open(vct->getComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

		  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, headerRhoe, nelem, MPI_BYTE, &status);
		  }

	      int err = MPI_File_set_view(fh, disp, MPI_FLOAT, EMf->getProcview(), "native", MPI_INFO_NULL);
	      if(err){
	          dprintf("Error in MPI_File_set_view\n");

	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }

	      err = MPI_File_write_all(fh, writebufferRhoe[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
	      if(err){
		      int tcount=0;
		      MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			  dprintf(" wrote %i",  tcount);
	          dprintf("Error in write1\n");
	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }
	      MPI_File_close(&fh);

	      nelem = strlen(headerRhoi);
	      disp  = nelem*charsize;

	      filename.str("");
	      filename << col->getSaveDirName() << "/" << col->getSimName() << "_rhoi_" << cycle << ".vtk";
	      MPI_File_open(vct->getComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	      MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
		  if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, headerRhoi, nelem, MPI_BYTE, &status);
		  }

	      err = MPI_File_set_view(fh, disp, MPI_FLOAT, EMf->getProcview(), "native", MPI_INFO_NULL);
		  if(err){
			  dprintf("Error in MPI_File_set_view\n");

			  int error_code = status.MPI_ERROR;
			  if (error_code != MPI_SUCCESS) {
				  char error_string[100];
				  int length_of_error_string, error_class;

				  MPI_Error_class(error_code, &error_class);
				  MPI_Error_string(error_class, error_string, &length_of_error_string);
				  dprintf("Error %s\n", error_string);
			  }
		  }

		  err = MPI_File_write_all(fh, writebufferRhoi[0][0],(nxn-3)*(nyn-3)*(nzn-3),MPI_FLOAT, &status);
		  if(err){
			  int tcount=0;
			  MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			  dprintf(" wrote %i",  tcount);
			  dprintf("Error in write1\n");
			  int error_code = status.MPI_ERROR;
			  if (error_code != MPI_SUCCESS) {
				  char error_string[100];
				  int length_of_error_string, error_class;

				  MPI_Error_class(error_code, &error_class);
				  MPI_Error_string(error_class, error_string, &length_of_error_string);
				  dprintf("Error %s\n", error_string);
			  }
		  }
		  MPI_File_close(&fh);
	}


/*


	 if(outputTag.find("rho", 0) != string::npos) {

		 double**** rhoc = EMf->getRHOcs().fetch_arr4();

		for(int is = 0;is<nspec;is++){

			grid->interpN2C(EMf->getRHOcs(), is, EMf->getRHOns());

			//Convert little Endian
			if(EMf->isLittleEndian()){
				 for(int ix= 1;ix<nxc-1;ix++)
					  for(int iy=1;iy<nyc-1;iy++)
						  for(int iz=1;iz<nzc-1;iz++){
							  //rhoc[is][ix][iy][iz]=1.23;
							  //dprintf("rhoc[is][ix][iy][iz]=%f",rhoc[is][ix][iy][iz]);
							  ByteSwap((unsigned char*) &(rhoc[is][ix][iy][iz]),8);}
			}

			//Write VTK header
			char   header[1024];
			sprintf(header, "# vtk DataFile Version 2.0\n"
			"Charge Density from iPIC3D\n"
			"BINARY\n"
			"DATASET STRUCTURED_POINTS\n"
			"DIMENSIONS %d %d %d\n"
			"ORIGIN 0 0 0\n"
			"SPACING %f %f %f\n"
			"POINT_DATA %d \n"
			"SCALARS rho%d double\n"
			"LOOKUP_TABLE default\n", dimX,dimY,dimZ, spaceX,spaceY,spaceZ, nPoints,is);

			int nelem = strlen(header);
			int charsize=sizeof(char);
			MPI_Offset disp = nelem*charsize;

			ostringstream filename;
			filename << col->getSaveDirName() << "/" << col->getSimName() << "_rho"<< is <<"_" << cycle << ".vtk";
			MPI_File_open(vct->getComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

			MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
			if (vct->getCartesian_rank()==0){
			  MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
			}
			former_MPI_Barrier(MPI_COMM_WORLD);


			int err = MPI_File_set_view(fh, disp, MPI_DOUBLE, EMf->getProcview(), "native", MPI_INFO_NULL);
			err = MPI_File_write_all(fh, rhoc[is][0][0],1,EMf->getGhostType(), &status);
			int tcount=0;
			MPI_Get_count(&status, MPI_DOUBLE, &tcount);
			dprintf(" wrote %i",  tcount);
			if(err){
				dprintf("Error in write1\n");
				int error_code = status.MPI_ERROR;
				if (error_code != MPI_SUCCESS) {
					char error_string[100];
					int length_of_error_string, error_class;
					MPI_Error_class(error_code, &error_class);
					MPI_Error_string(error_class, error_string, &length_of_error_string);
					dprintf("Error %s\n", error_string);
				}
			}
			MPI_File_close(&fh);
		}
	 }*/
}

void ByteSwap(unsigned char * b, int n)
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void WritePclsVTK(int nspec, Grid3DCU *grid, Particles3D *part, CollectiveIO *col, VCtopology3D *vct, const string & tag, int cycle){

}
