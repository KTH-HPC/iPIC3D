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

void WriteFieldsVTK(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){

		  //parallel VTK for B at cells
	      int  size[3], subsize[3], start[3], err;
	      const int nxc  =grid->getNXC(),nyc=grid->getNYC(),nzc=grid->getNZC();
	      const int xLen =vct->getXLEN(), yLen =vct->getYLEN(),zLen =vct->getZLEN();
	      MPI_Datatype  procview, ghosttype;
	      MPI_File      fh;
	      char   header[1024];
	      MPI_Status    status;

	      //convert little Endian to Big Endian
	      arr3_double bx = EMf->getBxTot();//EMf->getBxcWithGhost();
	      arr3_double by = EMf->getByTot();//EMf->getBycWithGhost();
	      arr3_double bz = EMf->getBzTot();//EMf->getBzcWithGhost();
	      const int len = (nxc-2)*(nyc-2)*(nzc-2)*3;
	      const double firstB = EMf->getBxTot(1,1,1);
	      double writebuffer[nzc-2][nyc-2][nxc-2][3];
	      int tmpid=0;
		  for(int iz=0;iz<nzc-2;iz++)
			  for(int iy=0;iy<nyc-2;iy++)
				  for(int ix= 0;ix<nxc-2;ix++){
					  //tmpid = ix*((nyc-2)*(nzc-2)*3) + iy*((nzc-2)*3) + iz*3;
					  writebuffer[iz][iy][ix][0]=EMf->getBxTot(ix+1,iy+1,iz+1);//(bx[ix+1][iy+1][iz+1]);
					  writebuffer[iz][iy][ix][1]=0;//(by[ix+1][iy+1][iz+1]);
					  writebuffer[iz][iy][ix][2]=0;//(bz[ix+1][iy+1][iz+1]);//dprintf("writebuffer[ix][iy][iz][2] %i %i %i tmpid%i %f",ix,iy,iz,tmpid,writebuffer[tmpid]);
				  }

		  dprintf("Finish initializing");

		  if(EMf->isLittleEndian()){
	//		  for(int id= 0;id<len;id++){
	//			   ByteSwap((unsigned char*) &writebuffer[id],8);
	//		   }
			  for(int iz=0;iz<nzc-2;iz++)
				  for(int iy=0;iy<nyc-2;iy++)
					  for(int ix= 0;ix<nxc-2;ix++){
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][0],8);
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][1],8);
						  ByteSwap((unsigned char*) &writebuffer[iz][iy][ix][2],8);
					  }
			   dprintf("islittleEndian, finish converting");
		   }

	      //create local subarray to exclude ghost cells
	      size[0] = nxc;size[1] = nyc;size[2] = nzc;
	      subsize[0] = nxc-2;subsize[1] = nyc-2;subsize[2] = nzc-2;
	      start[0]=1;start[1]=1;start[2]=1;
	      dprintf("ghost2dsubarr size %i %i %i \n", size[0], size[1], size[2]);
	      dprintf("ghost2dsubarr subsize %i %i %i \n", subsize[0], subsize[1], subsize[2]);
	      dprintf("ghost2dsubarr start %i %i %i \n", start[0], start[1], start[2]);

	      MPI_Datatype  testview;//=MPI_DOUBLE;
	      MPI_Type_contiguous(nxc*nyc*nzc*3,MPI_DOUBLE, &testview);
	      MPI_Type_commit(&testview);

	      MPI_Datatype  etype;//=MPI_DOUBLE;
	      MPI_Type_contiguous(3,MPI_DOUBLE, &etype);
	      MPI_Type_commit(&etype);
	      err=MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C, etype, &ghosttype);
	      if(err){
	          dprintf("Error in ghost2dsubarr\n");

	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }
	      MPI_Type_commit(&ghosttype);

	      //create process file view
	      size[0] = (nxc-2)*xLen;size[1] = (nyc-2)*yLen;size[2] = (nzc-2)*zLen;
	      start[0]= vct->getCoordinates(0)*subsize[0];
	      start[1]= vct->getCoordinates(1)*subsize[1];
	      start[2]= vct->getCoordinates(2)*subsize[2];
	      dprintf("proc2dsubarr size %i %i %i \n",    size[0], size[1], size[2]);
	      dprintf("proc2dsubarr subsize %i %i %i \n", subsize[0], subsize[1], subsize[2]);
	      dprintf("proc2dsubarr start %i %i %i \n",   start[0], start[1], start[2]);

	      err=MPI_Type_create_subarray(3, size, subsize, start,MPI_ORDER_C, etype, &procview);
	      if(err){
	          dprintf("Error in proc2dsubarr\n");

	          int error_code = status.MPI_ERROR;
	          if (error_code != MPI_SUCCESS) {
	              char error_string[100];
	              int length_of_error_string, error_class;

	              MPI_Error_class(error_code, &error_class);
	              MPI_Error_string(error_class, error_string, &length_of_error_string);
	              dprintf("Error %s\n", error_string);
	          }
	      }
	      MPI_Type_commit(&procview);

	      //Write VTK header B
	       sprintf(header, "# vtk DataFile Version 2.0\n"
	       "Magnetic Field from iPIC3D\n"
	       "BINARY\n"
	       "DATASET STRUCTURED_POINTS\n"
	       "DIMENSIONS %d %d %d\n"
	       "ORIGIN 0 0 0\n"
	       "SPACING %f %f %f\n"
	       "POINT_DATA %d \n"
	       "VECTORS B double\n", size[0],size[1],size[2],grid->getDX(),grid->getDY(),grid->getDZ(),size[0]*size[1]*size[2]);

	      int nelem = strlen(header);
	      int charsize=sizeof(char);
	      MPI_Offset disp = nelem*charsize;

	      /*
	      dprintf("disp %i nelem %i charsize %i",disp ,nelem, charsize);

	      MPI_Type_size(MPI_DOUBLE, &charsize);
	      dprintf("MPI_DOUBLE size %i", charsize);

	      MPI_Type_size(MPI_FLOAT, &charsize);
	      dprintf("MPI_FLOAT size %i", charsize);

	      MPI_Type_size(MPI_INT, &charsize);
	      dprintf("MPI_INT size %i", charsize);
	       */

	      MPI_Datatype  headertype;
	      MPI_Type_contiguous(disp, MPI_BYTE, &headertype);
	      MPI_Type_commit(&headertype);

	      ostringstream filename;
	      filename << col->getSaveDirName() << "/" << col->getSimName() << "_B_" << cycle << ".vtk";
	      MPI_File_open(vct->getComm(),filename.str().c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	      MPI_File_set_view(fh, 0, MPI_BYTE, headertype, "native", MPI_INFO_NULL);
	      if (vct->getCartesian_rank()==0){
	    	  MPI_File_write(fh, header, nelem, MPI_BYTE, &status);
	      }
	      former_MPI_Barrier(MPI_COMM_WORLD);


	      //err = MPI_File_set_view(fh, disp, etype, procview, "native", MPI_INFO_NULL);
	      //err = MPI_File_set_view(fh, disp, MPI_DOUBLE, testview, "native", MPI_INFO_NULL);
	      err = MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
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
	      //err = MPI_File_write_all(fh, writebuffer,len,MPI_DOUBLE, &status);//test->fetch_arr3()
	      err = MPI_File_write_all(fh, writebuffer[0][0][0],len,MPI_DOUBLE, &status);//test->fetch_arr3()
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

	      MPI_Type_free(&procview);
	      MPI_Type_free(&ghosttype);
	      MPI_Type_free(&headertype);
	      MPI_Type_free(&etype);

	      MPI_File_close(&fh);
	      dprintf("Exit parallel VTK");
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

void WritePartclVTK(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){

}
