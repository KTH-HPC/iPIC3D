
#include "mpi.h"
#include "MPIdata.h"
#include "iPic3D.h"
#include "TimeTasks.h"
#include "ipicdefs.h"
#include "debug.h"
#include "Parameters.h"
#include "ompdefs.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"
#include "Timing.h"
//
#ifndef NO_HDF5
#include "ParallelIO.h"
#include "WriteOutputParallel.h"
#include "OutputWrapperFPP.h"
#endif


#include <algorithm> //required for std::swap
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

#include <iostream>
#include <fstream>
#include <sstream>

#include "Moments.h" // for debugging

using namespace iPic3D;
//MPIdata* iPic3D::c_Solver::mpi=0;

c_Solver::~c_Solver()
{
  delete col; // configuration parameters ("collectiveIO")
  delete vct; // process topology
  delete grid; // grid
  delete EMf; // field
  delete outputWrapperFPP;

  // delete particles
  //
  if(part)
  {
    for (int i = 0; i < ns; i++)
    {
      // placement delete
      part[i].~Particles3D();
    }
    free(part);
  }

  delete [] Ke;
  delete [] momentum;
  delete [] Qremoved;
  delete my_clock;
}

int c_Solver::Init(int argc, char **argv) {
  #if defined(__MIC__)
  assert_eq(DVECWIDTH,8);
  #endif
  // get MPI data
  //
  // c_Solver is not a singleton, so the following line was pulled out.
  //MPIdata::init(&argc, &argv);
  //
  // initialized MPI environment
  // nprocs = number of processors
  // myrank = rank of tha process*/
  Parameters::init_parameters();
  //mpi = &MPIdata::instance();
  nprocs = MPIdata::get_nprocs();
  myrank = MPIdata::get_rank();

  col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  verbose = col->getVerbose();
  restart_cycle = col->getRestartOutputCycle();
  SaveDirName = col->getSaveDirName();
  RestartDirName = col->getRestartDirName();
  restart_status = col->getRestart_status();
  ns = col->getNs();            // get the number of particle species involved in simulation
  first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart
  // initialize the virtual cartesian topology 
  vct = new VCtopology3D(*col);
  // Check if we can map the processes into a matrix ordering defined in Collective.cpp
  if (nprocs != vct->getNprocs()) {
    if (myrank == 0) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << "x" << vct->getZLEN() << " matrix: Change XLEN,YLEN, ZLEN in method VCtopology3D.init()" << endl;
      MPIdata::instance().finalize_mpi();
      return (1);
    }
  }
  // We create a new communicator with a 3D virtual Cartesian topology
  vct->setup_vctopology(MPI_COMM_WORLD);
  {
    stringstream num_proc_ss;
    num_proc_ss << vct->getCartesian_rank();
    num_proc_str = num_proc_ss.str();
  }
  // initialize the central cell index

#ifdef BATSRUS
  // set index offset for each processor
  col->setGlobalStartIndex(vct);
#endif

  // Print the initial settings to stdout and a file
  if (myrank == 0) {
    MPIdata::instance().Print();
    vct->Print();
    col->Print();
    col->save();
  }
  // Create the local grid
  grid = new Grid3DCU(col, vct);  // Create the local grid
  EMf = new EMfields3D(col, grid, vct);  // Create Electromagnetic Fields Object

  if      (col->getCase()=="GEMnoPert") EMf->initGEMnoPert();
  else if (col->getCase()=="ForceFree") EMf->initForceFree();
  else if (col->getCase()=="GEM")       EMf->initGEM();
#ifdef BATSRUS
  else if (col->getCase()=="BATSRUS")   EMf->initBATSRUS();
#endif
  else if (col->getCase()=="Dipole")    EMf->initDipole();
  else if (col->getCase()=="RandomCase") {
    EMf->initRandomField();
    if (myrank==0) {
      cout << "Case is " << col->getCase() <<"\n";
      cout <<"total # of particle per cell is " << col->getNpcel(0) << "\n";
    }
  }
  else {
    if (myrank==0) {
      cout << " =========================================================== " << endl;
      cout << " WARNING: The case '" << col->getCase() << "' was not recognized. " << endl;
      cout << "          Runing simulation with the default initialization. " << endl;
      cout << " =========================================================== " << endl;
    }
    EMf->init();
  }

  // OpenBC
  EMf->updateInfoFields();

  // Allocation of particles
  // part = new Particles3D[ns];
  part = (Particles3D*) malloc(sizeof(Particles3D)*ns);
  for (int i = 0; i < ns; i++)
  {
    new(&part[i]) Particles3D(i,col,vct,grid);
    //part[i] = new Particles3D(i, col, vct, grid);
    //part[i].allocate(i, col, vct, grid);
  }

  // Initial Condition for PARTICLES if you are not starting from RESTART
  if (restart_status == 0) {
    // wave = new Planewave(col, EMf, grid, vct);
    // wave->Wave_Rotated(part); // Single Plane Wave
    for (int i = 0; i < ns; i++)
    {
      if      (col->getCase()=="ForceFree") part[i].force_free(EMf);
#ifdef BATSRUS
      else if (col->getCase()=="BATSRUS")   part[i].MaxwellianFromFluid(EMf,col,i);
#endif
      else                                  part[i].maxwellian(EMf);
      part[i].reserve_remaining_particle_IDs();
    }
  }

  if (col->getWriteMethod() == "shdf5")
  {
    #ifndef NO_HDF5
  outputWrapperFPP = new OutputWrapperFPP;
  fetch_outputWrapperFPP().init_output_files(col,vct,grid,EMf,part,ns);
    #endif
  }

  Ke = new double[ns];
  momentum = new double[ns];
  cq = SaveDirName + "/ConservedQuantities.txt";
  if (myrank == 0) {
    ofstream my_file(cq.c_str());
    my_file.close();
  }
  

  Qremoved = new double[ns];

  my_clock = new Timing(myrank);

  return 0;
}

void c_Solver::sortParticles() {
  
  for(int species_idx=0; species_idx<ns; species_idx++)
    part[species_idx].sort_particles_serial();
  
}

void c_Solver::CalculateMoments() {

  timeTasks_set_main_task(TimeTasks::MOMENTS);

  pad_particle_capacities();

  // vectorized assumes that particles are sorted by mesh cell
  if(Parameters::get_VECTORIZE_MOMENTS())
  {
    switch(Parameters::get_MOMENTS_TYPE())
    {
      case Parameters::SoA:
        // since particles are sorted,
        // we can vectorize interpolation of particles to grid
        convertParticlesToSoA();
        sortParticles();
        EMf->sumMoments_vectorized(part);
        break;
      case Parameters::AoS:
        convertParticlesToAoS();
        sortParticles();
        EMf->sumMoments_vectorized_AoS(part);
        break;
      default:
        unsupported_value_error(Parameters::get_MOMENTS_TYPE());
    }
  }
  else
  {
    if(Parameters::get_SORTING_PARTICLES())
      sortParticles();
    switch(Parameters::get_MOMENTS_TYPE())
    {
      case Parameters::SoA:
        EMf->setZeroPrimaryMoments();
        convertParticlesToSoA();
        EMf->sumMoments(part);
        break;
      case Parameters::AoS:
        EMf->setZeroPrimaryMoments();
        convertParticlesToAoS();
        EMf->sumMoments_AoS(part);
        break;
      case Parameters::AoSintr:
        EMf->setZeroPrimaryMoments();
        convertParticlesToAoS();
        EMf->sumMoments_AoS_intr(part);
        break;
      default:
        unsupported_value_error(Parameters::get_MOMENTS_TYPE());
    }
  }
  //for (int i = 0; i < ns; i++)
  //{
  //  EMf->sumMomentsOld(part[i], grid, vct);
  //}
  EMf->setZeroDerivedMoments();
  // sum all over the species
  EMf->sumOverSpecies();
  // Fill with constant charge the planet
  if (col->getCase()=="Dipole") {
    EMf->ConstantChargePlanet(col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
  }
  // Set a constant charge in the OpenBC boundaries
  //EMf->ConstantChargeOpenBC();
  // calculate densities on centers from nodes
  EMf->interpDensitiesN2C();
  // calculate the hat quantities for the implicit method
  EMf->calculateHatFunctions();
}

//! MAXWELL SOLVER for Efield
void c_Solver::CalculateField() {
  timeTasks_set_main_task(TimeTasks::FIELDS);
  EMf->updateInfoFields();
  // calculate the E field
  EMf->calculateE();
}

//! MAXWELL SOLVER for Bfield (assuming Efield has already been calculated)
void c_Solver::CalculateB() {
  timeTasks_set_main_task(TimeTasks::FIELDS);
  // calculate the B field
  EMf->calculateB();
}

/*  -------------- */
/*!  Particle mover */
/*  -------------- */
bool c_Solver::ParticlesMover()
{
  // move all species of particles
  {
    timeTasks_set_main_task(TimeTasks::PARTICLES);
    // Should change this to add background field
    EMf->set_fieldForPcls();

    pad_particle_capacities();

    #pragma omp parallel
    {
    for (int i = 0; i < ns; i++)  // move each species
    {
      // #pragma omp task inout(part[i]) in(grid) target_device(booster)
      // should merely pass EMf->get_fieldForPcls() rather than EMf.
      // use the Predictor Corrector scheme to move particles
      switch(Parameters::get_MOVER_TYPE())
      {
        case Parameters::SoA:
          part[i].mover_PC(EMf);
          break;
        //case Parameters::SoA_vec_resort:
        //  part[i].mover_PC_vectorized(EMf);
        //  break;
        case Parameters::AoS:
          part[i].mover_PC_AoS(EMf);
          break;
        case Parameters::AoS_Relativistic:
        	part[i].mover_PC_AoS_Relativistic(EMf);
        	break;
        case Parameters::AoSintr:
          part[i].mover_PC_AoS_vec_intr(EMf);
          break;
        case Parameters::AoSvec:
          part[i].mover_PC_AoS_vec(EMf);
          break;
        //case Parameters::AoS_vec_onesort:
        //  part[i].mover_PC_AoS_vec_onesort(EMf);
        //  break;
        default:
          unsupported_value_error(Parameters::get_MOVER_TYPE());
      }
      // overlap initial communication of electrons with moving of ions
      #pragma omp master
      {
    	  //Should integrate BC into separate_and_send_particles
    	  part[i].openbc_particles();
    	  part[i].separate_and_send_particles();
      }
    }
    }//End of omp parallel
    for (int i = 0; i < ns; i++)  // communicate each species
    {
      //part[i].communicate_particles();
      part[i].recommunicate_particles_until_done(1);
    }
  }

  /* -------------------------------------- */
  /* Repopulate the buffer zone at the edge */
  /* -------------------------------------- */

  for (int i=0; i < ns; i++) {
    if (col->getRHOinject(i)>0.0)
      part[i].repopulate_particles();
  }

  /* --------------------------------------- */
  /* Remove particles from depopulation area */
  /* --------------------------------------- */

  if (col->getCase()=="Dipole") {
    for (int i=0; i < ns; i++)
      Qremoved[i] = part[i].deleteParticlesInsideSphere(col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
  }
  return (false);
}

void c_Solver::WriteRestart(int cycle)
{
  bool do_WriteRestart = (cycle % restart_cycle == 0 && restart_cycle!=0);
  if(!do_WriteRestart)
    return;

  #ifndef NO_HDF5
  convertParticlesToSynched(); // hack
  // write the RESTART file
  if (restart_cycle!=0)
      writeRESTART(RestartDirName, myrank, cycle, ns, vct, col, grid, EMf, part, 0); // without ,0 add to restart file
  #endif
}

// write the conserved quantities
void c_Solver::WriteConserved(int cycle) {
  if(cycle % col->getDiagnosticsOutputCycle() == 0)
  {
    Eenergy = EMf->getEenergy();
    Benergy = EMf->getBenergy();
    TOTenergy = 0.0;
    TOTmomentum = 0.0;
    for (int is = 0; is < ns; is++) {
      Ke[is] = part[is].getKe();
      TOTenergy += Ke[is];
      momentum[is] = part[is].getP();
      TOTmomentum += momentum[is];
    }
    if (myrank == 0) {
      ofstream my_file(cq.c_str(), fstream::app);
      if(cycle == 0) my_file << "\t" << "\t" << "\t" << "Total_Energy" << "\t" << "Momentum" << "\t" << "Eenergy" << "\t" << "Benergy" << "\t" << "Kenergy" << endl;

      my_file << cycle << "\t" << "\t" << (Eenergy + Benergy + TOTenergy) << "\t" << TOTmomentum << "\t" << Eenergy << "\t" << Benergy << "\t" << TOTenergy << endl;
      my_file.close();
    }
  }
}

void c_Solver::WriteVelocityDistribution(int cycle)
{
  // Velocity distribution
  //if(cycle % col->getVelocityDistributionOutputCycle() == 0)
  {
    for (int is = 0; is < ns; is++) {
      double maxVel = part[is].getMaxVelocity();
      long long *VelocityDist = part[is].getVelocityDistribution(nDistributionBins, maxVel);
      if (myrank == 0) {
        ofstream my_file(ds.c_str(), fstream::app);
        my_file << cycle << "\t" << is << "\t" << maxVel;
        for (int i = 0; i < nDistributionBins; i++)
          my_file << "\t" << VelocityDist[i];
        my_file << endl;
        my_file.close();
      }
      delete [] VelocityDist;
    }
  }
}

// This seems to record values at a grid of sample points
//
void c_Solver::WriteVirtualSatelliteTraces()
{
  if(ns <= 2) return;
  assert_eq(ns,4);

  ofstream my_file(cqsat.c_str(), fstream::app);
  const int nx0 = grid->get_nxc_r();
  const int ny0 = grid->get_nyc_r();
  const int nz0 = grid->get_nzc_r();
  for (int isat = 0; isat < nsat; isat++) {
    for (int jsat = 0; jsat < nsat; jsat++) {
      for (int ksat = 0; ksat < nsat; ksat++) {
        int index1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
        int index2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
        int index3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
        my_file << EMf->getBx(index1, index2, index3) << "\t" << EMf->getBy(index1, index2, index3) << "\t" << EMf->getBz(index1, index2, index3) << "\t";
        my_file << EMf->getEx(index1, index2, index3) << "\t" << EMf->getEy(index1, index2, index3) << "\t" << EMf->getEz(index1, index2, index3) << "\t";
        my_file << EMf->getJxs(index1, index2, index3, 0) + EMf->getJxs(index1, index2, index3, 2) << "\t" << EMf->getJys(index1, index2, index3, 0) + EMf->getJys(index1, index2, index3, 2) << "\t" << EMf->getJzs(index1, index2, index3, 0) + EMf->getJzs(index1, index2, index3, 2) << "\t";
        my_file << EMf->getJxs(index1, index2, index3, 1) + EMf->getJxs(index1, index2, index3, 3) << "\t" << EMf->getJys(index1, index2, index3, 1) + EMf->getJys(index1, index2, index3, 3) << "\t" << EMf->getJzs(index1, index2, index3, 1) + EMf->getJzs(index1, index2, index3, 3) << "\t";
        my_file << EMf->getRHOns(index1, index2, index3, 0) + EMf->getRHOns(index1, index2, index3, 2) << "\t";
        my_file << EMf->getRHOns(index1, index2, index3, 1) + EMf->getRHOns(index1, index2, index3, 3) << "\t";
      }}}
  my_file << endl;
  my_file.close();
}

void c_Solver::WriteFields(int cycle) {
  if(col->field_output_is_off())
    return;
  #ifdef NO_HDF5
    eprintf("must compile with HDF5");
  #else
  if(cycle % (col->getFieldOutputCycle()) == 0 || cycle == first_cycle)
  {
    timeTasks_set_task(TimeTasks::WRITE_FIELDS);
    if (col->getWriteMethod() == "Parallel") {
        WriteOutputParallel(grid, EMf, col, vct, cycle);
    }
    else // OUTPUT to large file, called proc**
    {
        // Pressure tensor is available
        fetch_outputWrapperFPP().append_output(
          "Eall + Ball + rhos + Jsall + pressure", cycle);
    }
  }
  #endif
}

void c_Solver::WriteParticles(int cycle)
{
  if(col->particle_output_is_off())
    return;
  #ifdef NO_HDF5
    eprintf("NO_HDF5 requires OutputMethod=none")
  #else

  timeTasks_set_task(TimeTasks::WRITE_PARTICLES);

  // this is a hack
  convertParticlesToSynched();

  if (col->getWriteMethod() == "Parallel")
  {
    dprintf("pretending to write particles (not yet implemented)");
  }
  else
  {
    fetch_outputWrapperFPP().append_output(
      "position + velocity + q ", cycle, 0);
  }
  #endif // NO_HDF5
}

// This needs to be separated into methods that save particles
// and methods that save field data
//
void c_Solver::WriteOutput(int cycle) {

  // The quickest things should be written first.

  WriteConserved(cycle);
  // This should be invoked by user if desired
  // by means of a callback mechanism.
  //WriteVelocityDistribution(cycle);

  // mechanism to suppress output
  if(!Parameters::get_doWriteOutput())
    return;

  #ifndef NO_HDF5 // array output is only supported for HDF5
  // once phdf5 is properly supported,
  // let's change "Parallel" to mean "phdf5".
  if (col->getWriteMethod() == "H5hut"
   || col->getWriteMethod() == "Parallel")
  {
    /* -------------------------------------------- */
    /* Parallel HDF5 output using the H5hut library */
    /* -------------------------------------------- */

    if (!col->field_output_is_off() && cycle%(col->getFieldOutputCycle())==0)
      WriteFieldsH5hut(ns, grid, EMf, col, vct, cycle);
    if (!col->particle_output_is_off() && cycle%(col->getParticlesOutputCycle())==0)
      WritePartclH5hut(ns, grid, part, col, vct, cycle);
  }
  else if (col->getWriteMethod() == "phdf5")
  {
    /* -------------------------------------------- */
    /* Parallel output using basic hdf5 */
    /* -------------------------------------------- */

    if (!col->field_output_is_off() && cycle%(col->getFieldOutputCycle())==0)
      WriteOutputParallel(grid, EMf, part, col, vct, cycle);
    if (!col->particle_output_is_off() && cycle%(col->getParticlesOutputCycle())==0)
    {
      if(!MPIdata::get_rank())
        warning_printf("WriteParticlesParallel() is not yet implemented.");
      //WritePartclH5hut(ns, grid, part, col, vct, cycle);
    }
  }
  else if (col->getWriteMethod() == "shdf5")
  {
      /*original code
    // write fields-related data
    if (!col->field_output_is_off() && cycle%(col->getFieldOutputCycle())==0)
          WriteFields(cycle);
    // This should be invoked by user if desired
    // by means of a callback mechanism.
    //WriteVirtualSatelliteTraces();

    // write particles-related data
    //
    // this also writes field data...
    
        WriteRestart(cycle);
    
      if (!col->particle_output_is_off() && cycle%(col->getParticlesOutputCycle())==0)
        WriteParticles(cycle);
       */
      if (col->field_output_is_off() || cycle%(col->getFieldOutputCycle())!=0)
          return;

      dprintf("Enter parallel VTK");

      //parallel VTK for B at cells
      int  size[3], subsize[3], start[3], err;
      const int nxc =grid->getNXC(),nyc=grid->getNYC(),nzc=grid->getNZC();
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
      const double firstB = EMf->getBxTot(1,1,1);dprintf("firstB %f", firstB);
      double writebuffer[len];
      int tmpid=0;
	  for(int ix= 0;ix<nxc-2;ix++)
		  for(int iy=0;iy<nyc-2;iy++)
			  for(int iz=0;iz<nzc-2;iz++){
				  tmpid = ix*((nyc-2)*(nzc-2)*3) + iy*((nzc-2)*3) + iz*3;
				  writebuffer[tmpid + 0]=EMf->getBxTot(ix+1,iy+1,iz+1);//(bx[ix+1][iy+1][iz+1]);
				  writebuffer[tmpid + 1]=EMf->getBxTot(ix+1,iy+1,iz+1);//(by[ix+1][iy+1][iz+1]);
				  writebuffer[tmpid + 2]=EMf->getBxTot(ix+1,iy+1,iz+1);//(bz[ix+1][iy+1][iz+1]);//dprintf("writebuffer[ix][iy][iz][2] %i %i %i tmpid%i %f",ix,iy,iz,tmpid,writebuffer[tmpid]);
			  }

	  int TestEndian = 1;
	  int islittleEndian =*(char*)&TestEndian;
	  if(islittleEndian){
		  for(int id= 0;id<len;id++){
			   ByteSwap((unsigned char*) &writebuffer[id],8);
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
      err = MPI_File_write_all(fh, writebuffer,len,MPI_DOUBLE, &status);//test->fetch_arr3()
      int tcount=0;
      MPI_Get_count(&status, testview, &tcount);
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
  else if (col->getWriteMethod() == "default")
  {
    if(col->getParticlesOutputCycle()==1)
    {
      warning_printf(
        "ParticlesOutputCycle=1 now means output particles with every cycle.\n"
        "\tParticlesOutputCycle = 0 turns off particle output.");
    }
    eprintf("The new name for serial hdf5 output is shdf5.\n"
      "\tselect WriteMethod=shdf5.");
  }
  else
  {
    invalid_value_error(col->getWriteMethod().c_str());
  }
  #endif
}

void c_Solver::Finalize() {
  EMf->freeDataType();
  if (col->getCallFinalize())
  {
    #ifndef NO_HDF5
    convertParticlesToSynched();
    writeRESTART(RestartDirName, myrank, (col->getNcycles() + first_cycle) - 1, ns, vct, col, grid, EMf, part, 0);
    #endif
  }

  // stop profiling
  my_clock->stopTiming();
}

//void c_Solver::copyParticlesToSoA()
//{
//  for (int i = 0; i < ns; i++)
//    part[i].copyParticlesToSoA();
//}

void c_Solver::pad_particle_capacities()
{
  for (int i = 0; i < ns; i++)
    part[i].pad_capacities();
}

// convert particle to struct of arrays (assumed by I/O)
void c_Solver::convertParticlesToSoA()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSoA();
}

// convert particle to array of structs (used in computing)
void c_Solver::convertParticlesToAoS()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToAoS();
}

// convert particle to array of structs (used in computing)
void c_Solver::convertParticlesToSynched()
{
  for (int i = 0; i < ns; i++)
    part[i].convertParticlesToSynched();
}

int c_Solver::LastCycle() {
    return (col->getNcycles() + first_cycle);
}
