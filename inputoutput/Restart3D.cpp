#include <mpi.h>
#include "Restart3D.h"
#include "CollectiveIO.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "Particles3Dcomm.h"
#include "PSKOutput.h"
#include "PSKhdf5adaptor.h"
using std::stringstream;

// === methods to write restart files ===

/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;

  printf("%s/restart%s.hdf\n", SaveDirName.c_str(), ss.str().c_str());
  //cout << SaveDirName + "/restart" + ss.str() + ".hdf" << endl;

  hdf5_agent.open_append(SaveDirName + "/restart" + ss.str().c_str() + ".hdf");
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Eall + Ball + rhos", cycle);
  output_mgr.output("position + velocity + q ", cycle, 0);
  hdf5_agent.close();

}

/** this restart function writes the last restart with the last cycle */
void writeRESTART(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part, bool fool) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;
  hdf5_agent.open(SaveDirName + "/restart" + ss.str() + ".hdf");
  output_mgr.output("proc_topology ", cycle);
  output_mgr.output("Eall + Ball + rhos", cycle);
  output_mgr.output("position + velocity + q ", cycle, 0);
  output_mgr.output("last_cycle", cycle);
  hdf5_agent.close();

}


/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART_ES(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;

  printf("%s/restart%s.hdf\n", SaveDirName.c_str(), ss.str().c_str());
  //cout << SaveDirName + "/restart" + ss.str() + ".hdf" << endl;

  hdf5_agent.open_append(SaveDirName + "/restart" + ss.str().c_str() + ".hdf");
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Ex + rhos", cycle);
  output_mgr.output("x + u + q ", cycle, 0);
  hdf5_agent.close();

}

/** this restart function writes the last restart with the last cycle */
void writeRESTART_ES(const string& SaveDirName, int myrank, int cycle, int ns, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part, bool fool) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;
  hdf5_agent.open(SaveDirName + "/restart" + ss.str() + ".hdf");
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Ex + rhos", 0);
  output_mgr.output("x + u + q ", 0, 0);
  output_mgr.output("last_cycle", cycle);
  hdf5_agent.close();

}

// === methods to read restart files ===

// extracted from EMfields3D.cpp
//
void read_field_restart(
    const Collective* col,
    const VCtopology3D* vct,
    const Grid* grid,
    arr3_double Bxn, arr3_double Byn, arr3_double Bzn,
    arr3_double Ex, arr3_double Ey, arr3_double Ez,
    array4_double* rhons_, int ns)
{
    const int nxn = grid->getNXN();
    const int nyn = grid->getNYN();
    const int nzn = grid->getNZN();
    if (vct->getCartesian_rank() == 0)
    {
      printf("LOADING EM FIELD FROM RESTART FILE in %s/restart.hdf\n",
        col->getRestartDirName().c_str());
      //cout << "LOADING EM FIELD FROM RESTART FILE in " + col->getRestartDirName() + "/restart.hdf" << endl;
    }
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = col->getRestartDirName() + "/restart" + ss.str() + ".hdf";
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[3];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      eprintf("couldn't open file: %s\n"
        "\tRESTART NOT POSSIBLE", name_file.c_str());
      //cout << "couldn't open file: " << name_file << endl;
      //cout << "RESTART NOT POSSIBLE" << endl;
    }

    dataset_id = H5Dopen2(file_id, "/fields/Bx/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);



    // Bxn
    double *temp_storage = new double[dims_out[0] * dims_out[1] * dims_out[2]];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bxn[i][j][jj] = temp_storage[k++];


    status = H5Dclose(dataset_id);

    // Byn
    dataset_id = H5Dopen2(file_id, "/fields/By/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Byn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Bzn
    dataset_id = H5Dopen2(file_id, "/fields/Bz/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bzn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ex
    dataset_id = H5Dopen2(file_id, "/fields/Ex/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ex[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ey 
    dataset_id = H5Dopen2(file_id, "/fields/Ey/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ey[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // Ez 
    dataset_id = H5Dopen2(file_id, "/fields/Ez/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ez[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // open the charge density for species

    stringstream *species_name = new stringstream[ns];
    for (int is = 0; is < ns; is++) {
      species_name[is] << is;
      string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
      dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
      status = H5Dclose(dataset_id);
      array4_double& rhons = *rhons_;
      k = 0;
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int jj = 1; jj < nzn - 1; jj++)
            rhons[is][i][j][jj] = temp_storage[k++];
    }

    // close the hdf file
    status = H5Fclose(file_id);
    delete[]temp_storage;
    delete[]species_name;
}

// extracted from Particles3Dcomm.cpp
//
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
    vector_double& t)
{
    if (vct->getCartesian_rank() == 0 && species_number == 0)
    {
      printf("LOADING PARTICLES FROM RESTART FILE in %s/restart.hdf\n",
        col->getRestartDirName().c_str());
      //cout << "LOADING PARTICLES FROM RESTART FILE in " + col->getRestartDirName() + "/restart.hdf" << endl;
    }
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = col->getRestartDirName() + "/restart" + ss.str() + ".hdf";
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[1];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      eprintf("couldn't open file: %s\n"
        "\tRESTART NOT POSSIBLE", name_file.c_str());
      //cout << "couldn't open file: " << name_file << endl;
      //cout << "RESTART NOT POSSIBLE" << endl;
    }

    stringstream species_name;
    species_name << species_number;
    // the cycle of the last restart is set to 0
    string name_dataset = "/particles/species_" + species_name.str() + "/x/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id); /* dataspace handle */
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    // get how many particles there are on this processor for this species
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    const int nop = dims_out[0]; // number of particles in this process
    //Particles3Dcomm::resize_SoA(nop);
    {
      //
      // allocate space for particles including padding
      //
      const int padded_nop = roundup_to_multiple(nop,DVECWIDTH);
      u.reserve(padded_nop);
      v.reserve(padded_nop);
      w.reserve(padded_nop);
      q.reserve(padded_nop);
      x.reserve(padded_nop);
      y.reserve(padded_nop);
      z.reserve(padded_nop);
      t.reserve(padded_nop);
      //
      // define size of particle data
      //
      u.resize(nop);
      v.resize(nop);
      w.resize(nop);
      q.resize(nop);
      x.resize(nop);
      y.resize(nop);
      z.resize(nop);
      t.resize(nop);
    }
    // get x
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x[0]);
    // close the data set
    status = H5Dclose(dataset_id);

    // get y
    name_dataset = "/particles/species_" + species_name.str() + "/y/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y[0]);
    status = H5Dclose(dataset_id);

    // get z
    name_dataset = "/particles/species_" + species_name.str() + "/z/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z[0]);
    status = H5Dclose(dataset_id);

    // get u
    name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &u[0]);
    status = H5Dclose(dataset_id);
    // get v
    name_dataset = "/particles/species_" + species_name.str() + "/v/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0]);
    status = H5Dclose(dataset_id);
    // get w
    name_dataset = "/particles/species_" + species_name.str() + "/w/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &w[0]);
    status = H5Dclose(dataset_id);
    // get q
    name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &q[0]);
    status = H5Dclose(dataset_id);
    // ID 
    //if (TrackParticleID) {
    //  // herr_t (*old_func)(void*); // HDF 1.6
    //  H5E_auto2_t old_func;      // HDF 1.8.8
    //  void *old_client_data;
    //  H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);  // HDF 1.8.8
    //  /* Turn off error handling */
    //  // H5Eset_auto(NULL, NULL); // HDF 1.6
    //  H5Eset_auto2(H5E_DEFAULT, 0, 0); // HDF 1.8
    //  name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
    //  dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    //
    //  // H5Eset_auto(old_func, old_client_data); // HDF 1.6
    //  H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    //  if (dataset_id > 0)
    //    status = H5Dread(dataset_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParticleID);
    //  //else {
    //  //  for (int counter = 0; counter < nop; counter++)
    //  //    fetch_ParticleID(counter) = particleIDgenerator.get_ID();
    //  //}
    //}
    // close the hdf file
    status = H5Fclose(file_id);
}
