
#include "mpi.h"
#include "OutputWrapperFPP.h"

#include "CollectiveIO.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"

void OutputWrapperFPP::init_output_files(
    Collective    *col,
    VCtopology3D  *vct,
    Grid3DCU      *grid,
    EMfields3D    *EMf,
    Particles3D   *part,
    int ns)
{
    cartesian_rank = vct->getCartesian_rank();
    stringstream num_proc_ss;
    num_proc_ss << cartesian_rank;
    string num_proc_str = num_proc_ss.str();

    SaveDirName = col->getSaveDirName();
    RestartDirName = col->getRestartDirName();
    int restart_status = col->getRestart_status();
    output_file = SaveDirName + "/proc" + num_proc_str + ".hdf";
    // Initialize the output (simulation results and restart file)
    // PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    // myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
    hdf5_agent.set_simulation_pointers(EMf, grid, vct, col);
    for (int i = 0; i < ns; ++i)
      hdf5_agent.set_simulation_pointers_part(&part[i]);
    // Add the HDF5 output agent to the Output Manager's list
    output_mgr.push_back(&hdf5_agent);
    // I change "&" to "&&" in the following line. -eaj
    if (cartesian_rank == 0 && restart_status < 2) {
      hdf5_agent.open(SaveDirName + "/settings.hdf");
      output_mgr.output("collective + total_topology + proc_topology", 0);
      hdf5_agent.close();
      hdf5_agent.open(RestartDirName + "/settings.hdf");
      output_mgr.output("collective + total_topology + proc_topology", 0);
      hdf5_agent.close();
    }
    // Restart
    if (restart_status == 0) {           // new simulation from input file
      hdf5_agent.open(output_file);
      output_mgr.output("proc_topology ", 0);
      hdf5_agent.close();
    }
    // restart append the results to the previous simulation 
    else { 
      hdf5_agent.open_append(output_file);
      output_mgr.output("proc_topology ", 0);
      hdf5_agent.close();
    }
}

void OutputWrapperFPP::append_output(
  const char* tag, int cycle)
{
    hdf5_agent.open_append(output_file);
    output_mgr.output(tag, cycle);
    hdf5_agent.close();
}

void OutputWrapperFPP::append_output(
  const char* tag, int cycle, int sample)
{
    hdf5_agent.open_append(output_file);
    output_mgr.output(tag, cycle, sample);
    // Pressure tensor is available
    hdf5_agent.close();
}
