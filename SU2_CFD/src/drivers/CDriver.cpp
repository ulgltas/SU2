/*!
 * \file driver_structure.cpp
 * \brief The main subroutines for driving single or multi-zone problems.
 * \author T. Economon, H. Kline, R. Sanchez, F. Palacios
 * \version 7.0.2 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/drivers/CDriver.hpp"
#include "../../include/definition_structure.hpp"

#include "../../../Common/include/geometry/CDummyGeometry.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../Common/include/geometry/CMultiGridGeometry.hpp"

#include "../../include/solvers/CEulerSolver.hpp"
#include "../../include/solvers/CIncEulerSolver.hpp"
#include "../../include/solvers/CNSSolver.hpp"
#include "../../include/solvers/CIncNSSolver.hpp"
#include "../../include/solvers/CTurbSASolver.hpp"
#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../include/solvers/CTransLMSolver.hpp"
#include "../../include/solvers/CAdjEulerSolver.hpp"
#include "../../include/solvers/CAdjNSSolver.hpp"
#include "../../include/solvers/CAdjTurbSolver.hpp"
#include "../../include/solvers/CHeatSolver.hpp"
#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/solvers/CTemplateSolver.hpp"
#include "../../include/solvers/CDiscAdjSolver.hpp"
#include "../../include/solvers/CDiscAdjFEASolver.hpp"
#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"
#include "../../include/solvers/CFEM_DG_NSSolver.hpp"
#include "../../include/solvers/CDiscAdjMeshSolver.hpp"
#include "../../include/solvers/CMeshSolver.hpp"
#include "../../include/solvers/CRadSolver.hpp"
#include "../../include/solvers/CRadP1Solver.hpp"

#include "../../include/interfaces/cfd/CConservativeVarsInterface.hpp"
#include "../../include/interfaces/cfd/CMixingPlaneInterface.hpp"
#include "../../include/interfaces/cfd/CSlidingInterface.hpp"
#include "../../include/interfaces/cht/CConjugateHeatInterface.hpp"
#include "../../include/interfaces/fsi/CDisplacementsInterface.hpp"
#include "../../include/interfaces/fsi/CFlowTractionInterface.hpp"
#include "../../include/interfaces/fsi/CDiscAdjFlowTractionInterface.hpp"
#include "../../include/interfaces/fsi/CDisplacementsInterfaceLegacy.hpp"
#include "../../include/interfaces/fsi/CDiscAdjDisplacementsInterfaceLegacy.hpp"

#include "../../include/numerics/template.hpp"
#include "../../include/numerics/transition.hpp"
#include "../../include/numerics/radiation.hpp"
#include "../../include/numerics/heat.hpp"
#include "../../include/numerics/flow/convection/roe.hpp"
#include "../../include/numerics/flow/convection/fds.hpp"
#include "../../include/numerics/flow/convection/fvs.hpp"
#include "../../include/numerics/flow/convection/cusp.hpp"
#include "../../include/numerics/flow/convection/hllc.hpp"
#include "../../include/numerics/flow/convection/ausm_slau.hpp"
#include "../../include/numerics/flow/convection/centered.hpp"
#include "../../include/numerics/flow/flow_diffusion.hpp"
#include "../../include/numerics/flow/flow_sources.hpp"
#include "../../include/numerics/continuous_adjoint/adj_convection.hpp"
#include "../../include/numerics/continuous_adjoint/adj_diffusion.hpp"
#include "../../include/numerics/continuous_adjoint/adj_sources.hpp"
#include "../../include/numerics/turbulent/turb_convection.hpp"
#include "../../include/numerics/turbulent/turb_diffusion.hpp"
#include "../../include/numerics/turbulent/turb_sources.hpp"
#include "../../include/numerics/elasticity/CFEAElasticity.hpp"
#include "../../include/numerics/elasticity/CFEALinearElasticity.hpp"
#include "../../include/numerics/elasticity/CFEANonlinearElasticity.hpp"
#include "../../include/numerics/elasticity/nonlinear_models.hpp"

#include "../../include/integration/CSingleGridIntegration.hpp"
#include "../../include/integration/CMultiGridIntegration.hpp"
#include "../../include/integration/CStructuralIntegration.hpp"
#include "../../include/integration/CFEM_DG_Integration.hpp"

#include "../../../Common/include/omp_structure.hpp"

#include <cassert>

#ifdef VTUNEPROF
#include <ittnotify.h>
#endif
#include <fenv.h>

CDriver::CDriver(char* confFile,
                 unsigned short val_nZone,
                 SU2_Comm MPICommunicator, bool dummy_geo):config_file_name(confFile), StartTime(0.0), StopTime(0.0), UsedTime(0.0),
                 TimeIter(0), nZone(val_nZone), StopCalc(false), fsi(false), fem_solver(false), dry_run(dummy_geo) {

  /*--- Initialize Medipack (must also be here so it is initialized from python) ---*/
#ifdef HAVE_MPI
  #if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
    SU2_MPI::Init_AMPI();
  #endif
#endif

  unsigned short jZone;

  SU2_MPI::SetComm(MPICommunicator);

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  /*--- Start timer to track preprocessing for benchmarking. ---*/

#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif

  /*--- Initialize containers with null --- */

  SetContainers_Null();

  /*--- Preprocessing of the config files. In this routine, the config file is read
   and it is determined whether a problem is single physics or multiphysics. . ---*/

  Input_Preprocessing(config_container, driver_config);

  /*--- Retrieve dimension from mesh file ---*/

  nDim = CConfig::GetnDim(config_container[ZONE_0]->GetMesh_FileName(),
                          config_container[ZONE_0]->GetMesh_FileFormat());

  /*--- Output preprocessing ---*/

  Output_Preprocessing(config_container, driver_config, output_container, driver_output);


  for (iZone = 0; iZone < nZone; iZone++) {

    /*--- Read the number of instances for each zone ---*/

    nInst[iZone] = config_container[iZone]->GetnTimeInstances();

    geometry_container[iZone]    = new CGeometry**    [nInst[iZone]];
    iteration_container[iZone]   = new CIteration*    [nInst[iZone]];
    solver_container[iZone]      = new CSolver***     [nInst[iZone]];
    integration_container[iZone] = new CIntegration** [nInst[iZone]];
    numerics_container[iZone]    = new CNumerics****  [nInst[iZone]];
    grid_movement[iZone]         = new CVolumetricMovement* [nInst[iZone]];

    /*--- Allocate transfer and interpolation container --- */

    interface_container[iZone]    = new CInterface*[nZone];
    interpolator_container[iZone] = new CInterpolator*[nZone];

    for (jZone = 0; jZone < nZone; jZone++){
      interface_container[iZone][jZone]            = NULL;
      interpolator_container[iZone][jZone]         = NULL;
    }

    for (iInst = 0; iInst < nInst[iZone]; iInst++){

      config_container[iZone]->SetiInst(iInst);

      geometry_container[iZone][iInst]    = NULL;
      iteration_container[iZone][iInst]   = NULL;
      solver_container[iZone][iInst]      = NULL;
      integration_container[iZone][iInst] = NULL;
      grid_movement[iZone][iInst]         = NULL;

      /*--- Preprocessing of the geometry for all zones. In this routine, the edge-
       based data structure is constructed, i.e. node and cell neighbors are
       identified and linked, face areas and volumes of the dual mesh cells are
       computed, and the multigrid levels are created using an agglomeration procedure. ---*/

      Geometrical_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], dry_run);

      /*--- Definition of the solver class: solver_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS].
       The solver classes are specific to a particular set of governing equations,
       and they contain the subroutines with instructions for computing each spatial
       term of the PDE, i.e. loops over the edges to compute convective and viscous
       fluxes, loops over the nodes to compute source terms, and routines for
       imposing various boundary condition type for the PDE. ---*/

      Solver_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], solver_container[iZone][iInst]);

      /*--- Definition of the numerical method class:
       numerics_container[#ZONES][#INSTANCES][#MG_GRIDS][#EQ_SYSTEMS][#EQ_TERMS].
       The numerics class contains the implementation of the numerical methods for
       evaluating convective or viscous fluxes between any two nodes in the edge-based
       data structure (centered, upwind, galerkin), as well as any source terms
       (piecewise constant reconstruction) evaluated in each dual mesh volume. ---*/

      Numerics_Preprocessing(config_container[iZone], geometry_container[iZone][iInst],
                             solver_container[iZone][iInst], numerics_container[iZone][iInst]);

      /*--- Definition of the integration class: integration_container[#ZONES][#INSTANCES][#EQ_SYSTEMS].
       The integration class orchestrates the execution of the spatial integration
       subroutines contained in the solver class (including multigrid) for computing
       the residual at each node, R(U) and then integrates the equations to a
       steady state or time-accurately. ---*/

      Integration_Preprocessing(config_container[iZone], integration_container[iZone][iInst]);

      /*--- Instantiate the type of physics iteration to be executed within each zone. For
       example, one can execute the same physics across multiple zones (mixing plane),
       different physics in different zones (fluid-structure interaction), or couple multiple
       systems tightly within a single zone by creating a new iteration class (e.g., RANS). ---*/

      Iteration_Preprocessing(config_container[iZone], iteration_container[iZone][iInst]);

      /*--- Dynamic mesh processing.  ---*/

      DynamicMesh_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], solver_container[iZone][iInst],
                                 iteration_container[iZone][iInst], grid_movement[iZone][iInst], surface_movement[iZone]);
      /*--- Static mesh processing.  ---*/

      StaticMesh_Preprocessing(config_container[iZone], geometry_container[iZone][iInst], surface_movement[iZone]);

    }

  }

  /*--- Definition of the interface and transfer conditions between different zones.
   *--- The transfer container is defined for zones paired one to one.
   *--- This only works for a multizone FSI problem (nZone > 1).
   *--- Also, at the moment this capability is limited to two zones (nZone < 3).
   *--- This will change in the future. ---*/

  if ( nZone > 1 ) {
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Multizone Interface Preprocessing -------------------" << endl;

    Interface_Preprocessing(config_container, solver_container, geometry_container,
                            interface_types, interface_container, interpolator_container);
  }

  if(fsi && (config_container[ZONE_0]->GetRestart() || config_container[ZONE_0]->GetDiscrete_Adjoint())){
    if (rank == MASTER_NODE)cout << endl <<"Restarting Fluid and Structural Solvers." << endl;

    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        Solver_Restart(solver_container[iZone][iInst], geometry_container[iZone][iInst],
                       config_container[iZone], true);
      }
    }
  }

  if (config_container[ZONE_0]->GetBoolTurbomachinery()){
    if (rank == MASTER_NODE)cout << endl <<"---------------------- Turbomachinery Preprocessing ---------------------" << endl;
    Turbomachinery_Preprocessing(config_container, geometry_container, solver_container, interface_container);
  }


  PythonInterface_Preprocessing(config_container, geometry_container, solver_container);

  /*--- Open the FSI convergence history file ---*/

//  if (fsi){
//      if (rank == MASTER_NODE) cout << endl <<"Opening FSI history file." << endl;
//      unsigned short ZONE_FLOW = 0, ZONE_STRUCT = 1;
//      output_container->SpecialOutput_FSI(&FSIHist_file, geometry_container, solver_container,
//                                config_container, integration_container, 0,
//                                ZONE_FLOW, ZONE_STRUCT, true);
//  }

  /*--- Preprocessing time is reported now, but not included in the next compute portion. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc    = UsedTime;
  UsedTimeCompute    = 0.0;
  UsedTimeOutput     = 0.0;
  IterCount          = 0;
  OutputCount        = 0;
  MDOFs              = 0.0;
  MDOFsDomain        = 0.0;
  Mpoints            = 0.0;
  MpointsDomain      = 0.0;
  for (iZone = 0; iZone < nZone; iZone++) {
    Mpoints       +=(su2double)geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPoint()/(1.0e6);
    MpointsDomain +=(su2double)geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPointDomain()/(1.0e6);
    MDOFs         += (su2double)DOFsPerPoint*(su2double)geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPoint()/(1.0e6);
    MDOFsDomain   += (su2double)DOFsPerPoint*(su2double)geometry_container[iZone][INST_0][MESH_0]->GetGlobal_nPointDomain()/(1.0e6);
  }

  /*--- Reset timer for compute/output performance benchmarking. ---*/
#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif

  /*--- Compute/print the total time for performance benchmarking. ---*/

  UsedTime = StopTime-StartTime;
  UsedTimePreproc = UsedTime;

  /*--- Reset timer for compute performance benchmarking. ---*/
#ifndef HAVE_MPI
  StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StartTime = MPI_Wtime();
#endif

}

void CDriver::SetContainers_Null(){

  /*--- Create pointers to all of the classes that may be used throughout
   the SU2_CFD code. In general, the pointers are instantiated down a
   hierarchy over all zones, multigrid levels, equation sets, and equation
   terms as described in the comments below. ---*/

  ConvHist_file                  = NULL;
  iteration_container            = NULL;
  output_container               = NULL;
  integration_container          = NULL;
  geometry_container             = NULL;
  solver_container               = NULL;
  numerics_container             = NULL;
  config_container               = NULL;
  surface_movement               = NULL;
  grid_movement                  = NULL;
  FFDBox                         = NULL;
  interpolator_container         = NULL;
  interface_container            = NULL;
  interface_types                = NULL;
  nInst                          = NULL;


  /*--- Definition and of the containers for all possible zones. ---*/

  iteration_container            = new CIteration**[nZone];
  solver_container               = new CSolver****[nZone];
  integration_container          = new CIntegration***[nZone];
  numerics_container             = new CNumerics*****[nZone];
  config_container               = new CConfig*[nZone];
  geometry_container             = new CGeometry***[nZone];
  surface_movement               = new CSurfaceMovement*[nZone];
  grid_movement                  = new CVolumetricMovement**[nZone];
  FFDBox                         = new CFreeFormDefBox**[nZone];
  interpolator_container         = new CInterpolator**[nZone];
  interface_container            = new CInterface**[nZone];
  interface_types                = new unsigned short*[nZone];
  output_container               = new COutput*[nZone];
  nInst                          = new unsigned short[nZone];
  driver_config                  = NULL;
  driver_output                  = NULL;


  for (iZone = 0; iZone < nZone; iZone++) {
    solver_container[iZone]               = NULL;
    integration_container[iZone]          = NULL;
    numerics_container[iZone]             = NULL;
    config_container[iZone]               = NULL;
    geometry_container[iZone]             = NULL;
    surface_movement[iZone]               = NULL;
    grid_movement[iZone]                  = NULL;
    FFDBox[iZone]                         = NULL;
    interpolator_container[iZone]         = NULL;
    interface_container[iZone]            = NULL;
    interface_types[iZone]                = new unsigned short[nZone];
    output_container[iZone]               = NULL;
    nInst[iZone]                          = 1;
  }

  strcpy(runtime_file_name, "runtime.dat");

}


void CDriver::Postprocessing() {

  bool isBinary = config_container[ZONE_0]->GetWrt_Binary_Restart();
  bool wrt_perf = config_container[ZONE_0]->GetWrt_Performance();

    /*--- Output some information to the console. ---*/

  if (rank == MASTER_NODE) {

    /*--- Print out the number of non-physical points and reconstructions ---*/

    if (config_container[ZONE_0]->GetNonphysical_Points() > 0)
      cout << "Warning: there are " << config_container[ZONE_0]->GetNonphysical_Points() << " non-physical points in the solution." << endl;
    if (config_container[ZONE_0]->GetNonphysical_Reconstr() > 0)
      cout << "Warning: " << config_container[ZONE_0]->GetNonphysical_Reconstr() << " reconstructed states for upwinding are non-physical." << endl;
  }

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Numerics_Postprocessing(numerics_container[iZone], solver_container[iZone][iInst],
          geometry_container[iZone][iInst], config_container[iZone], iInst);
    }
    delete [] numerics_container[iZone];
  }
  delete [] numerics_container;
  if (rank == MASTER_NODE) cout << "Deleted CNumerics container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Integration_Postprocessing(integration_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    }
    delete [] integration_container[iZone];
  }
  delete [] integration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIntegration container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      Solver_Postprocessing(solver_container[iZone],
          geometry_container[iZone][iInst],
          config_container[iZone],
          iInst);
    }
    delete [] solver_container[iZone];
  }
  delete [] solver_container;
  if (rank == MASTER_NODE) cout << "Deleted CSolver container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      delete iteration_container[iZone][iInst];
    delete [] iteration_container[iZone];
  }
  delete [] iteration_container;
  if (rank == MASTER_NODE) cout << "Deleted CIteration container." << endl;

  if (interpolator_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (interpolator_container[iZone] != NULL) {
        for (unsigned short jZone = 0; jZone < nZone; jZone++)
          if (interpolator_container[iZone][jZone] != NULL)
            delete interpolator_container[iZone][jZone];
        delete [] interpolator_container[iZone];
      }
    }
    delete [] interpolator_container;
    if (rank == MASTER_NODE) cout << "Deleted CInterpolator container." << endl;
  }

  if (interface_container != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (interface_container[iZone] != NULL) {
        for (unsigned short jZone = 0; jZone < nZone; jZone++)
          if (interface_container[iZone][jZone] != NULL)
            delete interface_container[iZone][jZone];
        delete [] interface_container[iZone];
      }
    }
    delete [] interface_container;
    if (rank == MASTER_NODE) cout << "Deleted CInterface container." << endl;
  }

  if (interface_types != NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (interface_types[iZone] != NULL)
      delete [] interface_types[iZone];
    }
    delete [] interface_types;
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    if (geometry_container[iZone] != NULL) {
      for (iInst = 0; iInst < nInst[iZone]; iInst++){
        for (unsigned short iMGlevel = 0; iMGlevel < config_container[iZone]->GetnMGLevels()+1; iMGlevel++) {
          if (geometry_container[iZone][iInst][iMGlevel] != NULL) delete geometry_container[iZone][iInst][iMGlevel];
        }
        if (geometry_container[iZone][iInst] != NULL) delete [] geometry_container[iZone][iInst];
      }
      delete [] geometry_container[iZone];
    }
  }
  delete [] geometry_container;
  if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete [] FFDBox[iZone];
  }
  delete [] FFDBox;
  if (rank == MASTER_NODE) cout << "Deleted CFreeFormDefBox class." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    delete surface_movement[iZone];
  }
  delete [] surface_movement;
  if (rank == MASTER_NODE) cout << "Deleted CSurfaceMovement class." << endl;

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++){
      if (grid_movement[iZone][iInst] != NULL) delete grid_movement[iZone][iInst];
    }
    if (grid_movement[iZone] != NULL) delete [] grid_movement[iZone];
  }
  delete [] grid_movement;
  if (rank == MASTER_NODE) cout << "Deleted CVolumetricMovement class." << endl;

  /*--- Output profiling information ---*/
  // Note that for now this is called only by a single thread, but all
  // necessary variables have been made thread private for safety (tick/tock)!!

  config_container[ZONE_0]->SetProfilingCSV();
  config_container[ZONE_0]->GEMMProfilingCSV();

  /*--- Deallocate config container ---*/
  if (config_container!= NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (config_container[iZone] != NULL) {
        delete config_container[iZone];
      }
    }
    delete [] config_container;
  }
  if (driver_config != NULL) delete driver_config;
  if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;

  if (nInst != NULL) delete [] nInst;
  if (rank == MASTER_NODE) cout << "Deleted nInst container." << endl;

  /*--- Deallocate output container ---*/

  if (output_container!= NULL) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if (output_container[iZone] != NULL) {
        delete output_container[iZone];
      }
    }
    delete [] output_container;
  }

  if(driver_output != NULL){
    delete driver_output;
  }

  if (rank == MASTER_NODE) cout << "Deleted COutput class." << endl;

  if (rank == MASTER_NODE) cout << "-------------------------------------------------------------------------" << endl;


  /*--- Stop the timer and output the final performance summary. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  UsedTime = StopTime-StartTime;
  UsedTimeCompute += UsedTime;

  if ((rank == MASTER_NODE) && (wrt_perf)) {
    su2double TotalTime = UsedTimePreproc + UsedTimeCompute + UsedTimeOutput;
    cout.precision(6);
    cout << endl << endl <<"-------------------------- Performance Summary --------------------------" << endl;
    cout << "Simulation totals:" << endl;
    cout << setw(25) << "Wall-clock time (hrs):" << setw(12) << (TotalTime)/(60.0*60.0) << " | ";
    cout << setw(20) << "Core-hrs:" << setw(12) << (su2double)size*(TotalTime)/(60.0*60.0) << endl;
    cout << setw(25) << "Cores:" << setw(12) << size << " | ";
    cout << setw(20) << "DOFs/point:" << setw(12) << (su2double)DOFsPerPoint << endl;
    cout << setw(25) << "Points/core:" << setw(12) << 1.0e6*MpointsDomain/(su2double)size << " | ";
    cout << setw(20) << "Ghost points/core:" << setw(12) << 1.0e6*(Mpoints-MpointsDomain)/(su2double)size << endl;
    cout << setw(25) << "Ghost/Owned Point Ratio:" << setw(12) << (Mpoints-MpointsDomain)/MpointsDomain << " | " << endl;
    cout << endl;
    cout << "Preprocessing phase:" << endl;
    cout << setw(25) << "Preproc. Time (s):"  << setw(12)<< UsedTimePreproc << " | ";
    cout << setw(20) << "Preproc. Time (%):" << setw(12)<< ((UsedTimePreproc * 100.0) / (TotalTime)) << endl;
    cout << endl;
    cout << "Compute phase:" << endl;
    cout << setw(25) << "Compute Time (s):"  << setw(12)<< UsedTimeCompute << " | ";
    cout << setw(20) << "Compute Time (%):" << setw(12)<< ((UsedTimeCompute * 100.0) / (TotalTime)) << endl;
    cout << setw(25) << "Iteration count:"  << setw(12)<< IterCount << " | ";
    if (IterCount != 0) {
      cout << setw(20) << "Avg. s/iter:" << setw(12)<< UsedTimeCompute/(su2double)IterCount << endl;
      cout << setw(25) << "Core-s/iter/Mpoints:" << setw(12)<< (su2double)size*UsedTimeCompute/(su2double)IterCount/Mpoints << " | ";
      cout << setw(20) << "Mpoints/s:" << setw(12)<< Mpoints*(su2double)IterCount/UsedTimeCompute << endl;
    } else cout << endl;
    cout << endl;
    cout << "Output phase:" << endl;
    cout << setw(25) << "Output Time (s):"  << setw(12)<< UsedTimeOutput << " | ";
    cout << setw(20) << "Output Time (%):" << setw(12)<< ((UsedTimeOutput * 100.0) / (TotalTime)) << endl;
    cout << setw(25) << "Output count:" << setw(12)<< OutputCount << " | ";
    if (OutputCount != 0) {
      cout << setw(20)<< "Avg. s/output:" << setw(12)<< UsedTimeOutput/(su2double)OutputCount << endl;
      if (isBinary) {
        cout << setw(25)<< "Restart Aggr. BW (MB/s):" << setw(12)<< BandwidthSum/(su2double)OutputCount << " | ";
        cout << setw(20)<< "MB/s/core:" << setw(12)<< BandwidthSum/(su2double)OutputCount/(su2double)size << endl;
      }
    } else cout << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    cout << endl;
  }

  /*--- Exit the solver cleanly ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------- Exit Success (SU2_CFD) ------------------------" << endl << endl;

}


void CDriver::Input_Preprocessing(CConfig **&config, CConfig *&driver_config) {

  char zone_file_name[MAX_STRING_SIZE];

  /*--- Initialize the configuration of the driver ---*/

  driver_config = new CConfig(config_file_name, SU2_CFD, false);

  for (iZone = 0; iZone < nZone; iZone++) {

    if (rank == MASTER_NODE){
      cout  << endl << "Parsing config file for zone " << iZone << endl;
    }
    /*--- Definition of the configuration option class for all zones. In this
     constructor, the input configuration file is parsed and all options are
     read and stored. ---*/

    if (driver_config->GetnConfigFiles() > 0){

      strcpy(zone_file_name, driver_config->GetConfigFilename(iZone).c_str());
      config[iZone] = new CConfig(driver_config, zone_file_name, SU2_CFD, iZone, nZone, true);
    }
    else{
      config[iZone] = new CConfig(driver_config, config_file_name, SU2_CFD, iZone, nZone, true);
    }

    /*--- Set the MPI communicator ---*/

    config[iZone]->SetMPICommunicator(SU2_MPI::GetComm());
  }


  /*--- Set the multizone part of the problem. ---*/
  if (driver_config->GetMultizone_Problem()){
    for (iZone = 0; iZone < nZone; iZone++) {
      /*--- Set the interface markers for multizone ---*/
      config_container[iZone]->SetMultizone(driver_config, config_container);
    }
  }

  /*--- Determine whether or not the FEM solver is used, which decides the
 type of geometry classes that are instantiated. Only adapted for single-zone problems ---*/

  fem_solver = ((config_container[ZONE_0]->GetKind_Solver() == FEM_EULER)          ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_NAVIER_STOKES)  ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_RANS)           ||
                (config_container[ZONE_0]->GetKind_Solver() == FEM_LES)            ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_EULER) ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_NS)    ||
                (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_FEM_RANS));

  fsi = config_container[ZONE_0]->GetFSI_Simulation();
}

void CDriver::Geometrical_Preprocessing(CConfig* config, CGeometry **&geometry, bool dummy){

  if (!dummy){
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Geometry Preprocessing ( Zone " << config->GetiZone() <<" ) -------------------" << endl;

    if( fem_solver ) {
      switch( config->GetKind_FEM_Flow() ) {
        case DG: {
            Geometrical_Preprocessing_DGFEM(config, geometry);
            break;
          }
      }
    }
    else {
      Geometrical_Preprocessing_FVM(config, geometry);
    }
  } else {
    if (rank == MASTER_NODE)
      cout << endl <<"-------------------------- Using Dummy Geometry -------------------------" << endl;

    unsigned short iMGlevel;

    geometry = new CGeometry*[config->GetnMGLevels()+1];

    if (!fem_solver){
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        geometry[iMGlevel] = new CDummyGeometry(config);
      }
    } else {
      geometry[ZONE_0] = new CDummyMeshFEM_DG(config);
    }

    nDim = geometry[ZONE_0]->GetnDim();
  }

  /*--- Computation of wall distances for turbulence modeling ---*/

  if ((config->GetKind_Solver() == RANS) ||
      (config->GetKind_Solver() == INC_RANS) ||
      (config->GetKind_Solver() == ADJ_RANS) ||
      (config->GetKind_Solver() == DISC_ADJ_INC_RANS) ||
      (config->GetKind_Solver() == DISC_ADJ_RANS) ||
      (config->GetKind_Solver() == FEM_RANS) ||
      (config->GetKind_Solver() == FEM_LES) ) {

    if (rank == MASTER_NODE)
      cout << "Computing wall distances." << endl;

    geometry[MESH_0]->ComputeWall_Distance(config);
  }

  /*--- Computation of positive surface area in the z-plane which is used for
     the calculation of force coefficient (non-dimensionalization). ---*/

  geometry[MESH_0]->SetPositive_ZArea(config);

  /*--- Set the near-field, interface and actuator disk boundary conditions, if necessary. ---*/

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    geometry[iMesh]->MatchNearField(config);
    geometry[iMesh]->MatchActuator_Disk(config);
  }

  /*--- If we have any periodic markers in this calculation, we must
       match the periodic points found on both sides of the periodic BC.
       Note that the current implementation requires a 1-to-1 matching of
       periodic points on the pair of periodic faces after the translation
       or rotation is taken into account. ---*/

  if ((config->GetnMarker_Periodic() != 0) && !fem_solver) {
    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

      /*--- Note that we loop over pairs of periodic markers individually
           so that repeated nodes on adjacent periodic faces are properly
           accounted for in multiple places. ---*/

      for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
        geometry[iMesh]->MatchPeriodic(config, iPeriodic);
      }

      /*--- Initialize the communication framework for the periodic BCs. ---*/
      geometry[iMesh]->PreprocessPeriodicComms(geometry[iMesh], config);

    }
  }

  /*--- If activated by the compile directive, perform a partition analysis. ---*/
#if PARTITION
  if (!dummy){
    if( fem_solver ) Partition_Analysis_FEM(geometry[MESH_0], config);
    else Partition_Analysis(geometry[MESH_0], config);
  }
#endif

  /*--- Check if Euler & Symmetry markers are straight/plane. This information
        is used in the Euler & Symmetry boundary routines. ---*/
  if((config_container[iZone]->GetnMarker_Euler() != 0 ||
     config_container[iZone]->GetnMarker_SymWall() != 0) &&
     !fem_solver) {

    if (rank == MASTER_NODE)
      cout << "Checking if Euler & Symmetry markers are straight/plane:" << endl;

    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++)
      geometry_container[iZone][iInst][iMesh]->ComputeSurf_Straightness(config_container[iZone], (iMesh==MESH_0) );

  }

}

void CDriver::Geometrical_Preprocessing_FVM(CConfig *config, CGeometry **&geometry) {

  unsigned short iZone = config->GetiZone(), iMGlevel;
  unsigned short requestedMGlevels = config->GetnMGLevels();
  unsigned long iPoint;
  bool fea = false;

  /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/

  CGeometry *geometry_aux = NULL;

  /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

  geometry_aux = new CPhysicalGeometry(config, iZone, nZone);

  /*--- Set the dimension --- */

  nDim = geometry_aux->GetnDim();

  /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

  geometry_aux->SetColorGrid_Parallel(config);

  /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

  geometry = NULL;
  geometry = new CGeometry *[config->GetnMGLevels()+1];

  /*--- Build the grid data structures using the ParMETIS coloring. ---*/

  geometry[MESH_0] = new CPhysicalGeometry(geometry_aux, config);

  /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

  delete geometry_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetBoundaries(config);

  fea = ((config->GetKind_Solver() == FEM_ELASTICITY) ||
         (config->GetKind_Solver() == DISC_ADJ_FEM));

  /*--- Compute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Setting point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();

  /*--- Renumbering points using Reverse Cuthill McKee ordering ---*/

  if (rank == MASTER_NODE) cout << "Renumbering points (Reverse Cuthill McKee Ordering)." << endl;
  geometry[MESH_0]->SetRCM_Ordering(config);

  /*--- recompute elements surrounding points, points surrounding points ---*/

  if (rank == MASTER_NODE) cout << "Recomputing point connectivity." << endl;
  geometry[MESH_0]->SetPoint_Connectivity();

  /*--- Compute elements surrounding elements ---*/

  if (rank == MASTER_NODE) cout << "Setting element connectivity." << endl;
  geometry[MESH_0]->SetElement_Connectivity();

  /*--- Check the orientation before computing geometrical quantities ---*/

  geometry[MESH_0]->SetBoundVolume();
  if (config->GetReorientElements()) {
    if (rank == MASTER_NODE) cout << "Checking the numerical grid orientation." << endl;
    geometry[MESH_0]->Check_IntElem_Orientation(config);
    geometry[MESH_0]->Check_BoundElem_Orientation(config);
  }

  /*--- Create the edge structure ---*/

  if (rank == MASTER_NODE) cout << "Identifying edges and vertices." << endl;
  geometry[MESH_0]->SetEdges();
  geometry[MESH_0]->SetVertex(config);

  /*--- Compute cell center of gravity ---*/

  if ((rank == MASTER_NODE) && (!fea)) cout << "Computing centers of gravity." << endl;
  geometry[MESH_0]->SetCoord_CG();

  /*--- Create the control volume structures ---*/

  if ((rank == MASTER_NODE) && (!fea)) cout << "Setting the control volume structure." << endl;
  geometry[MESH_0]->SetControlVolume(config, ALLOCATE);
  geometry[MESH_0]->SetBoundControlVolume(config, ALLOCATE);

  /*--- Visualize a dual control volume if requested ---*/

  if ((config->GetVisualize_CV() >= 0) &&
      (config->GetVisualize_CV() < (long)geometry[MESH_0]->GetnPointDomain()))
    geometry[MESH_0]->VisualizeControlVolume(config, UPDATE);

  /*--- Identify closest normal neighbor ---*/

  if (rank == MASTER_NODE) cout << "Searching for the closest normal neighbors to the surfaces." << endl;
  geometry[MESH_0]->FindNormal_Neighbor(config);

  /*--- Store the global to local mapping. ---*/

  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
  geometry[MESH_0]->SetGlobal_to_Local_Point();

  /*--- Compute the surface curvature ---*/

  if ((rank == MASTER_NODE) && (!fea)) cout << "Compute the surface curvature." << endl;
  geometry[MESH_0]->ComputeSurf_Curvature(config);

  /*--- Check for periodicity and disable MG if necessary. ---*/

  if (rank == MASTER_NODE) cout << "Checking for periodicity." << endl;
  geometry[MESH_0]->Check_Periodicity(config);

  /*--- Compute mesh quality statistics on the fine grid. ---*/

  if (!fea) {
    if (rank == MASTER_NODE)
      cout << "Computing mesh quality statistics for the dual control volumes." << endl;
    geometry[MESH_0]->ComputeMeshQualityStatistics(config);
  }

  geometry[MESH_0]->SetMGLevel(MESH_0);
  if ((config->GetnMGLevels() != 0) && (rank == MASTER_NODE))
    cout << "Setting the multigrid structure." << endl;

  /*--- Loop over all the new grid ---*/

  for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Create main agglomeration structure ---*/

    geometry[iMGlevel] = new CMultiGridGeometry(geometry, config, iMGlevel);

    /*--- Compute points surrounding points. ---*/

    geometry[iMGlevel]->SetPoint_Connectivity(geometry[iMGlevel-1]);

    /*--- Create the edge structure ---*/

    geometry[iMGlevel]->SetEdges();
    geometry[iMGlevel]->SetVertex(geometry[iMGlevel-1], config);

    /*--- Create the control volume structures ---*/

    geometry[iMGlevel]->SetControlVolume(config, geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetBoundControlVolume(config, geometry[iMGlevel-1], ALLOCATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGlevel-1]);

    /*--- Find closest neighbor to a surface point ---*/

    geometry[iMGlevel]->FindNormal_Neighbor(config);

    /*--- Store our multigrid index. ---*/

    geometry[iMGlevel]->SetMGLevel(iMGlevel);

    /*--- Protect against the situation that we were not able to complete
       the agglomeration for this level, i.e., there weren't enough points.
       We need to check if we changed the total number of levels and delete
       the incomplete CMultiGridGeometry object. ---*/

    if (config->GetnMGLevels() != requestedMGlevels) {
      delete geometry[iMGlevel];
      break;
    }

  }

  /*--- For unsteady simulations, initialize the grid volumes
   and coordinates for previous solutions. Loop over all zones/grids ---*/

  if (config->GetTime_Marching() && config->GetGrid_Movement()) {
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      for (iPoint = 0; iPoint < geometry[iMGlevel]->GetnPoint(); iPoint++) {

        /*--- Update cell volume ---*/

        geometry[iMGlevel]->node[iPoint]->SetVolume_n();
        geometry[iMGlevel]->node[iPoint]->SetVolume_nM1();

        /*--- Update point coordinates ---*/
        geometry[iMGlevel]->node[iPoint]->SetCoord_n();
        geometry[iMGlevel]->node[iPoint]->SetCoord_n1();

      }
    }
  }


  /*--- Create the data structure for MPI point-to-point communications. ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
    geometry[iMGlevel]->PreprocessP2PComms(geometry[iMGlevel], config);


  /*--- Perform a few preprocessing routines and communications. ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Compute the max length. ---*/

    if ((rank == MASTER_NODE) && (!fea) && (iMGlevel == MESH_0)) cout << "Finding max control volume width." << endl;
    geometry[iMGlevel]->SetMaxLength(config);

    /*--- Communicate the number of neighbors. This is needed for
         some centered schemes and for multigrid in parallel. ---*/

    if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (!fea) && (iMGlevel == MESH_0)) cout << "Communicating number of neighbors." << endl;
    geometry[iMGlevel]->InitiateComms(geometry[iMGlevel], config, NEIGHBORS);
    geometry[iMGlevel]->CompleteComms(geometry[iMGlevel], config, NEIGHBORS);
  }

}

void CDriver::Geometrical_Preprocessing_DGFEM(CConfig* config, CGeometry **&geometry) {

  /*--- Definition of the geometry class to store the primal grid in the
     partitioning process. ---*/

  CGeometry *geometry_aux = NULL;

  /*--- All ranks process the grid and call ParMETIS for partitioning ---*/

  geometry_aux = new CPhysicalGeometry(config, iZone, nZone);

  /*--- Set the dimension --- */

  nDim = geometry_aux->GetnDim();

  /*--- For the FEM solver with time-accurate local time-stepping, use
       a dummy solver class to retrieve the initial flow state. ---*/

  CSolver *solver_aux = NULL;
  solver_aux = new CFEM_DG_EulerSolver(config, nDim, MESH_0);

  /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/

  geometry_aux->SetColorFEMGrid_Parallel(config);

  /*--- Allocate the memory of the current domain, and divide the grid
     between the ranks. ---*/

  geometry = NULL;
  geometry = new CGeometry *[config->GetnMGLevels()+1];

  geometry[MESH_0] = new CMeshFEM_DG(geometry_aux, config);

  /*--- Deallocate the memory of geometry_aux and solver_aux ---*/

  delete geometry_aux;
  if (solver_aux != NULL) delete solver_aux;

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetSendReceive(config);

  /*--- Add the Send/Receive boundaries ---*/
  geometry[MESH_0]->SetBoundaries(config);

  /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
       define all virtual functions in the base class CGeometry. ---*/
  CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry[MESH_0]);

  /*--- Determine the standard elements for the volume elements. ---*/
  if (rank == MASTER_NODE) cout << "Creating standard volume elements." << endl;
  DGMesh->CreateStandardVolumeElements(config);

  /*--- Create the face information needed to compute the contour integral
       for the elements in the Discontinuous Galerkin formulation. ---*/
  if (rank == MASTER_NODE) cout << "Creating face information." << endl;
  DGMesh->CreateFaces(config);

  /*--- Compute the metric terms of the volume elements. ---*/
  if (rank == MASTER_NODE) cout << "Computing metric terms volume elements." << endl;
  DGMesh->MetricTermsVolumeElements(config);

  /*--- Compute the metric terms of the surface elements. ---*/
  if (rank == MASTER_NODE) cout << "Computing metric terms surface elements." << endl;
  DGMesh->MetricTermsSurfaceElements(config);

  /*--- Compute a length scale of the volume elements. ---*/
  if (rank == MASTER_NODE) cout << "Computing length scale volume elements." << endl;
  DGMesh->LengthScaleVolumeElements();

  /*--- Compute the coordinates of the integration points. ---*/
  if (rank == MASTER_NODE) cout << "Computing coordinates of the integration points." << endl;
  DGMesh->CoordinatesIntegrationPoints();

  /*--- Compute the coordinates of the location of the solution DOFs. This is different
            from the grid points when a different polynomial degree is used to represent the
            geometry and solution. ---*/
  if (rank == MASTER_NODE) cout << "Computing coordinates of the solution DOFs." << endl;
  DGMesh->CoordinatesSolDOFs();

  /*--- Perform the preprocessing tasks when wall functions are used. ---*/
  if (rank == MASTER_NODE) cout << "Preprocessing for the wall functions. " << endl;
  DGMesh->WallFunctionPreprocessing(config);

  /*--- Store the global to local mapping. ---*/
  if (rank == MASTER_NODE) cout << "Storing a mapping from global to local DOF index." << endl;
  geometry[MESH_0]->SetGlobal_to_Local_Point();


  /*--- Loop to create the coarser grid levels. ---*/

  for(unsigned short iMGlevel=1; iMGlevel<=config->GetnMGLevels(); iMGlevel++) {

    SU2_MPI::Error("Geometrical_Preprocessing_DGFEM: Coarse grid levels not implemented yet.",
                   CURRENT_FUNCTION);
  }

}

void CDriver::Solver_Preprocessing(CConfig* config, CGeometry** geometry, CSolver ***&solver) {

  unsigned short iSol;

  if (rank == MASTER_NODE)
    cout << endl <<"-------------------- Solver Preprocessing ( Zone " << config->GetiZone() <<" ) --------------------" << endl;

  solver = new CSolver**[config->GetnMGLevels()+1];

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    solver[iMesh] = NULL;

  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    solver[iMesh] = new CSolver* [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      solver[iMesh][iSol] = NULL;
  }

  unsigned short iMGlevel;
  bool euler, ns, turbulent,
      fem_euler, fem_ns, fem_turbulent, fem_transition,
      adj_euler, adj_ns, adj_turb,
      heat, disc_adj_heat,
      fem, disc_adj_fem,
      spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
      template_solver, disc_adj, disc_adj_turb,
      fem_dg_flow, fem_dg_shock_persson,
      e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  /*--- Count the number of DOFs per solution point. ---*/

  DOFsPerPoint = 0;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent     = false;
  fem_euler        = false;  fem_ns          = false;  fem_turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb      = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat             = false;  disc_adj_heat    = false;
  transition       = false;  fem_transition   = false;
  template_solver  = false;
  fem_dg_flow      = false;  fem_dg_shock_persson = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  bool compressible   = false;
  bool incompressible = false;

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : euler = true; compressible = true; break;
    case NAVIER_STOKES: ns = true; compressible = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case RANS : ns = true; turbulent = true; compressible = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case INC_EULER : euler = true; incompressible = true; break;
    case INC_NAVIER_STOKES: ns = true; incompressible = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case INC_RANS : ns = true; turbulent = true; incompressible = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; compressible = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; compressible = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; compressible = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case FEM_LES : fem_ns = true; compressible = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; compressible = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); compressible = true; adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; compressible = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: euler = true; disc_adj = true; compressible = true; break;
    case DISC_ADJ_NAVIER_STOKES: ns = true; disc_adj = true; compressible = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: ns = true; turbulent = true; disc_adj = true; compressible = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; incompressible = true; break;
    case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; incompressible = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; incompressible = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; compressible = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; compressible = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; compressible = true; if(config->GetKind_Trans_Model() == LM) fem_transition = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; compressible = true; break;
    case DISC_ADJ_HEAT: heat = true; disc_adj_heat = true; break;
  }

  /*--- Determine the kind of FEM solver used for the flow. ---*/

  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }

  /*--- Determine the kind of shock capturing method for FEM DG solver. ---*/

  switch( config->GetKind_FEM_DG_Shock() ) {
    case PERSSON: fem_dg_shock_persson = true; break;
  }

  /*--- Assign turbulence model booleans ---*/

  if (turbulent || fem_turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:        spalart_allmaras = true;        break;
      case SA_NEG:    neg_spalart_allmaras = true;    break;
      case SA_E:      e_spalart_allmaras = true;      break;
      case SA_COMP:   comp_spalart_allmaras = true;   break;
      case SA_E_COMP: e_comp_spalart_allmaras = true; break;
      case SST:       menter_sst = true;              break;
      case SST_SUST:  menter_sst = true;              break;
      default: SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION); break;
    }

  /*--- Definition of the Class for the solution: solver[DOMAIN][INSTANCE][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- Allocate solution for a template problem ---*/

    if (template_solver) {
      solver[iMGlevel][TEMPLATE_SOL] = new CTemplateSolver(geometry[iMGlevel], config);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][TEMPLATE_SOL]->GetnVar();
    }

    /*--- Allocate solution for direct problem, and run the preprocessing and postprocessing ---*/

    if (euler) {
      if (compressible) {
        solver[iMGlevel][FLOW_SOL] = new CEulerSolver(geometry[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (incompressible) {
        solver[iMGlevel][FLOW_SOL] = new CIncEulerSolver(geometry[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][FLOW_SOL]->GetnVar();
    }
    if (ns) {
      if (compressible) {
        solver[iMGlevel][FLOW_SOL] = new CNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        solver[iMGlevel][FLOW_SOL] = new CIncNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][FLOW_SOL]->GetnVar();
    }
    if (turbulent) {
      if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras || neg_spalart_allmaras) {
        solver[iMGlevel][TURB_SOL] = new CTurbSASolver(geometry[iMGlevel], config, iMGlevel, solver[iMGlevel][FLOW_SOL]->GetFluidModel() );
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel);
      }
      else if (menter_sst) {
        solver[iMGlevel][TURB_SOL] = new CTurbSSTSolver(geometry[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        solver[iMGlevel][TURB_SOL]->Postprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel);
        solver[iMGlevel][FLOW_SOL]->Preprocessing(geometry[iMGlevel], solver[iMGlevel], config, iMGlevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][TURB_SOL]->GetnVar();
      if (transition) {
        solver[iMGlevel][TRANS_SOL] = new CTransLMSolver(geometry[iMGlevel], config, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][TRANS_SOL]->GetnVar();
      }
    }
    if (fem_euler) {
      if( fem_dg_flow ) {
        if( fem_dg_shock_persson ) {
          solver[iMGlevel][FLOW_SOL] = new CFEM_DG_NSSolver(geometry[iMGlevel], config, iMGlevel);
        }
        else {
          solver[iMGlevel][FLOW_SOL] = new CFEM_DG_EulerSolver(geometry[iMGlevel], config, iMGlevel);
        }
      }
    }
    if (fem_ns) {
      if( fem_dg_flow )
        solver[iMGlevel][FLOW_SOL] = new CFEM_DG_NSSolver(geometry[iMGlevel], config, iMGlevel);
    }
    if (fem_turbulent) {
      SU2_MPI::Error("Finite element turbulence model not yet implemented.", CURRENT_FUNCTION);

      if(fem_transition)
        SU2_MPI::Error("Finite element transition model not yet implemented.", CURRENT_FUNCTION);
    }
    if (heat) {
      solver[iMGlevel][HEAT_SOL] = new CHeatSolver(geometry[iMGlevel], config, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][HEAT_SOL]->GetnVar();
    }
    if (fem) {
      solver[iMGlevel][FEA_SOL] = new CFEASolver(geometry[iMGlevel], config);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][FEA_SOL]->GetnVar();
    }
    if (config->AddRadiation()) {
      solver[iMGlevel][RAD_SOL] = new CRadP1Solver(geometry[iMGlevel], config);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][RAD_SOL]->GetnVar();
    }

    /*--- Allocate solution for adjoint problem ---*/

    if (adj_euler) {
      if (compressible) {
        solver[iMGlevel][ADJFLOW_SOL] = new CAdjEulerSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        SU2_MPI::Error("Continuous adjoint for the incompressible solver is not currently available.", CURRENT_FUNCTION);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFLOW_SOL]->GetnVar();
    }
    if (adj_ns) {
      if (compressible) {
        solver[iMGlevel][ADJFLOW_SOL] = new CAdjNSSolver(geometry[iMGlevel], config, iMGlevel);
      }
      if (incompressible) {
        SU2_MPI::Error("Continuous adjoint for the incompressible solver is not currently available.", CURRENT_FUNCTION);
      }
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFLOW_SOL]->GetnVar();
    }
    if (adj_turb) {
      solver[iMGlevel][ADJTURB_SOL] = new CAdjTurbSolver(geometry[iMGlevel], config, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJTURB_SOL]->GetnVar();
    }

    if (disc_adj) {
      solver[iMGlevel][ADJFLOW_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][FLOW_SOL], RUNTIME_FLOW_SYS, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFLOW_SOL]->GetnVar();
      if (disc_adj_turb) {
        solver[iMGlevel][ADJTURB_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][TURB_SOL], RUNTIME_TURB_SYS, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJTURB_SOL]->GetnVar();
      }
      if (heat) {
        solver[iMGlevel][ADJHEAT_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][HEAT_SOL], RUNTIME_HEAT_SYS, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJHEAT_SOL]->GetnVar();
      }
      if (config->AddRadiation()){
        solver[iMGlevel][ADJRAD_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][RAD_SOL], RUNTIME_RADIATION_SYS, iMGlevel);
        if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJRAD_SOL]->GetnVar();
      }
    }

    if (disc_adj_fem) {
      solver[iMGlevel][ADJFEA_SOL] = new CDiscAdjFEASolver(geometry[iMGlevel], config, solver[iMGlevel][FEA_SOL], RUNTIME_FEA_SYS, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJFEA_SOL]->GetnVar();
    }

    if (disc_adj_heat) {
      solver[iMGlevel][ADJHEAT_SOL] = new CDiscAdjSolver(geometry[iMGlevel], config, solver[iMGlevel][HEAT_SOL], RUNTIME_HEAT_SYS, iMGlevel);
      if (iMGlevel == MESH_0) DOFsPerPoint += solver[iMGlevel][ADJHEAT_SOL]->GetnVar();
    }
  }

  /*--- Preprocess the mesh solver for dynamic meshes. ---*/
  /*--- This needs to be done before solver restart so the old coordinates are stored. ---*/
  if (config->GetDeform_Mesh()){
    solver[MESH_0][MESH_SOL] = new CMeshSolver(geometry[MESH_0], config);

    if (config->GetDiscrete_Adjoint())
      solver[MESH_0][ADJMESH_SOL] = new CDiscAdjMeshSolver(geometry[MESH_0], config, solver[MESH_0][MESH_SOL]);

  }

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool update_geo = true;
  if (config->GetFSI_Simulation()) update_geo = false;

  Solver_Restart(solver, geometry, config, update_geo);

  /*--- Set up any necessary inlet profiles ---*/

  Inlet_Preprocessing(solver, geometry, config);

}

void CDriver::Inlet_Preprocessing(CSolver ***solver, CGeometry **geometry,
                                  CConfig *config) {

  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat,
  fem,
  template_solver, disc_adj, disc_adj_fem, disc_adj_turb;
  int val_iter = 0;
  unsigned short iMesh;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  heat             = false;  disc_adj_turb    = false;
  template_solver  = false;

  /*--- Adjust iteration number for unsteady restarts. ---*/

  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());

  if (dual_time) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else if (config->GetTime_Marching() == DT_STEPPING_1ST)
      val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    else val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
  }

  if (time_stepping) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
  }

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case INC_EULER: euler = true; break;
    case NAVIER_STOKES: case INC_NAVIER_STOKES: ns = true; break;
    case RANS : case INC_RANS: ns = true; turbulent = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; break;
    case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
  }


  /*--- Load inlet profile files for any of the active solver containers.
   Note that these routines fill the fine grid data structures for the markers
   and restrict values down to all coarser MG levels. ---*/

  if (config->GetInlet_Profile_From_File()) {

    /*--- Use LoadInletProfile() routines for the particular solver. ---*/

    if (rank == MASTER_NODE) {
      cout << endl;
      cout << "Reading inlet profile from file: ";
      cout << config->GetInlet_FileName() << endl;
    }

    bool no_profile = false;

    if (euler || ns || adj_euler || adj_ns || disc_adj) {
      solver[MESH_0][FLOW_SOL]->LoadInletProfile(geometry, solver, config, val_iter, FLOW_SOL, INLET_FLOW);
    }
    if (turbulent || adj_turb || disc_adj_turb) {
      solver[MESH_0][TURB_SOL]->LoadInletProfile(geometry, solver, config, val_iter, TURB_SOL, INLET_FLOW);
    }

    if (template_solver) {
      no_profile = true;
    }
    if (heat) {
      no_profile = true;
    }
    if (fem) {
      no_profile = true;
    }
    if (disc_adj_fem) {
      no_profile = true;
    }

    /*--- Exit if profiles were requested for a solver that is not available. ---*/

    if (no_profile) {
      SU2_MPI::Error(string("Inlet profile specification via file (C++) has not been \n") +
                     string("implemented yet for this solver.\n") +
                     string("Please set SPECIFIED_INLET_PROFILE= NO and try again."), CURRENT_FUNCTION);
    }

  } else {

    /*--- Uniform inlets or python-customized inlets ---*/

    /* --- Initialize quantities for inlet boundary
     * This routine does not check if they python wrapper is being used to
     * set custom boundary conditions.  This is intentional; the
     * default values for python custom BCs are initialized with the default
     * values specified in the config (avoiding non physical values) --- */

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (euler || ns || adj_euler || adj_ns || disc_adj)
          solver[iMesh][FLOW_SOL]->SetUniformInlet(config, iMarker);
        if (turbulent)
          solver[iMesh][TURB_SOL]->SetUniformInlet(config, iMarker);
      }
    }

  }

}

void CDriver::Solver_Restart(CSolver ***solver, CGeometry **geometry,
                             CConfig *config, bool update_geo) {

  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat, fem, fem_euler, fem_ns, fem_dg_flow,
  template_solver, disc_adj, disc_adj_fem, disc_adj_turb, disc_adj_heat;
  int val_iter = 0;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns           = false;  turbulent   = false;
  adj_euler        = false;  adj_ns       = false;  adj_turb    = false;
  fem_euler        = false;  fem_ns       = false;  fem_dg_flow = false;
  disc_adj         = false;
  fem              = false;  disc_adj_fem     = false;
  disc_adj_turb    = false;
  heat             = false;  disc_adj_heat    = false;
  template_solver  = false;

  /*--- Check for restarts and use the LoadRestart() routines. ---*/

  bool restart      = config->GetRestart();
  bool restart_flow = config->GetRestart_Flow();
  bool no_restart   = false;

  /*--- Adjust iteration number for unsteady restarts. ---*/

  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;
  bool adjoint = (config->GetDiscrete_Adjoint() || config->GetContinuous_Adjoint());
  bool time_domain = (config->GetTime_Domain()); // Dynamic simulation (FSI).

  if (dual_time) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else if (config->GetTime_Marching() == DT_STEPPING_1ST)
      val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    else val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
  }

  if (time_stepping) {
    if (adjoint) val_iter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
    else val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
  }

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case INC_EULER: euler = true; break;
    case NAVIER_STOKES: case INC_NAVIER_STOKES: ns = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case RANS : case INC_RANS: ns = true; turbulent = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; break;
    case FEM_LES : fem_ns = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat = true; disc_adj_heat = true; break;

  }

  /*--- Determine the kind of FEM solver used for the flow. ---*/

  switch( config->GetKind_FEM_Flow() ) {
    case DG: fem_dg_flow = true; break;
  }

  /*--- Load restarts for any of the active solver containers. Note that
   these restart routines fill the fine grid and interpolate to all MG levels. ---*/

  if (restart || restart_flow) {
    if (euler || ns) {
      solver[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (turbulent) {
      solver[MESH_0][TURB_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (config->AddRadiation()) {
      solver[MESH_0][RAD_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (fem) {
      if (time_domain) val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      solver[MESH_0][FEA_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (fem_euler || fem_ns) {
      if (fem_dg_flow)
        solver[MESH_0][FLOW_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (heat) {
      solver[MESH_0][HEAT_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
  }

  if (restart) {
    if (template_solver) {
      no_restart = true;
    }
    if (heat) {
      solver[MESH_0][HEAT_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (adj_euler || adj_ns) {
      solver[MESH_0][ADJFLOW_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (adj_turb) {
      no_restart = true;
    }
    if (disc_adj) {
      solver[MESH_0][ADJFLOW_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
      if (disc_adj_turb)
        solver[MESH_0][ADJTURB_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
      if (disc_adj_heat)
        solver[MESH_0][ADJHEAT_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
      if (config->AddRadiation())
        solver[MESH_0][ADJRAD_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (disc_adj_fem) {
        if (time_domain) val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
        solver[MESH_0][ADJFEA_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
    if (disc_adj_heat) {
      solver[MESH_0][ADJHEAT_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
    }
  }

  if ((restart || restart_flow) && config->GetDeform_Mesh() && update_geo){
    /*--- Always restart with the last state ---*/
    val_iter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
    solver[MESH_0][MESH_SOL]->LoadRestart(geometry, solver, config, val_iter, update_geo);
  }

  /*--- Exit if a restart was requested for a solver that is not available. ---*/

  if (no_restart) {
    SU2_MPI::Error(string("A restart capability has not been implemented yet for this solver.\n") +
                   string("Please set RESTART_SOL= NO and try again."), CURRENT_FUNCTION);
  }

  /*--- Think about calls to pre / post-processing here, plus realizability checks. ---*/


}

void CDriver::Solver_Postprocessing(CSolver ****solver, CGeometry **geometry,
                                    CConfig *config, unsigned short val_iInst) {
  unsigned short iMGlevel;
  bool euler, ns, turbulent,
  adj_euler, adj_ns, adj_turb,
  heat, fem,
  spalart_allmaras, neg_spalart_allmaras, menter_sst, transition,
  template_solver, disc_adj, disc_adj_turb, disc_adj_fem, disc_adj_heat,
  e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras;

  /*--- Initialize some useful booleans ---*/

  euler            = false;  ns              = false;  turbulent = false;
  adj_euler        = false;  adj_ns          = false;  adj_turb  = false;
  spalart_allmaras = false;  menter_sst      = false;  disc_adj_turb = false;
  neg_spalart_allmaras = false;
  disc_adj        = false;
  fem              = false;  disc_adj_fem    = false;
  heat             = false;   disc_adj_heat   = false;
  transition       = false;
  template_solver  = false;
  e_spalart_allmaras = false; comp_spalart_allmaras = false; e_comp_spalart_allmaras = false;

  /*--- Assign booleans ---*/

  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case INC_EULER: euler = true; break;
    case NAVIER_STOKES: case INC_NAVIER_STOKES: ns = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case RANS : case INC_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : euler = true; break;
    case FEM_NAVIER_STOKES:
    case FEM_LES: ns = true; break;
    case FEM_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER: case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: ns = true; turbulent = true; disc_adj = true; disc_adj_turb = (!config->GetFrozen_Visc_Disc()); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat = true; disc_adj_heat = true; break;
  }

  /*--- Assign turbulence model booleans ---*/

  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
    case SA:        spalart_allmaras = true;        break;
    case SA_NEG:    neg_spalart_allmaras = true;    break;
    case SA_E:      e_spalart_allmaras = true;      break;
    case SA_COMP:   comp_spalart_allmaras = true;   break;
    case SA_E_COMP: e_comp_spalart_allmaras = true; break;
    case SST:       menter_sst = true;              break;
    case SST_SUST:  menter_sst = true;              break;
    default: SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION); break;
    }

  /*--- Definition of the Class for the solution: solver_container[DOMAIN][MESH_LEVEL][EQUATION]. Note that euler, ns
   and potential are incompatible, they use the same position in sol container ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    /*--- DeAllocate solution for a template problem ---*/

    if (template_solver) {
      delete solver[val_iInst][iMGlevel][TEMPLATE_SOL];
    }

    /*--- DeAllocate solution for adjoint problem ---*/

    if (adj_euler || adj_ns || disc_adj) {
      delete solver[val_iInst][iMGlevel][ADJFLOW_SOL];
      if (disc_adj_turb || adj_turb) {
        delete solver[val_iInst][iMGlevel][ADJTURB_SOL];
      }
      if (heat) {
        delete solver[val_iInst][iMGlevel][ADJHEAT_SOL];
      }
    }

    if (disc_adj_heat) {
      delete solver[val_iInst][iMGlevel][ADJHEAT_SOL];
    }

    /*--- DeAllocate solution for direct problem ---*/

    if (euler || ns) {
      delete solver[val_iInst][iMGlevel][FLOW_SOL];
    }

    if (turbulent) {
      if (spalart_allmaras || neg_spalart_allmaras || menter_sst || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
        delete solver[val_iInst][iMGlevel][TURB_SOL];
      }
      if (transition) {
        delete solver[val_iInst][iMGlevel][TRANS_SOL];
      }
    }
    if (heat) {
      delete solver[val_iInst][iMGlevel][HEAT_SOL];
    }
    if (fem) {
      delete solver[val_iInst][iMGlevel][FEA_SOL];
    }
    if (disc_adj_fem) {
      delete solver[val_iInst][iMGlevel][ADJFEA_SOL];
    }
    if (config->AddRadiation()) {
      delete solver[val_iInst][iMGlevel][RAD_SOL];
      if (disc_adj) delete solver[val_iInst][iMGlevel][ADJRAD_SOL];
    }

    if (iMGlevel == 0){
      if (config->GetDeform_Mesh()){
        delete solver[val_iInst][MESH_0][MESH_SOL];
        if (config->GetDiscrete_Adjoint())
          delete solver[val_iInst][MESH_0][ADJMESH_SOL];
      }
    }

    delete [] solver[val_iInst][iMGlevel];

  }

  delete [] solver[val_iInst];

}

void CDriver::Integration_Preprocessing(CConfig *config, CIntegration **&integration) {

  unsigned short iSol;

  if (rank == MASTER_NODE)
    cout << endl <<"----------------- Integration Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  integration = new CIntegration* [MAX_SOLS];
  for (iSol = 0; iSol < MAX_SOLS; iSol++)
    integration[iSol] = NULL;

  bool euler, adj_euler, ns, adj_ns, turbulent, adj_turb, fem,
      fem_euler, fem_ns, fem_turbulent,
      heat, template_solver, transition, disc_adj, disc_adj_fem, disc_adj_heat;

  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  disc_adj         = false;
  fem_euler        = false;
  fem_ns           = false;
  fem_turbulent    = false;
  heat             = false; disc_adj_heat    = false;
  fem              = false; disc_adj_fem     = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case INC_EULER: euler = true; break;
    case NAVIER_STOKES: case INC_NAVIER_STOKES: ns = true;  heat = config->GetWeakly_Coupled_Heat(); break;
    case RANS : case INC_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; break;
    case FEM_LES :  fem_ns = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER : case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS : case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat = true; disc_adj_heat = true; break;
  }

  /*--- Allocate solution for a template problem ---*/
  if (template_solver) integration[TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- Allocate solution for direct problem ---*/
  if (euler) integration[FLOW_SOL] = new CMultiGridIntegration(config);
  if (ns) integration[FLOW_SOL] = new CMultiGridIntegration(config);
  if (turbulent) integration[TURB_SOL] = new CSingleGridIntegration(config);
  if (transition) integration[TRANS_SOL] = new CSingleGridIntegration(config);
  if (heat) integration[HEAT_SOL] = new CSingleGridIntegration(config);
  if (fem) integration[FEA_SOL] = new CStructuralIntegration(config);
  if (config->AddRadiation()) integration[RAD_SOL] = new CSingleGridIntegration(config);

  /*--- Allocate integration container for finite element flow solver. ---*/

  if (fem_euler) integration[FLOW_SOL] = new CFEM_DG_Integration(config);
  if (fem_ns)    integration[FLOW_SOL] = new CFEM_DG_Integration(config);
  //if (fem_turbulent) integration[FEM_TURB_SOL] = new CSingleGridIntegration(config);

  if (fem_turbulent)
    SU2_MPI::Error("No turbulent FEM solver yet", CURRENT_FUNCTION);

  /*--- Allocate solution for adjoint problem ---*/
  if (adj_euler) integration[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_ns) integration[ADJFLOW_SOL] = new CMultiGridIntegration(config);
  if (adj_turb) integration[ADJTURB_SOL] = new CSingleGridIntegration(config);

  if (disc_adj) integration[ADJFLOW_SOL] = new CIntegration(config);
  if (disc_adj_fem) integration[ADJFEA_SOL] = new CIntegration(config);
  if (disc_adj_heat) integration[ADJHEAT_SOL] = new CIntegration(config);

}

void CDriver::Integration_Postprocessing(CIntegration ***integration, CGeometry **geometry, CConfig *config, unsigned short val_iInst) {
  bool euler, adj_euler, ns, adj_ns, turbulent, adj_turb, fem,
      fem_euler, fem_ns, fem_turbulent,
      heat, template_solver, transition, disc_adj, disc_adj_fem, disc_adj_heat;

  /*--- Initialize some useful booleans ---*/
  euler            = false; adj_euler        = false;
  ns               = false; adj_ns           = false;
  turbulent        = false; adj_turb         = false;
  disc_adj         = false;
  fem_euler        = false;
  fem_ns           = false;
  fem_turbulent    = false;
  heat             = false; disc_adj_heat    = false;
  fem              = false; disc_adj_fem     = false;
  transition       = false;
  template_solver  = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER: template_solver = true; break;
    case EULER : case INC_EULER: euler = true; break;
    case NAVIER_STOKES: case INC_NAVIER_STOKES: ns = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case RANS : case INC_RANS: ns = true; turbulent = true; if (config->GetKind_Trans_Model() == LM) transition = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case FEM_EULER : fem_euler = true; break;
    case FEM_NAVIER_STOKES: fem_ns = true; break;
    case FEM_RANS : fem_ns = true; fem_turbulent = true; break;
    case FEM_LES :  fem_ns = true; break;
    case HEAT_EQUATION: heat = true; break;
    case FEM_ELASTICITY: fem = true; break;
    case ADJ_EULER : euler = true; adj_euler = true; break;
    case ADJ_NAVIER_STOKES : ns = true; turbulent = (config->GetKind_Turb_Model() != NONE); adj_ns = true; break;
    case ADJ_RANS : ns = true; turbulent = true; adj_ns = true; adj_turb = (!config->GetFrozen_Visc_Cont()); break;
    case DISC_ADJ_EULER : case DISC_ADJ_INC_EULER: euler = true; disc_adj = true; break;
    case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_NAVIER_STOKES: ns = true; disc_adj = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_RANS : case DISC_ADJ_INC_RANS: ns = true; turbulent = true; disc_adj = true; heat = config->GetWeakly_Coupled_Heat(); break;
    case DISC_ADJ_FEM_EULER: fem_euler = true; disc_adj = true; break;
    case DISC_ADJ_FEM_NS: fem_ns = true; disc_adj = true; break;
    case DISC_ADJ_FEM_RANS: fem_ns = true; fem_turbulent = true; disc_adj = true; break;
    case DISC_ADJ_FEM: fem = true; disc_adj_fem = true; break;
    case DISC_ADJ_HEAT: heat = true; disc_adj_heat = true; break;
  }

  /*--- DeAllocate solution for a template problem ---*/
  if (template_solver) integration[val_iInst][TEMPLATE_SOL] = new CSingleGridIntegration(config);

  /*--- DeAllocate solution for direct problem ---*/
  if (euler || ns) delete integration[val_iInst][FLOW_SOL];
  if (turbulent) delete integration[val_iInst][TURB_SOL];
  if (transition) delete integration[val_iInst][TRANS_SOL];
  if (heat) delete integration[val_iInst][HEAT_SOL];
  if (fem) delete integration[val_iInst][FEA_SOL];
  if (disc_adj_fem) delete integration[val_iInst][ADJFEA_SOL];
  if (disc_adj_heat) delete integration[val_iInst][ADJHEAT_SOL];
  if (config->AddRadiation()) delete integration[val_iInst][RAD_SOL];

  /*--- DeAllocate solution for adjoint problem ---*/
  if (adj_euler || adj_ns || disc_adj) delete integration[val_iInst][ADJFLOW_SOL];
  if (adj_turb) delete integration[val_iInst][ADJTURB_SOL];

  /*--- DeAllocate integration container for finite element flow solver. ---*/
  if (fem_euler || fem_ns) delete integration[val_iInst][FLOW_SOL];
  //if (fem_turbulent)     delete integration_container[val_iInst][FEM_TURB_SOL];

  if (fem_turbulent)
    SU2_MPI::Error("No turbulent FEM solver yet", CURRENT_FUNCTION);

  delete [] integration[val_iInst];
}

void CDriver::Numerics_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CNumerics ****&numerics) {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Numerics Preprocessing ( Zone " << config->GetiZone() <<" ) -------------------" << endl;

  unsigned short iMGlevel, iSol,

  nVar_Template         = 0,
  nVar_Flow             = 0,
  nVar_Trans            = 0,
  nVar_Turb             = 0,
  nVar_Adj_Flow         = 0,
  nVar_Adj_Turb         = 0,
  nVar_FEM              = 0,
  nVar_Rad              = 0,
  nVar_Heat             = 0;

  numerics = new CNumerics***[config->GetnMGLevels()+1];

  const su2double *constants = nullptr;
  su2double kine_Inf = 0.0, omega_Inf = 0.0;

  bool compressible = false;
  bool incompressible = false;
  bool ideal_gas = (config->GetKind_FluidModel() == STANDARD_AIR) || (config->GetKind_FluidModel() == IDEAL_GAS);
  bool roe_low_dissipation = (config->GetKind_RoeLowDiss() != NO_ROELOWDISS);

  /*--- Initialize some useful booleans ---*/
  bool euler, ns, turbulent, adj_euler, adj_ns, adj_turb, fem_euler, fem_ns, fem_turbulent;
  bool spalart_allmaras, neg_spalart_allmaras, e_spalart_allmaras, comp_spalart_allmaras, e_comp_spalart_allmaras, menter_sst;
  bool fem, heat, transition, template_solver;

  euler = ns = turbulent = adj_euler = adj_ns = adj_turb = fem_euler = fem_ns = fem_turbulent = false;
  spalart_allmaras = neg_spalart_allmaras = e_spalart_allmaras = comp_spalart_allmaras = e_comp_spalart_allmaras = menter_sst = false;
  fem = heat = transition = template_solver = false;

  /*--- Assign booleans ---*/
  switch (config->GetKind_Solver()) {
    case TEMPLATE_SOLVER:
      template_solver = true; break;

    case EULER :
    case DISC_ADJ_EULER:
      euler = compressible = true; break;

    case NAVIER_STOKES:
    case DISC_ADJ_NAVIER_STOKES:
      ns = compressible = true; break;

    case RANS:
    case DISC_ADJ_RANS:
      ns = compressible = turbulent = true;
      transition = (config->GetKind_Trans_Model() == LM); break;

    case INC_EULER:
    case DISC_ADJ_INC_EULER:
      euler = incompressible = true; break;

    case INC_NAVIER_STOKES:
    case DISC_ADJ_INC_NAVIER_STOKES:
      ns = incompressible = true;
      heat = config->GetWeakly_Coupled_Heat(); break;

    case INC_RANS:
    case DISC_ADJ_INC_RANS:
      ns = incompressible = turbulent = true;
      heat = config->GetWeakly_Coupled_Heat();
      transition = (config->GetKind_Trans_Model() == LM); break;

    case FEM_EULER:
    case DISC_ADJ_FEM_EULER:
      fem_euler = compressible = true; break;

    case FEM_NAVIER_STOKES:
    case DISC_ADJ_FEM_NS:
      fem_ns = compressible = true; break;

    case FEM_RANS:
    case DISC_ADJ_FEM_RANS:
      fem_ns = compressible = fem_turbulent = true; break;

    case FEM_LES:
      fem_ns = compressible = true; break;

    case HEAT_EQUATION:
    case DISC_ADJ_HEAT:
      heat = true; break;

    case FEM_ELASTICITY:
    case DISC_ADJ_FEM:
      fem = true; break;

    case ADJ_EULER:
      adj_euler = euler = compressible = true; break;

    case ADJ_NAVIER_STOKES:
      adj_ns = ns = compressible = true;
      turbulent = (config->GetKind_Turb_Model() != NONE); break;

    case ADJ_RANS:
      adj_ns = ns = compressible = turbulent = true;
      adj_turb = !config->GetFrozen_Visc_Cont(); break;

  }

  /*--- Assign turbulence model booleans ---*/

  if (turbulent || fem_turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA:        spalart_allmaras = true;        break;
      case SA_NEG:    neg_spalart_allmaras = true;    break;
      case SA_E:      e_spalart_allmaras = true;      break;
      case SA_COMP:   comp_spalart_allmaras = true;   break;
      case SA_E_COMP: e_comp_spalart_allmaras = true; break;
      case SST:       menter_sst = true;              break;
      case SST_SUST:  menter_sst = true;              break;
      default:
        SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION);
        break;
    }

  /*--- If the Menter SST model is used, store the constants of the model and determine the
        free stream values of the turbulent kinetic energy and dissipation rate. ---*/

  if (menter_sst) {
    constants = solver[MESH_0][TURB_SOL]->GetConstants();
    kine_Inf  = solver[MESH_0][TURB_SOL]->GetTke_Inf();
    omega_Inf = solver[MESH_0][TURB_SOL]->GetOmega_Inf();
  }

  /*--- Number of variables for the template ---*/

  if (template_solver) nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();

  /*--- Number of variables for direct problem ---*/

  if (euler)        nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (ns)           nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (turbulent)    nVar_Turb = solver[MESH_0][TURB_SOL]->GetnVar();
  if (transition)   nVar_Trans = solver[MESH_0][TRANS_SOL]->GetnVar();

  if (fem_euler)    nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  if (fem_ns)       nVar_Flow = solver[MESH_0][FLOW_SOL]->GetnVar();
  //if (fem_turbulent)    nVar_Turb = solver_container[MESH_0][FEM_TURB_SOL]->GetnVar();

  if (fem)          nVar_FEM = solver[MESH_0][FEA_SOL]->GetnVar();
  if (heat)     nVar_Heat = solver[MESH_0][HEAT_SOL]->GetnVar();

  if (config->AddRadiation())    nVar_Rad = solver[MESH_0][RAD_SOL]->GetnVar();

  /*--- Number of variables for adjoint problem ---*/

  if (adj_euler)    nVar_Adj_Flow = solver[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_ns)       nVar_Adj_Flow = solver[MESH_0][ADJFLOW_SOL]->GetnVar();
  if (adj_turb)     nVar_Adj_Turb = solver[MESH_0][ADJTURB_SOL]->GetnVar();

  /*--- Definition of the Class for the numerical method: numerics_container[INSTANCE_LEVEL][MESH_LEVEL][EQUATION][EQ_TERM] ---*/

  for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
    numerics[iMGlevel] = new CNumerics** [MAX_SOLS];
    for (iSol = 0; iSol < MAX_SOLS; iSol++)
      numerics[iMGlevel][iSol] = new CNumerics* [MAX_TERMS*omp_get_max_threads()]();
  }

  /*--- Instantiate one numerics object per thread for each required term. ---*/

  for (int thread = 0; thread < omp_get_max_threads(); ++thread)
  {
  const int offset = thread * MAX_TERMS;

  const int conv_term = CONV_TERM + offset;
  const int visc_term = VISC_TERM + offset;

  const int source_first_term = SOURCE_FIRST_TERM + offset;
  const int source_second_term = SOURCE_SECOND_TERM + offset;

  const int conv_bound_term = CONV_BOUND_TERM + offset;
  const int visc_bound_term = VISC_BOUND_TERM + offset;

  const int fea_term = FEA_TERM + offset;

  /*--- Solver definition for the template problem ---*/
  if (template_solver) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Template()) {
      case SPACE_CENTERED : case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][TEMPLATE_SOL][conv_term] = new CConvective_Template(nDim, nVar_Template, config);
        break;
      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Convective scheme not implemented (template_solver).", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][TEMPLATE_SOL][visc_term] = new CViscous_Template(nDim, nVar_Template, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][TEMPLATE_SOL][source_first_term] = new CSource_Template(nDim, nVar_Template, config);

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TEMPLATE_SOL][conv_bound_term] = new CConvective_Template(nDim, nVar_Template, config);
    }

  }

  /*--- Solver definition for the Potential, Euler, Navier-Stokes problems ---*/
  if ((euler) || (ns)) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Flow()) {
      case NO_CONVECTIVE :
        SU2_OMP_MASTER
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_FLOW option.", CURRENT_FUNCTION);
        break;

      case SPACE_CENTERED :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentLax_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentJST_Flow(nDim, nVar_Flow, config); break;
            case JST_KE : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentJST_KE_Flow(nDim, nVar_Flow, config); break;
            default:
              SU2_OMP_MASTER
              SU2_MPI::Error("Invalid centered scheme or not implemented.", CURRENT_FUNCTION);
              break;
          }

          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_term] = new CCentLax_Flow(nDim, nVar_Flow, config);

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);

        }
        if (incompressible) {
          /*--- Incompressible flow, use preconditioning method ---*/
          switch (config->GetKind_Centered_Flow()) {
            case LAX : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentLaxInc_Flow(nDim, nVar_Flow, config); break;
            case JST : numerics[MESH_0][FLOW_SOL][conv_term] = new CCentJSTInc_Flow(nDim, nVar_Flow, config); break;
            default:
              SU2_OMP_MASTER
              SU2_MPI::Error("Invalid centered scheme or not implemented.\n Currently, only JST and LAX-FRIEDRICH are available for incompressible flows.", CURRENT_FUNCTION);
              break;
          }
          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_term] = new CCentLaxInc_Flow(nDim, nVar_Flow, config);

          /*--- Definition of the boundary condition method ---*/
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);

        }
        break;
      case SPACE_UPWIND :
        if (compressible) {
          /*--- Compressible flow ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case ROE:
              if (ideal_gas) {

                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwRoe_Flow(nDim, nVar_Flow, config, false);
                }
              } else {

                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwGeneralRoe_Flow(nDim, nVar_Flow, config);
                }
              }
              break;

            case AUSM:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
              }
              break;

            case AUSMPLUSUP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP_Flow(nDim, nVar_Flow, config);
              }
              break;

            case AUSMPLUSUP2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSMPLUSUP2_Flow(nDim, nVar_Flow, config);
              }
              break;

            case TURKEL:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
              }
              break;

            case L2ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwL2Roe_Flow(nDim, nVar_Flow, config);
              }
              break;
            case LMROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwLMRoe_Flow(nDim, nVar_Flow, config);
              }
              break;

            case SLAU:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwSLAU_Flow(nDim, nVar_Flow, config, false);
              }
              break;

            case SLAU2:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, roe_low_dissipation);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwSLAU2_Flow(nDim, nVar_Flow, config, false);
              }
              break;

            case HLLC:
              if (ideal_gas) {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              else {
                for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                  numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                  numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwGeneralHLLC_Flow(nDim, nVar_Flow, config);
                }
              }
              break;

            case MSW:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
              }
              break;

            case CUSP:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
              }
              break;

            default:
              SU2_OMP_MASTER
              SU2_MPI::Error("Invalid upwind scheme or not implemented.", CURRENT_FUNCTION);
              break;
          }

        }
        if (incompressible) {
          /*--- Incompressible flow, use artificial compressibility method ---*/
          switch (config->GetKind_Upwind_Flow()) {
            case FDS:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
                numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwFDSInc_Flow(nDim, nVar_Flow, config);
              }
              break;
            default:
              SU2_OMP_MASTER
              SU2_MPI::Error("Invalid upwind scheme or not implemented.\n Currently, only FDS is available for incompressible flows.", CURRENT_FUNCTION);
              break;
          }
        }
        break;

      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Invalid convective scheme for the Euler / Navier-Stokes equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    if (compressible) {
      if (ideal_gas) {

        /*--- Compressible flow Ideal gas ---*/
        numerics[MESH_0][FLOW_SOL][visc_term] = new CAvgGrad_Flow(nDim, nVar_Flow, true, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_term] = new CAvgGrad_Flow(nDim, nVar_Flow, false, config);

        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CAvgGrad_Flow(nDim, nVar_Flow, false, config);

      } else {

        /*--- Compressible flow Real gas ---*/
        numerics[MESH_0][FLOW_SOL][visc_term] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, true, config);
        for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_term] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, false, config);

        /*--- Definition of the boundary condition method ---*/
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CGeneralAvgGrad_Flow(nDim, nVar_Flow, false, config);

      }
    }
    if (incompressible) {
      /*--- Incompressible flow, use preconditioning method ---*/
      numerics[MESH_0][FLOW_SOL][visc_term] = new CAvgGradInc_Flow(nDim, nVar_Flow, true, config);
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics[iMGlevel][FLOW_SOL][visc_term] = new CAvgGradInc_Flow(nDim, nVar_Flow, false, config);

      /*--- Definition of the boundary condition method ---*/
      for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
        numerics[iMGlevel][FLOW_SOL][visc_bound_term] = new CAvgGradInc_Flow(nDim, nVar_Flow, false, config);
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      if (config->GetBody_Force() == YES) {
        if (incompressible)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncBodyForce(nDim, nVar_Flow, config);
        else
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceBodyForce(nDim, nVar_Flow, config);
      }
      else if (incompressible && (config->GetKind_DensityModel() == BOUSSINESQ)) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceBoussinesq(nDim, nVar_Flow, config);
      }
      else if (config->GetRotating_Frame() == YES) {
        if (incompressible)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncRotatingFrame_Flow(nDim, nVar_Flow, config);
        else
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceRotatingFrame_Flow(nDim, nVar_Flow, config);
      }
      else if (config->GetAxisymmetric() == YES) {
        if (incompressible)
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceIncAxisymmetric_Flow(nDim, nVar_Flow, config);
        else
          numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceAxisymmetric_Flow(nDim, nVar_Flow, config);
      }
      else if (config->GetGravityForce() == YES) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceGravity(nDim, nVar_Flow, config);
      }
      else if (config->GetWind_Gust() == YES) {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceWindGust(nDim, nVar_Flow, config);
      }
      else {
        numerics[iMGlevel][FLOW_SOL][source_first_term] = new CSourceNothing(nDim, nVar_Flow, config);
      }

      /*--- At the moment it is necessary to have the RHT equation in order to have a volumetric heat source. ---*/
      if (config->AddRadiation())
        numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceRadiation(nDim, nVar_Flow, config);
      else
        numerics[iMGlevel][FLOW_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Flow, config);
    }

  }

  /*--- Riemann solver definition for the Euler, Navier-Stokes problems for the FEM discretization. ---*/
  if ((fem_euler) || (fem_ns)) {

    switch (config->GetRiemann_Solver_FEM()) {
      case ROE:
      case LAX_FRIEDRICH:
        /* Hard coded optimized implementation is used in the DG solver. No need to allocate the
           corresponding entry in numerics. */
        break;

      case AUSM:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwAUSM_Flow(nDim, nVar_Flow, config);
        }
        break;

      case TURKEL:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwTurkel_Flow(nDim, nVar_Flow, config);
        }
        break;

      case HLLC:
          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
            numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
            numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwHLLC_Flow(nDim, nVar_Flow, config);
          }
        break;

      case MSW:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwMSW_Flow(nDim, nVar_Flow, config);
        }
        break;

      case CUSP:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][FLOW_SOL][conv_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
          numerics[iMGlevel][FLOW_SOL][conv_bound_term] = new CUpwCUSP_Flow(nDim, nVar_Flow, config);
        }
        break;

      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Riemann solver not implemented.", CURRENT_FUNCTION);
        break;
    }

  }

  /*--- Solver definition for the turbulent model problem ---*/

  if (turbulent) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/

    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NO_UPWIND:
        SU2_OMP_MASTER
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_TURB option.", CURRENT_FUNCTION);
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          if (spalart_allmaras || neg_spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras ) {
            numerics[iMGlevel][TURB_SOL][conv_term] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
          }
          else if (menter_sst) numerics[iMGlevel][TURB_SOL][conv_term] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        }
        break;
      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Invalid convective scheme for the turbulence equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras){
        numerics[iMGlevel][TURB_SOL][visc_term] = new CAvgGrad_TurbSA(nDim, nVar_Turb, true, config);
      }
      else if (neg_spalart_allmaras) numerics[iMGlevel][TURB_SOL][visc_term] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, true, config);
      else if (menter_sst) numerics[iMGlevel][TURB_SOL][visc_term] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, true, config);
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras) numerics[iMGlevel][TURB_SOL][source_first_term] = new CSourcePieceWise_TurbSA(nDim, nVar_Turb, config);
      else if (e_spalart_allmaras) numerics[iMGlevel][TURB_SOL][source_first_term] = new CSourcePieceWise_TurbSA_E(nDim, nVar_Turb, config);
      else if (comp_spalart_allmaras) numerics[iMGlevel][TURB_SOL][source_first_term] = new CSourcePieceWise_TurbSA_COMP(nDim, nVar_Turb, config);
      else if (e_comp_spalart_allmaras) numerics[iMGlevel][TURB_SOL][source_first_term] = new CSourcePieceWise_TurbSA_E_COMP(nDim, nVar_Turb, config);
      else if (neg_spalart_allmaras) numerics[iMGlevel][TURB_SOL][source_first_term] = new CSourcePieceWise_TurbSA_Neg(nDim, nVar_Turb, config);
      else if (menter_sst) numerics[iMGlevel][TURB_SOL][source_first_term] = new CSourcePieceWise_TurbSST(nDim, nVar_Turb, constants, kine_Inf, omega_Inf, config);
      numerics[iMGlevel][TURB_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Turb, config);
    }

    /*--- Definition of the boundary condition method ---*/

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      if (spalart_allmaras || e_spalart_allmaras || comp_spalart_allmaras || e_comp_spalart_allmaras) {
        numerics[iMGlevel][TURB_SOL][conv_bound_term] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics[iMGlevel][TURB_SOL][visc_bound_term] = new CAvgGrad_TurbSA(nDim, nVar_Turb, false, config);
      }
      else if (neg_spalart_allmaras) {
        numerics[iMGlevel][TURB_SOL][conv_bound_term] = new CUpwSca_TurbSA(nDim, nVar_Turb, config);
        numerics[iMGlevel][TURB_SOL][visc_bound_term] = new CAvgGrad_TurbSA_Neg(nDim, nVar_Turb, false, config);
      }
      else if (menter_sst) {
        numerics[iMGlevel][TURB_SOL][conv_bound_term] = new CUpwSca_TurbSST(nDim, nVar_Turb, config);
        numerics[iMGlevel][TURB_SOL][visc_bound_term] = new CAvgGrad_TurbSST(nDim, nVar_Turb, constants, false, config);
      }
    }
  }

  /*--- Solver definition for the transition model problem ---*/
  if (transition) {

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_Turb()) {
      case NO_UPWIND:
        SU2_OMP_MASTER
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_TURB option.", CURRENT_FUNCTION);
        break;
      case SPACE_UPWIND:
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
          numerics[iMGlevel][TRANS_SOL][conv_term] = new CUpwSca_TransLM(nDim, nVar_Trans, config);
        }
        break;
      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Invalid convective scheme for the transition equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TRANS_SOL][visc_term] = new CAvgGradCorrected_TransLM(nDim, nVar_Trans, config);
    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TRANS_SOL][source_first_term] = new CSourcePieceWise_TransLM(nDim, nVar_Trans, config);
      numerics[iMGlevel][TRANS_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Trans, config);
    }

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][TRANS_SOL][conv_bound_term] = new CUpwLin_TransLM(nDim, nVar_Trans, config);
    }
  }

  /*--- Solver definition of the finite volume heat solver  ---*/
  if (heat) {

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      numerics[iMGlevel][HEAT_SOL][visc_term] = new CAvgGradCorrected_Heat(nDim, nVar_Heat, config);
      numerics[iMGlevel][HEAT_SOL][visc_bound_term] = new CAvgGrad_Heat(nDim, nVar_Heat, config);

      switch (config->GetKind_ConvNumScheme_Heat()) {

        case SPACE_UPWIND :
          numerics[iMGlevel][HEAT_SOL][conv_term] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          numerics[iMGlevel][HEAT_SOL][conv_bound_term] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          break;

        case SPACE_CENTERED :
          numerics[iMGlevel][HEAT_SOL][conv_term] = new CCentSca_Heat(nDim, nVar_Heat, config);
          numerics[iMGlevel][HEAT_SOL][conv_bound_term] = new CUpwSca_Heat(nDim, nVar_Heat, config);
          break;

        default:
          SU2_OMP_MASTER
          SU2_MPI::Error("Invalid convective scheme for the heat transfer equations.", CURRENT_FUNCTION);
          break;
      }
    }
  }

  /*--- Solver definition for the radiation model problem ---*/

  if (config->AddRadiation()) {
    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    numerics[MESH_0][RAD_SOL][VISC_TERM] = new CAvgGradCorrected_P1(nDim, nVar_Rad, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    numerics[MESH_0][RAD_SOL][SOURCE_FIRST_TERM] = new CSourceP1(nDim, nVar_Rad, config);

    /*--- Definition of the boundary condition method ---*/
    numerics[MESH_0][RAD_SOL][VISC_BOUND_TERM] = new CAvgGradCorrected_P1(nDim, nVar_Rad, config);
  }

  /*--- Solver definition for the flow adjoint problem ---*/

  if (adj_euler || adj_ns) {

    if (incompressible)
      SU2_OMP_MASTER
      SU2_MPI::Error("Convective schemes not implemented for incompressible continuous adjoint.", CURRENT_FUNCTION);

    /*--- Definition of the convective scheme for each equation and mesh level ---*/

    switch (config->GetKind_ConvNumScheme_AdjFlow()) {
      case NO_CONVECTIVE:
        SU2_OMP_MASTER
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_ADJFLOW option.", CURRENT_FUNCTION);
        break;

      case SPACE_CENTERED :

        if (compressible) {

          /*--- Compressible flow ---*/

          switch (config->GetKind_Centered_AdjFlow()) {
            case LAX : numerics[MESH_0][ADJFLOW_SOL][conv_term] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            case JST : numerics[MESH_0][ADJFLOW_SOL][conv_term] = new CCentJST_AdjFlow(nDim, nVar_Adj_Flow, config); break;
            default:
              SU2_OMP_MASTER
              SU2_MPI::Error("Centered scheme not implemented.", CURRENT_FUNCTION);
              break;
          }

          for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][ADJFLOW_SOL][conv_term] = new CCentLax_AdjFlow(nDim, nVar_Adj_Flow, config);

          for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
            numerics[iMGlevel][ADJFLOW_SOL][conv_bound_term] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);

        }
        break;

      case SPACE_UPWIND :

        if (compressible) {

          /*--- Compressible flow ---*/

          switch (config->GetKind_Upwind_AdjFlow()) {
            case ROE:
              for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
                numerics[iMGlevel][ADJFLOW_SOL][conv_term] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
                numerics[iMGlevel][ADJFLOW_SOL][conv_bound_term] = new CUpwRoe_AdjFlow(nDim, nVar_Adj_Flow, config);
              }
              break;
            default:
              SU2_OMP_MASTER
              SU2_MPI::Error("Upwind scheme not implemented.", CURRENT_FUNCTION);
              break;
          }
        }
        break;

      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Invalid convective scheme for the continuous adjoint Euler / Navier-Stokes equations.", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/

    if (compressible) {

      /*--- Compressible flow ---*/

      numerics[MESH_0][ADJFLOW_SOL][visc_term] = new CAvgGradCorrected_AdjFlow(nDim, nVar_Adj_Flow, config);
      numerics[MESH_0][ADJFLOW_SOL][visc_bound_term] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);

      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
        numerics[iMGlevel][ADJFLOW_SOL][visc_term] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
        numerics[iMGlevel][ADJFLOW_SOL][visc_bound_term] = new CAvgGrad_AdjFlow(nDim, nVar_Adj_Flow, config);
      }

    }

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/

    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

      /*--- Note that RANS is incompatible with Axisymmetric or Rotational (Fix it!) ---*/

      if (compressible) {

        if (adj_ns) {

          numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceViscous_AdjFlow(nDim, nVar_Adj_Flow, config);

          if (config->GetRotating_Frame() == YES)
            numerics[iMGlevel][ADJFLOW_SOL][source_second_term] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics[iMGlevel][ADJFLOW_SOL][source_second_term] = new CSourceConservative_AdjFlow(nDim, nVar_Adj_Flow, config);

        }

        else {

          if (config->GetRotating_Frame() == YES)
            numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceRotatingFrame_AdjFlow(nDim, nVar_Adj_Flow, config);
          else if (config->GetAxisymmetric() == YES)
            numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceAxisymmetric_AdjFlow(nDim, nVar_Adj_Flow, config);
          else
            numerics[iMGlevel][ADJFLOW_SOL][source_first_term] = new CSourceNothing(nDim, nVar_Adj_Flow, config);

          numerics[iMGlevel][ADJFLOW_SOL][source_second_term] = new CSourceNothing(nDim, nVar_Adj_Flow, config);

        }

      }

    }

  }

  /*--- Solver definition for the turbulent adjoint problem ---*/
  if (adj_turb) {

    if (!spalart_allmaras)
      SU2_OMP_MASTER
      SU2_MPI::Error("Only the SA turbulence model can be used with the continuous adjoint solver.", CURRENT_FUNCTION);

    /*--- Definition of the convective scheme for each equation and mesh level ---*/
    switch (config->GetKind_ConvNumScheme_AdjTurb()) {
      case NO_CONVECTIVE:
        SU2_OMP_MASTER
        SU2_MPI::Error("Config file is missing the CONV_NUM_METHOD_ADJTURB option.", CURRENT_FUNCTION);
        break;
      case SPACE_UPWIND :
        for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
          numerics[iMGlevel][ADJTURB_SOL][conv_term] = new CUpwSca_AdjTurb(nDim, nVar_Adj_Turb, config);
        break;
      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Convective scheme not implemented (adjoint turbulence).", CURRENT_FUNCTION);
        break;
    }

    /*--- Definition of the viscous scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][ADJTURB_SOL][visc_term] = new CAvgGradCorrected_AdjTurb(nDim, nVar_Adj_Turb, config);

    /*--- Definition of the source term integration scheme for each equation and mesh level ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {
      numerics[iMGlevel][ADJTURB_SOL][source_first_term] = new CSourcePieceWise_AdjTurb(nDim, nVar_Adj_Turb, config);
      numerics[iMGlevel][ADJTURB_SOL][source_second_term] = new CSourceConservative_AdjTurb(nDim, nVar_Adj_Turb, config);
    }

    /*--- Definition of the boundary condition method ---*/
    for (iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++)
      numerics[iMGlevel][ADJTURB_SOL][conv_bound_term] = new CUpwLin_AdjTurb(nDim, nVar_Adj_Turb, config);

  }

  /*--- Numerics definition for FEM-like problems. ---*/

  if (fem) {
    /*--- Initialize the container for FEA_TERM. This will be the only one for most of the cases. ---*/
    switch (config->GetGeometricConditions()) {
      case SMALL_DEFORMATIONS:
        switch (config->GetMaterialModel()) {
          case LINEAR_ELASTIC:
            numerics[MESH_0][FEA_SOL][fea_term] = new CFEALinearElasticity(nDim, nVar_FEM, config);
            break;
          case NEO_HOOKEAN:
            SU2_OMP_MASTER
            SU2_MPI::Error("Material model does not correspond to geometric conditions.", CURRENT_FUNCTION);
            break;
          default:
            SU2_OMP_MASTER
            SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            break;
        }
        break;
      case LARGE_DEFORMATIONS :
        switch (config->GetMaterialModel()) {
          case LINEAR_ELASTIC:
            SU2_OMP_MASTER
            SU2_MPI::Error("Material model does not correspond to geometric conditions.", CURRENT_FUNCTION);
            break;
          case NEO_HOOKEAN:
            if (config->GetMaterialCompressibility() == COMPRESSIBLE_MAT) {
              numerics[MESH_0][FEA_SOL][fea_term] = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config);
            } else {
              SU2_OMP_MASTER
              SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            }
            break;
          case KNOWLES:
            if (config->GetMaterialCompressibility() == NEARLY_INCOMPRESSIBLE_MAT) {
              numerics[MESH_0][FEA_SOL][fea_term] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config);
            } else {
              SU2_OMP_MASTER
              SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            }
            break;
          case IDEAL_DE:
            if (config->GetMaterialCompressibility() == NEARLY_INCOMPRESSIBLE_MAT) {
              numerics[MESH_0][FEA_SOL][fea_term] = new CFEM_IdealDE(nDim, nVar_FEM, config);
            } else {
              SU2_OMP_MASTER
              SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            }
            break;
          default:
            SU2_OMP_MASTER
            SU2_MPI::Error("Material model not implemented.", CURRENT_FUNCTION);
            break;
        }
        break;
      default:
        SU2_OMP_MASTER
        SU2_MPI::Error("Solver not implemented.", CURRENT_FUNCTION);
        break;
    }

    /*--- The following definitions only make sense if we have a non-linear solution. ---*/
    if (config->GetGeometricConditions() == LARGE_DEFORMATIONS) {

      /*--- This allocates a container for electromechanical effects. ---*/

      bool de_effects = config->GetDE_Effects();
      if (de_effects)
        numerics[MESH_0][FEA_SOL][DE_TERM+offset] = new CFEM_DielectricElastomer(nDim, nVar_FEM, config);

      ifstream properties_file;

      string filename = config->GetFEA_FileName();
      if (nZone > 1)
        filename = config->GetMultizone_FileName(filename, iZone, ".dat");

      properties_file.open(filename.data(), ios::in);

      /*--- In case there is a properties file, containers are allocated for a number of material models. ---*/

      if (!(properties_file.fail())) {
        numerics[MESH_0][FEA_SOL][MAT_NHCOMP+offset]  = new CFEM_NeoHookean_Comp(nDim, nVar_FEM, config);
        numerics[MESH_0][FEA_SOL][MAT_IDEALDE+offset] = new CFEM_IdealDE(nDim, nVar_FEM, config);
        numerics[MESH_0][FEA_SOL][MAT_KNOWLES+offset] = new CFEM_Knowles_NearInc(nDim, nVar_FEM, config);
      }
    }
  }

  /*--- Instantiate the numerics for the mesh solver. ---*/
  if (config->GetDeform_Mesh())
    numerics[MESH_0][MESH_SOL][fea_term] = new CFEAMeshElasticity(nDim, nDim, geometry[MESH_0]->GetnElem(), config);

  } // end "per-thread" allocation loop

}

void CDriver::Numerics_Postprocessing(CNumerics *****numerics, CSolver***, CGeometry**,
                                      CConfig *config, unsigned short val_iInst) {

  for (unsigned short iMGlevel = 0; iMGlevel <= config->GetnMGLevels(); iMGlevel++) {

    for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++) {

      for (unsigned int iTerm = 0; iTerm < MAX_TERMS*omp_get_max_threads(); iTerm++) {

        delete numerics[val_iInst][iMGlevel][iSol][iTerm];
      }
      delete [] numerics[val_iInst][iMGlevel][iSol];
    }
    delete[] numerics[val_iInst][iMGlevel];
  }
  delete[] numerics[val_iInst];

}

void CDriver::Iteration_Preprocessing(CConfig* config, CIteration *&iteration) {

  if (rank == MASTER_NODE)
    cout << endl <<"------------------- Iteration Preprocessing ( Zone " << config->GetiZone() <<" ) ------------------" << endl;

  /*--- Loop over all zones and instantiate the physics iteration. ---*/

  switch (config->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
    case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
      if(config->GetBoolTurbomachinery()){
        if (rank == MASTER_NODE)
          cout << "Euler/Navier-Stokes/RANS turbomachinery fluid iteration." << endl;
        iteration = new CTurboIteration(config);

      }
      else{
        if (rank == MASTER_NODE)
          cout << "Euler/Navier-Stokes/RANS fluid iteration." << endl;
        iteration = new CFluidIteration(config);
      }
      break;

    case FEM_EULER: case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
      if (rank == MASTER_NODE)
        cout << "Finite element Euler/Navier-Stokes/RANS/LES flow iteration." << endl;
      iteration = new CFEMFluidIteration(config);
      break;

    case HEAT_EQUATION:
      if (rank == MASTER_NODE)
        cout << "Heat iteration (finite volume method)." << endl;
      iteration = new CHeatIteration(config);
      break;

    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << "FEM iteration." << endl;
      iteration = new CFEAIteration(config);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << "Adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration = new CAdjFluidIteration(config);
      break;

    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration = new CDiscAdjFluidIteration(config);
      break;

    case DISC_ADJ_FEM_EULER : case DISC_ADJ_FEM_NS : case DISC_ADJ_FEM_RANS :
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint finite element Euler/Navier-Stokes/RANS fluid iteration." << endl;
      iteration = new CDiscAdjFluidIteration(config);
      break;

    case DISC_ADJ_FEM:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint FEM structural iteration." << endl;
      iteration = new CDiscAdjFEAIteration(config);
      break;

    case DISC_ADJ_HEAT:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint heat iteration." << endl;
      iteration = new CDiscAdjHeatIteration(config);
      break;
  }
}

void CDriver::DynamicMesh_Preprocessing(CConfig *config, CGeometry **geometry, CSolver ***solver, CIteration* iteration,
                                        CVolumetricMovement *&grid_movement, CSurfaceMovement *&surface_movement){

  /*--- Instantiate the geometry movement classes for the solution of unsteady
   flows on dynamic meshes, including rigid mesh transformations, dynamically
   deforming meshes, and preprocessing of harmonic balance. ---*/

  if (!fem_solver && (config->GetGrid_Movement() ||
                      (config->GetDirectDiff() == D_DESIGN)) && !config->GetSurface_Movement(FLUID_STRUCTURE_STATIC)) {
    if (rank == MASTER_NODE)
      cout << "Setting dynamic mesh structure for zone "<< iZone + 1<<"." << endl;
    grid_movement = new CVolumetricMovement(geometry[MESH_0], config);

    surface_movement = new CSurfaceMovement();
    surface_movement->CopyBoundary(geometry[MESH_0], config);
    if (config->GetTime_Marching() == HARMONIC_BALANCE){
      if (rank == MASTER_NODE) cout << endl <<  "Instance "<< iInst + 1 <<":" << endl;
      iteration->SetGrid_Movement(geometry, surface_movement, grid_movement,  solver, config, 0, iInst);
    }
  }

  if (config->GetDirectDiff() == D_DESIGN) {
    if (rank == MASTER_NODE)
      cout << "Setting surface/volume derivatives." << endl;

    /*--- Set the surface derivatives, i.e. the derivative of the surface mesh nodes with respect to the design variables ---*/

    surface_movement->SetSurface_Derivative(geometry[MESH_0],config);

    /*--- Call the volume deformation routine with derivative mode enabled.
       This computes the derivative of the volume mesh with respect to the surface nodes ---*/


    grid_movement->SetVolume_Deformation(geometry[MESH_0],config, true, true);

    /*--- Update the multi-grid structure to propagate the derivative information to the coarser levels ---*/

    geometry[MESH_0]->UpdateGeometry(geometry,config);

    /*--- Set the derivative of the wall-distance with respect to the surface nodes ---*/

    if ( (config->GetKind_Solver() == RANS) ||
         (config->GetKind_Solver() == ADJ_RANS) ||
         (config->GetKind_Solver() == DISC_ADJ_RANS) ||
         (config->GetKind_Solver() == INC_RANS) ||
         (config->GetKind_Solver() == DISC_ADJ_INC_RANS))
      geometry[MESH_0]->ComputeWall_Distance(config);
  }


  if (config->GetSurface_Movement(FLUID_STRUCTURE_STATIC)){
    if (rank == MASTER_NODE)
      cout << "Setting moving mesh structure for FSI problems." << endl;
    /*--- Instantiate the container for the grid movement structure ---*/
    grid_movement = new CElasticityMovement(geometry[MESH_0], config);
  }

}

void CDriver::Interface_Preprocessing(CConfig **config, CSolver***** solver, CGeometry**** geometry,
                                      unsigned short** interface_types, CInterface ***&interface,
                                      CInterpolator ***&interpolation) {

  unsigned short donorZone, targetZone;
  unsigned short nVar, nVarTransfer;

  unsigned short nMarkerTarget, iMarkerTarget, nMarkerDonor, iMarkerDonor;

  /*--- Initialize some useful booleans ---*/
  bool fluid_donor, structural_donor, heat_donor;
  bool fluid_target, structural_target, heat_target;

  bool discrete_adjoint = config[ZONE_0]->GetDiscrete_Adjoint();

  int markDonor, markTarget, Donor_check, Target_check, iMarkerInt, nMarkerInt;

#ifdef HAVE_MPI
  int *Buffer_Recv_mark = NULL, iRank, nProcessor = size;

  if (rank == MASTER_NODE)
    Buffer_Recv_mark = new int[nProcessor];
#endif

  /*--- Coupling between zones ---*/
  // There's a limit here, the interface boundary must connect only 2 zones

  /*--- Loops over all target and donor zones to find which ones are connected through
   *--- an interface boundary (fsi or sliding mesh) ---*/
  for (targetZone = 0; targetZone < nZone; targetZone++) {

    for (donorZone = 0; donorZone < nZone; donorZone++) {

      interface_types[donorZone][targetZone] = NO_TRANSFER;

      if ( donorZone == targetZone ) {
        interface_types[donorZone][targetZone] = ZONES_ARE_EQUAL;
        // We're processing the same zone, so skip the following
        continue;
      }

      nMarkerInt = (int) ( config[donorZone]->GetMarker_n_ZoneInterface() / 2 );

      /*--- Loops on Interface markers to find if the 2 zones are sharing the boundary and to
       *--- determine donor and target marker tag ---*/
      for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

        markDonor  = -1;
        markTarget = -1;

        /*--- On the donor side ---*/
        nMarkerDonor = config[donorZone]->GetnMarker_All();

        for (iMarkerDonor = 0; iMarkerDonor < nMarkerDonor; iMarkerDonor++) {

          /*--- If the tag GetMarker_All_ZoneInterface(iMarker) equals the index we are looping at ---*/
          if ( config[donorZone]->GetMarker_All_ZoneInterface(iMarkerDonor) == iMarkerInt ) {
            /*--- We have identified the identifier for the interface marker ---*/
            markDonor = iMarkerDonor;

            break;
          }
        }

        /*--- On the target side ---*/
        nMarkerTarget = config[targetZone]->GetnMarker_All();

      for (iMarkerTarget = 0; iMarkerTarget < nMarkerTarget; iMarkerTarget++) {

          /*--- If the tag GetMarker_All_ZoneInterface(iMarker) equals the index we are looping at ---*/
        if ( config[targetZone]->GetMarker_All_ZoneInterface(iMarkerTarget) == iMarkerInt ) {
            /*--- We have identified the identifier for the interface marker ---*/
            markTarget = iMarkerTarget;

            break;
        }
        }

#ifdef HAVE_MPI

      Donor_check  = -1;
      Target_check = -1;

        /*--- We gather a vector in MASTER_NODE that determines if the boundary is not on the processor because
         * of the partition or because the zone does not include it ---*/

        SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

      if (rank == MASTER_NODE) {
        for (iRank = 0; iRank < nProcessor; iRank++) {
          if( Buffer_Recv_mark[iRank] != -1 ) {
              Donor_check = Buffer_Recv_mark[iRank];

              break;
            }
          }
        }

        SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

        SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

      if (rank == MASTER_NODE){
        for (iRank = 0; iRank < nProcessor; iRank++){
          if( Buffer_Recv_mark[iRank] != -1 ){
              Target_check = Buffer_Recv_mark[iRank];

              break;
            }
          }
        }

        SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

#else
      Donor_check  = markDonor;
      Target_check = markTarget;
#endif

      /* --- Check ifzones are actually sharing the interface boundary, if not skip ---*/
      if(Target_check == -1 || Donor_check == -1) {
        interface_types[donorZone][targetZone] = NO_COMMON_INTERFACE;
        continue;
      }

        /*--- Set some boolean to properly allocate data structure later ---*/
      fluid_target      = false;
      structural_target = false;

      fluid_donor       = false;
      structural_donor  = false;

      heat_donor        = false;
      heat_target       = false;

      switch ( config[targetZone]->GetKind_Solver() ) {

        case EULER : case NAVIER_STOKES: case RANS:
        case INC_EULER : case INC_NAVIER_STOKES: case INC_RANS:
        case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          fluid_target  = true;
          break;

        case FEM_ELASTICITY: case DISC_ADJ_FEM:
          structural_target = true;
          break;

        case HEAT_EQUATION: case DISC_ADJ_HEAT:
          heat_target = true;
          break;
      }

      switch ( config[donorZone]->GetKind_Solver() ) {

        case EULER : case NAVIER_STOKES: case RANS:
        case INC_EULER : case INC_NAVIER_STOKES: case INC_RANS:
        case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          fluid_donor  = true;
          break;

        case FEM_ELASTICITY: case DISC_ADJ_FEM:
          structural_donor = true;
          break;

        case HEAT_EQUATION : case DISC_ADJ_HEAT:
          heat_donor = true;
          break;
      }

      /*--- Begin the creation of the communication pattern among zones ---*/

      /*--- Retrieve the number of conservative variables (for problems not involving structural analysis ---*/
      if (fluid_donor && fluid_target)
        nVar = solver[donorZone][INST_0][MESH_0][FLOW_SOL]->GetnVar();
      else
        /*--- If at least one of the components is structural ---*/
        nVar = nDim;

      if (rank == MASTER_NODE) cout << "From zone " << donorZone << " to zone " << targetZone << ": ";

        /*--- Match Zones ---*/
      if (rank == MASTER_NODE) cout << "Setting coupling ";

          bool conservative_interp = config[donorZone]->GetConservativeInterpolation();

          /*--- Conditions for conservative interpolation are not met, we cannot fallback on the consistent approach
                because CFlowTractionInterface relies on the information in config to be correct. ---*/
          if ( conservative_interp && targetZone == 0 && structural_target )
            SU2_MPI::Error("Conservative interpolation assumes the structural model mesh is evaluated second, "
                           "somehow this has not happened.",CURRENT_FUNCTION);

        switch (config[donorZone]->GetKindInterpolation()) {

          case NEAREST_NEIGHBOR:
            if ( conservative_interp && targetZone > 0 && structural_target ) {
              interpolation[donorZone][targetZone] = new CMirror(geometry, config, donorZone, targetZone);
              if (rank == MASTER_NODE) cout << "using a mirror approach: matching coefficients "
                                               "from opposite mesh." << endl;
            }
            else {
            interpolation[donorZone][targetZone] = new CNearestNeighbor(geometry, config, donorZone, targetZone);
            if (rank == MASTER_NODE) cout << "using a nearest-neighbor approach." << endl;
            }
            break;

          case ISOPARAMETRIC:
            if ( conservative_interp && targetZone > 0 && structural_target ) {
              interpolation[donorZone][targetZone] = new CMirror(geometry, config, donorZone, targetZone);
              if (rank == MASTER_NODE) cout << "using a mirror approach: matching coefficients "
                                               "from opposite mesh." << endl;
            }
            else {
            interpolation[donorZone][targetZone] = new CIsoparametric(geometry, config, donorZone, targetZone);
            if (rank == MASTER_NODE) cout << "using an isoparametric approach." << endl;
            }
            break;

        case WEIGHTED_AVERAGE:
          interpolation[donorZone][targetZone] = new CSlidingMesh(geometry, config, donorZone, targetZone);
          if (rank == MASTER_NODE) cout << "using an sliding mesh approach." << endl;

          break;

          case RADIAL_BASIS_FUNCTION:
            if ( conservative_interp && targetZone > 0 && structural_target ) {
                interpolation[donorZone][targetZone] = new CMirror(geometry, config, donorZone, targetZone);
                if (rank == MASTER_NODE) cout << "using a mirror approach: matching coefficients "
                                                 "from opposite mesh." << endl;
              }
              else {
                interpolation[donorZone][targetZone] = new CRadialBasisFunction(geometry, config,
                                                                                donorZone, targetZone);
                if (rank == MASTER_NODE) cout << "using a radial basis function approach." << endl;
              }
            break;
            }

        /*--- Initialize the appropriate transfer strategy ---*/
      if (rank == MASTER_NODE) cout << "Transferring ";

      if (fluid_donor && structural_target) {
        interface_types[donorZone][targetZone] = FLOW_TRACTION;
        nVarTransfer = 2;
        if(!discrete_adjoint) {
          interface[donorZone][targetZone] = new CFlowTractionInterface(nVar, nVarTransfer, config[donorZone]);
        } else {
          interface[donorZone][targetZone] = new CDiscAdjFlowTractionInterface(nVar, nVarTransfer, config[donorZone]);
        }
        if (rank == MASTER_NODE) cout << "flow tractions. " << endl;
      }
      else if (structural_donor && fluid_target) {
        /*--- If we are using the new mesh solver, we transfer the total boundary displacements (not incremental) --*/
        if (solver_container[targetZone][INST_0][MESH_0][MESH_SOL] != NULL) {
          interface_types[donorZone][targetZone] = BOUNDARY_DISPLACEMENTS;
          nVarTransfer = 0;
          interface[donorZone][targetZone] = new CDisplacementsInterface(nVar, nVarTransfer, config[donorZone]);
          if (rank == MASTER_NODE) cout << "boundary displacements from the structural solver. " << endl;
        }
        /*--- We keep the legacy method temporarily until FSI-adjoint has been adapted ---*/
        /// TODO: LEGACY CLEANUP remove the "else" part and every class and enum referenced there,
        ///       add a check above to make sure MESH_SOL has been instantiated.
        else {
          nVarTransfer = 0;
          if(!discrete_adjoint) {
            interface_types[donorZone][targetZone] = STRUCTURAL_DISPLACEMENTS_LEGACY;
            interface[donorZone][targetZone] = new CDisplacementsInterfaceLegacy(nVar, nVarTransfer, config[donorZone]);
          } else {
            interface_types[donorZone][targetZone] = STRUCTURAL_DISPLACEMENTS_DISC_ADJ;
            interface[donorZone][targetZone] = new CDiscAdjDisplacementsInterfaceLegacy(nVar, nVarTransfer, config[donorZone]);
          }
          if (rank == MASTER_NODE) cout << "structural displacements (legacy). " << endl;
        }
      }
      else if (fluid_donor && fluid_target) {
        interface_types[donorZone][targetZone] = SLIDING_INTERFACE;
        nVarTransfer = 0;
        nVar = solver[donorZone][INST_0][MESH_0][FLOW_SOL]->GetnPrimVar();
        interface[donorZone][targetZone] = new CSlidingInterface(nVar, nVarTransfer, config[donorZone]);
        if (rank == MASTER_NODE) cout << "sliding interface. " << endl;
      }
      else if (fluid_donor && heat_target) {
        nVarTransfer = 0;
        nVar = 4;
        if(config[donorZone]->GetEnergy_Equation() || (config[donorZone]->GetKind_Regime() == COMPRESSIBLE))
          interface_types[donorZone][targetZone] = CONJUGATE_HEAT_FS;
        else if (config[donorZone]->GetWeakly_Coupled_Heat())
          interface_types[donorZone][targetZone] = CONJUGATE_HEAT_WEAKLY_FS;
        else { }
        interface[donorZone][targetZone] = new CConjugateHeatInterface(nVar, nVarTransfer, config[donorZone]);
        if (rank == MASTER_NODE) cout << "conjugate heat variables. " << endl;
      }
      else if (heat_donor && fluid_target) {
        nVarTransfer = 0;
        nVar = 4;
        if(config[targetZone]->GetEnergy_Equation() || (config[targetZone]->GetKind_Regime() == COMPRESSIBLE))
          interface_types[donorZone][targetZone] = CONJUGATE_HEAT_SF;
        else if (config[targetZone]->GetWeakly_Coupled_Heat())
          interface_types[donorZone][targetZone] = CONJUGATE_HEAT_WEAKLY_SF;
        else { }
        interface[donorZone][targetZone] = new CConjugateHeatInterface(nVar, nVarTransfer, config[donorZone]);
        if (rank == MASTER_NODE) cout << "conjugate heat variables. " << endl;
      }
      else if (heat_donor && heat_target) {
        SU2_MPI::Error("Conjugate heat transfer between solids not implemented yet.", CURRENT_FUNCTION);
      }
      else {
        interface_types[donorZone][targetZone] = CONSERVATIVE_VARIABLES;
        nVarTransfer = 0;
        interface[donorZone][targetZone] = new CConservativeVarsInterface(nVar, nVarTransfer, config[donorZone]);
        if (rank == MASTER_NODE) cout << "generic conservative variables. " << endl;
      }

      break;

      }

      if (config[donorZone]->GetBoolMixingPlaneInterface()){
        interface_types[donorZone][targetZone] = MIXING_PLANE;
        nVarTransfer = 0;
        nVar = solver[donorZone][INST_0][MESH_0][FLOW_SOL]->GetnVar();
        interface[donorZone][targetZone] = new CMixingPlaneInterface(nVar, nVarTransfer,
                                                                     config[donorZone], config[targetZone]);
        if (rank == MASTER_NODE) cout << "Set mixing-plane interface from donor zone "<< donorZone
                                      << " to target zone " << targetZone <<"."<<endl;
      }

    }

  }

#ifdef HAVE_MPI
  if (rank == MASTER_NODE)
  delete [] Buffer_Recv_mark;
#endif
}

void CDriver::StaticMesh_Preprocessing(CConfig *config, CGeometry** geometry, CSurfaceMovement* surface_movement){

  unsigned short iMGlevel, iMGfine;
  unsigned short Kind_Grid_Movement;

  unsigned short iZone = config->GetiZone();

  Kind_Grid_Movement = config->GetKind_GridMovement();

  if (!fem_solver) {

    switch (Kind_Grid_Movement) {

      case ROTATING_FRAME:

        /*--- Steadily rotating frame: set the grid velocities just once
         before the first iteration flow solver. ---*/

        if (rank == MASTER_NODE) {
          cout << endl << " Setting rotating frame grid velocities";
          cout << " for zone " << iZone << "." << endl;
        }

        /*--- Set the grid velocities on all multigrid levels for a steadily
           rotating reference frame. ---*/

        for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++){
          geometry[iMGlevel]->SetRotationalVelocity(config, true);
          geometry[iMGlevel]->SetShroudVelocity(config);
        }

        break;

      case STEADY_TRANSLATION:

        /*--- Set the translational velocity and hold the grid fixed during
         the calculation (similar to rotating frame, but there is no extra
         source term for translation). ---*/

        if (rank == MASTER_NODE)
          cout << endl << " Setting translational grid velocities." << endl;

        /*--- Set the translational velocity on all grid levels. ---*/

        for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++)
          geometry_container[iZone][INST_0][iMGlevel]->SetTranslationalVelocity(config, true);

        break;

      default:
        break;
    }

    if ((config->GetnMarker_Moving() > 0) && !config->GetSurface_Movement(FLUID_STRUCTURE_STATIC)) {

      /*--- Fixed wall velocities: set the grid velocities only one time
       before the first iteration flow solver. ---*/
      if (rank == MASTER_NODE)
        cout << endl << " Setting the moving wall velocities." << endl;

      assert(surface_movement != NULL && "A surface_movement was not instantiated.");
      surface_movement->Moving_Walls(geometry[MESH_0], config, iZone, 0);

      /*--- Update the grid velocities on the coarser multigrid levels after
        setting the moving wall velocities for the finest mesh. ---*/
      for (iMGlevel = 1; iMGlevel <= config->GetnMGLevels(); iMGlevel++){
        iMGfine = iMGlevel-1;
        geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine], config);
      }
    }
  } else {

    /*--- Carry out a dynamic cast to CMeshFEM_DG, such that it is not needed to
         define all virtual functions in the base class CGeometry. ---*/
    CMeshFEM_DG *DGMesh = dynamic_cast<CMeshFEM_DG *>(geometry[MESH_0]);

    /*--- Initialize the static mesh movement, if necessary. ---*/
    const unsigned short Kind_Grid_Movement = config->GetKind_GridMovement();
    const bool initStaticMovement = (config->GetGrid_Movement() &&
                                     (Kind_Grid_Movement == MOVING_WALL    ||
                                      Kind_Grid_Movement == ROTATING_FRAME ||
                                      Kind_Grid_Movement == STEADY_TRANSLATION));

    if(initStaticMovement){
      if (rank == MASTER_NODE) cout << "Initialize Static Mesh Movement" << endl;
      DGMesh->InitStaticMeshMovement(config, Kind_Grid_Movement, iZone);
    }
  }

}

void CDriver::Output_Preprocessing(CConfig **config, CConfig *driver_config, COutput **&output, COutput *&driver_output){

  /*--- Definition of the output class (one for each zone). The output class
   manages the writing of all restart, volume solution, surface solution,
   surface comma-separated value, and convergence history files (both in serial
   and in parallel). ---*/

  for (iZone = 0; iZone < nZone; iZone++){

    if (rank == MASTER_NODE)
      cout << endl <<"-------------------- Output Preprocessing ( Zone " << iZone <<" ) --------------------" << endl;

    /*--- Loop over all zones and instantiate the physics iteration. ---*/

    switch (config[iZone]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
      if (rank == MASTER_NODE)
        cout << "Euler/Navier-Stokes/RANS output structure." << endl;
      output[iZone] = new CFlowCompOutput(config[iZone], nDim);
      break;
    case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
      if (rank == MASTER_NODE)
        cout << "Euler/Navier-Stokes/RANS output structure." << endl;
      output[iZone] = new CFlowIncOutput(config[iZone], nDim);
      break;
    case HEAT_EQUATION:
      if (rank == MASTER_NODE)
        cout << "Heat output structure." << endl;
      output[iZone] = new CHeatOutput(config[iZone], nDim);
      break;
    case FEM_ELASTICITY:
      if (rank == MASTER_NODE)
        cout << "FEM output structure." << endl;
      output[iZone] = new CElasticityOutput(config[iZone], nDim);
      break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
      if (rank == MASTER_NODE)
        cout << "Adjoint Euler/Navier-Stokes/RANS output structure." << endl;
      output[iZone] = new CAdjFlowCompOutput(config[iZone], nDim);
      break;
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
      if (rank == MASTER_NODE)
        cout << "Adjoint Euler/Navier-Stokes/RANS output structure." << endl;
      output[iZone] = new CAdjFlowIncOutput(config[iZone], nDim);
      break;
    case DISC_ADJ_FEM:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint FEA output structure." << endl;
      output[iZone] = new CAdjElasticityOutput(config[iZone], nDim);
      break;

    case DISC_ADJ_HEAT:
      if (rank == MASTER_NODE)
        cout << "Discrete adjoint heat output structure." << endl;
      output[iZone] = new CAdjHeatOutput(config[iZone], nDim);
      break;

    case FEM_EULER: case FEM_LES: case FEM_RANS: case FEM_NAVIER_STOKES:
      if (rank == MASTER_NODE)
        cout << "FEM output structure." << endl;
      output[iZone] = new CFlowCompFEMOutput(config[iZone], nDim);
      break;

    default:
      if (rank == MASTER_NODE)
        cout << "Default output structure." << endl;
      output[iZone] = new COutput(config[iZone], nDim, false);
      break;
    }

    /*--- If dry-run is used, do not open/overwrite history file. ---*/
    output[iZone]->PreprocessHistoryOutput(config[iZone], !dry_run);

    output[iZone]->PreprocessVolumeOutput(config[iZone]);

  }

  if (driver_config->GetMultizone_Problem()){
    if (rank == MASTER_NODE)
      cout << endl <<"------------------- Output Preprocessing ( Multizone ) ------------------" << endl;

    driver_output = new CMultizoneOutput(driver_config, config, nDim);
    driver_output->PreprocessMultizoneHistoryOutput(output, config, driver_config, !dry_run);
  }


  /*--- Check for an unsteady restart. Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetTime_Domain() && config_container[ZONE_0]->GetRestart())
    TimeIter = config_container[ZONE_0]->GetRestart_Iter();

  /*--- Check for a dynamic restart (structural analysis). Update ExtIter if necessary. ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == FEM_ELASTICITY
      && config_container[ZONE_0]->GetWrt_Dynamic() && config_container[ZONE_0]->GetRestart())
    TimeIter = config_container[ZONE_0]->GetRestart_Iter();


}


void CDriver::Turbomachinery_Preprocessing(CConfig** config, CGeometry**** geometry, CSolver***** solver,
                                           CInterface*** interface){

  unsigned short donorZone,targetZone, nMarkerInt, iMarkerInt;
  unsigned short nSpanMax = 0;
  bool restart   = (config[ZONE_0]->GetRestart() || config[ZONE_0]->GetRestart_Flow());
  mixingplane = config[ZONE_0]->GetBoolMixingPlaneInterface();
  bool discrete_adjoint = config[ZONE_0]->GetDiscrete_Adjoint();
  su2double areaIn, areaOut, nBlades, flowAngleIn, flowAngleOut;

  /*--- Create turbovertex structure ---*/
  if (rank == MASTER_NODE) cout<<endl<<"Initialize Turbo Vertex Structure." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config[iZone]->GetBoolTurbomachinery()){
      geometry[iZone][INST_0][MESH_0]->ComputeNSpan(config[iZone], iZone, INFLOW, true);
      geometry[iZone][INST_0][MESH_0]->ComputeNSpan(config[iZone], iZone, OUTFLOW, true);
      if (rank == MASTER_NODE) cout <<"Number of span-wise sections in Zone "<< iZone<<": "<< config[iZone]->GetnSpanWiseSections() <<"."<< endl;
      if (config[iZone]->GetnSpanWiseSections() > nSpanMax){
        nSpanMax = config[iZone]->GetnSpanWiseSections();
      }

      config[ZONE_0]->SetnSpan_iZones(config[iZone]->GetnSpanWiseSections(), iZone);

      geometry[iZone][INST_0][MESH_0]->SetTurboVertex(config[iZone], iZone, INFLOW, true);
      geometry[iZone][INST_0][MESH_0]->SetTurboVertex(config[iZone], iZone, OUTFLOW, true);
    }
  }

  /*--- Set maximum number of Span among all zones ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if (config[iZone]->GetBoolTurbomachinery()){
      config[iZone]->SetnSpanMaxAllZones(nSpanMax);
    }
  }
  if (rank == MASTER_NODE) cout<<"Max number of span-wise sections among all zones: "<< nSpanMax<<"."<< endl;


  if (rank == MASTER_NODE) cout<<"Initialize solver containers for average and performance quantities." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->InitTurboContainers(geometry[iZone][INST_0][MESH_0],config[iZone]);
  }

//TODO(turbo) make it general for turbo HB
  if (rank == MASTER_NODE) cout<<"Compute inflow and outflow average geometric quantities." << endl;
  for (iZone = 0; iZone < nZone; iZone++) {
    geometry[iZone][INST_0][MESH_0]->SetAvgTurboValue(config[iZone], iZone, INFLOW, true);
    geometry[iZone][INST_0][MESH_0]->SetAvgTurboValue(config[iZone],iZone, OUTFLOW, true);
    geometry[iZone][INST_0][MESH_0]->GatherInOutAverageValues(config[iZone], true);
  }


  if(mixingplane){
    if (rank == MASTER_NODE) cout << "Set span-wise sections between zones on Mixing-Plane interface." << endl;
    for (donorZone = 0; donorZone < nZone; donorZone++) {
      for (targetZone = 0; targetZone < nZone; targetZone++) {
        if (targetZone != donorZone){
          interface[donorZone][targetZone]->SetSpanWiseLevels(config[donorZone], config[targetZone]);
        }
      }
    }
  }

  if (rank == MASTER_NODE) cout << "Transfer average geometric quantities to zone 0." << endl;
  for (iZone = 1; iZone < nZone; iZone++) {
    interface[iZone][ZONE_0]->GatherAverageTurboGeoValues(geometry[iZone][INST_0][MESH_0],geometry[ZONE_0][INST_0][MESH_0], iZone);
  }

  /*--- Transfer number of blade to ZONE_0 to correctly compute turbo performance---*/
  for (iZone = 1; iZone < nZone; iZone++) {
    nBlades = config[iZone]->GetnBlades(iZone);
    config[ZONE_0]->SetnBlades(iZone, nBlades);
  }

  if (rank == MASTER_NODE){
    for (iZone = 0; iZone < nZone; iZone++) {
    areaIn  = geometry[iZone][INST_0][MESH_0]->GetSpanAreaIn(iZone, config[iZone]->GetnSpanWiseSections());
    areaOut = geometry[iZone][INST_0][MESH_0]->GetSpanAreaOut(iZone, config[iZone]->GetnSpanWiseSections());
    nBlades = config[iZone]->GetnBlades(iZone);
    cout << "Inlet area for Row "<< iZone + 1<< ": " << areaIn*10000.0 <<" cm^2."  <<endl;
    cout << "Oulet area for Row "<< iZone + 1<< ": " << areaOut*10000.0 <<" cm^2."  <<endl;
    cout << "Recomputed number of blades for Row "<< iZone + 1 << ": " << nBlades<<"."  <<endl;
    }
  }


  if(mixingplane){
    if (rank == MASTER_NODE) cout<<"Preprocessing of the Mixing-Plane Interface." << endl;
    for (donorZone = 0; donorZone < nZone; donorZone++) {
      nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;
      for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
        for (targetZone = 0; targetZone < nZone; targetZone++) {
          if (targetZone != donorZone){
            interface[donorZone][targetZone]->PreprocessAverage(geometry[donorZone][INST_0][MESH_0], geometry[targetZone][INST_0][MESH_0],
                config[donorZone], config[targetZone],
                iMarkerInt);
          }
        }
      }
    }
  }

  if(!restart && !discrete_adjoint){
    if (rank == MASTER_NODE) cout<<"Initialize turbomachinery solution quantities." << endl;
    for(iZone = 0; iZone < nZone; iZone++) {
      solver[iZone][INST_0][MESH_0][FLOW_SOL]->SetFreeStream_TurboSolution(config[iZone]);
    }
  }

  if (rank == MASTER_NODE) cout<<"Initialize inflow and outflow average solution quantities." << endl;
  for(iZone = 0; iZone < nZone; iZone++) {
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->PreprocessAverage(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],INFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->PreprocessAverage(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],OUTFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->TurboAverageProcess(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],INFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->TurboAverageProcess(solver[iZone][INST_0][MESH_0], geometry[iZone][INST_0][MESH_0],config[iZone],OUTFLOW);
    solver[iZone][INST_0][MESH_0][FLOW_SOL]->GatherInOutAverageValues(config[iZone], geometry[iZone][INST_0][MESH_0]);
    if (rank == MASTER_NODE){
      flowAngleIn = solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityIn(iZone, config[iZone]->GetnSpanWiseSections())[1];
      flowAngleIn /= solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityIn(iZone, config[iZone]->GetnSpanWiseSections())[0];
      flowAngleIn = atan(flowAngleIn)*180.0/PI_NUMBER;
      cout << "Inlet flow angle for Row "<< iZone + 1<< ": "<< flowAngleIn <<"°."  <<endl;
      flowAngleOut = solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityOut(iZone, config[iZone]->GetnSpanWiseSections())[1];
      flowAngleOut /= solver[iZone][INST_0][MESH_0][FLOW_SOL]->GetTurboVelocityOut(iZone, config[iZone]->GetnSpanWiseSections())[0];
      flowAngleOut = atan(flowAngleOut)*180.0/PI_NUMBER;
      cout << "Outlet flow angle for Row "<< iZone + 1<< ": "<< flowAngleOut <<"°."  <<endl;

    }
  }

}




CDriver::~CDriver(void) {}

CFluidDriver::CFluidDriver(char* confFile, unsigned short val_nZone, SU2_Comm MPICommunicator) : CDriver(confFile, val_nZone, MPICommunicator, false) {
  Max_Iter = config_container[ZONE_0]->GetnInner_Iter();
}

CFluidDriver::~CFluidDriver(void) { }

void CFluidDriver::StartSolver(){

#ifdef VTUNEPROF
  __itt_resume();
#endif

  /*--- Main external loop of the solver. Within this loop, each iteration ---*/

  if (rank == MASTER_NODE)
    cout << endl <<"------------------------------ Begin Solver -----------------------------" << endl;

  unsigned long Iter = 0;
  while ( Iter < Max_Iter ) {

    /*--- Perform some external iteration preprocessing. ---*/

    Preprocess(Iter);

    /*--- Perform a dynamic mesh update if required. ---*/
    /*--- For the Disc.Adj. of a case with (rigidly) moving grid, the appropriate
          mesh cordinates are read from the restart files. ---*/
    if (!fem_solver &&
        !(config_container[ZONE_0]->GetGrid_Movement() && config_container[ZONE_0]->GetDiscrete_Adjoint())) {
      DynamicMeshUpdate(Iter);
    }

    /*--- Run a single iteration of the problem (fluid, elasticity, heat, ...). ---*/

    Run();

    /*--- Update the solution for dual time stepping strategy ---*/

    Update();

    /*--- Terminate the simulation if only the Jacobian must be computed. ---*/
    if (config_container[ZONE_0]->GetJacobian_Spatial_Discretization_Only()) break;

    /*--- Monitor the computations after each iteration. ---*/

    Monitor(Iter);

    /*--- Output the solution in files. ---*/

    Output(Iter);

    /*--- If the convergence criteria has been met, terminate the simulation. ---*/

    if (StopCalc) break;

    Iter++;

  }
#ifdef VTUNEPROF
  __itt_pause();
#endif
}


void CFluidDriver::Preprocess(unsigned long Iter) {

  /*--- Set the value of the external iteration and physical time. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    config_container[iZone]->SetInnerIter(Iter);
    if (config_container[iZone]->GetTime_Marching())
      config_container[iZone]->SetPhysicalTime(static_cast<su2double>(Iter)*config_container[iZone]->GetDelta_UnstTimeND());
    else
      config_container[iZone]->SetPhysicalTime(0.0);

  }


//  /*--- Read the target pressure ---*/

//  if (config_container[ZONE_0]->GetInvDesign_Cp() == YES)
//    output[ZONE_0]->SetCp_InverseDesign(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL],
//        geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], ExtIter);

//  /*--- Read the target heat flux ---*/

//  if (config_container[ZONE_0]->GetInvDesign_HeatFlux() == YES)
//    output[ZONE_0]->SetHeatFlux_InverseDesign(solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL],
//        geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0], ExtIter);

  /*--- Set the initial condition for EULER/N-S/RANS and for a non FSI simulation ---*/

  if(!fsi) {
    for (iZone = 0; iZone < nZone; iZone++) {
      if ((config_container[iZone]->GetKind_Solver() ==  EULER) ||
          (config_container[iZone]->GetKind_Solver() ==  NAVIER_STOKES) ||
          (config_container[iZone]->GetKind_Solver() ==  RANS) ||
          (config_container[iZone]->GetKind_Solver() ==  INC_EULER) ||
          (config_container[iZone]->GetKind_Solver() ==  INC_NAVIER_STOKES) ||
          (config_container[iZone]->GetKind_Solver() ==  INC_RANS)) {
        for (iInst = 0; iInst < nInst[iZone]; iInst++)
          solver_container[iZone][iInst][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[iZone][INST_0], solver_container[iZone][iInst], config_container[iZone], Iter);
      }
    }
  }

}

void CFluidDriver::Run() {

  unsigned short iZone, jZone, checkConvergence;
  unsigned long IntIter, nIntIter;
  bool unsteady;

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  unsteady = (config_container[MESH_0]->GetTime_Marching() == DT_STEPPING_1ST) || (config_container[MESH_0]->GetTime_Marching() == DT_STEPPING_2ND);

  /*--- Zone preprocessing ---*/

  for (iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);

  /*--- Updating zone interface communication patterns,
   needed only for unsteady simulation since for steady problems
   this is done once in the interpolator_container constructor
   at the beginning of the computation ---*/

  if ( unsteady ) {
    for (iZone = 0; iZone < nZone; iZone++) {
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && interpolator_container[iZone][jZone] != NULL)
        interpolator_container[iZone][jZone]->Set_TransferCoeff(config_container);
    }
  }

  /*--- Begin Unsteady pseudo-time stepping internal loop, if not unsteady it does only one step --*/

  if (unsteady)
    nIntIter = config_container[MESH_0]->GetnInner_Iter();
  else
    nIntIter = 1;

  for (IntIter = 0; IntIter < nIntIter; IntIter++) {

    /*--- At each pseudo time-step updates transfer data ---*/
    for (iZone = 0; iZone < nZone; iZone++)
      for (jZone = 0; jZone < nZone; jZone++)
        if(jZone != iZone && interface_container[iZone][jZone] != NULL)
          Transfer_Data(iZone, jZone);

    /*--- For each zone runs one single iteration ---*/

    for (iZone = 0; iZone < nZone; iZone++) {
      config_container[iZone]->SetInnerIter(IntIter);
      iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, iZone, INST_0);
    }

    /*--- Check convergence in each zone --*/

    checkConvergence = 0;
    for (iZone = 0; iZone < nZone; iZone++)
    checkConvergence += (int) integration_container[iZone][INST_0][FLOW_SOL]->GetConvergence();

    /*--- If convergence was reached in every zone --*/

  if (checkConvergence == nZone) break;
  }

}

void CFluidDriver::Transfer_Data(unsigned short donorZone, unsigned short targetZone) {

  interface_container[donorZone][targetZone]->BroadcastData(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
      geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
      config_container[donorZone], config_container[targetZone]);
  if (config_container[targetZone]->GetKind_Solver() == RANS)
    interface_container[donorZone][targetZone]->BroadcastData(solver_container[donorZone][INST_0][MESH_0][TURB_SOL],solver_container[targetZone][INST_0][MESH_0][TURB_SOL],
        geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
        config_container[donorZone], config_container[targetZone]);

}

void CFluidDriver::Update() {

  for(iZone = 0; iZone < nZone; iZone++)
    iteration_container[iZone][INST_0]->Update(output_container[iZone], integration_container, geometry_container,
         solver_container, numerics_container, config_container,
         surface_movement, grid_movement, FFDBox, iZone, INST_0);
}

void CFluidDriver::DynamicMeshUpdate(unsigned long TimeIter) {

  bool harmonic_balance;

  for (iZone = 0; iZone < nZone; iZone++) {
   harmonic_balance = (config_container[iZone]->GetTime_Marching() == HARMONIC_BALANCE);
    /*--- Dynamic mesh update ---*/
    if ((config_container[iZone]->GetGrid_Movement()) && (!harmonic_balance)) {
      iteration_container[iZone][INST_0]->SetGrid_Movement(geometry_container[iZone][INST_0], surface_movement[iZone], grid_movement[iZone][INST_0], solver_container[iZone][INST_0], config_container[iZone], 0, TimeIter );
    }
  }

}
bool CFluidDriver::Monitor(unsigned long ExtIter) {

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  IterCount++;
  UsedTime = (StopTime - StartTime) + UsedTimeCompute;


  /*--- Check if there is any change in the runtime parameters ---*/

  CConfig *runtime = NULL;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetTimeIter(ExtIter);
  delete runtime;

  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
    case EULER: case NAVIER_STOKES: case RANS:
      StopCalc = integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence(); break;
    case HEAT_EQUATION:
      StopCalc = integration_container[ZONE_0][INST_0][HEAT_SOL]->GetConvergence(); break;
    case FEM_ELASTICITY:
      StopCalc = integration_container[ZONE_0][INST_0][FEA_SOL]->GetConvergence(); break;
    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS:
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
    case DISC_ADJ_FEM_EULER: case DISC_ADJ_FEM_NS: case DISC_ADJ_FEM_RANS:
      StopCalc = integration_container[ZONE_0][INST_0][ADJFLOW_SOL]->GetConvergence(); break;
  }

  /*--- Set StopCalc to true if max. number of iterations has been reached ---*/

  StopCalc = StopCalc || (ExtIter == Max_Iter - 1);

  return StopCalc;

}


void CFluidDriver::Output(unsigned long InnerIter) {

  for (iZone = 0; iZone < nZone; iZone++) {
    const auto inst = config_container[iZone]->GetiInst();

    for (iInst = 0; iInst < nInst[iZone]; ++iInst) {
      config_container[iZone]->SetiInst(iInst);
      output_container[iZone]->SetResult_Files(geometry_container[iZone][iInst][MESH_0],
                                               config_container[iZone],
                                               solver_container[iZone][iInst][MESH_0],
                                               InnerIter, StopCalc);
    }
    config_container[iZone]->SetiInst(inst);
  }

}


CTurbomachineryDriver::CTurbomachineryDriver(char* confFile, unsigned short val_nZone,
                                             SU2_Comm MPICommunicator):
                                             CFluidDriver(confFile, val_nZone, MPICommunicator) {

  output_legacy = new COutputLegacy(config_container[ZONE_0]);

  /*--- LEGACY OUTPUT (going to be removed soon) --- */

  /*--- Open the convergence history file ---*/
  ConvHist_file = NULL;
  ConvHist_file = new ofstream*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    ConvHist_file[iZone] = NULL;
    if (rank == MASTER_NODE){
      ConvHist_file[iZone] = new ofstream[nInst[iZone]];
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        output_legacy->SetConvHistory_Header(&ConvHist_file[iZone][iInst], config_container[iZone], iZone, iInst);
        config_container[iZone]->SetHistFile(&ConvHist_file[iZone][INST_0]);
      }
    }
  }

  if (nZone > 1){
    Max_Iter = config_container[ZONE_0]->GetnOuter_Iter();
  }
}

CTurbomachineryDriver::~CTurbomachineryDriver(void) {
  if (rank == MASTER_NODE){
    /*--- Close the convergence history file. ---*/
    for (iZone = 0; iZone < nZone; iZone++) {
      for (iInst = 0; iInst < 1; iInst++) {
        ConvHist_file[iZone][iInst].close();
      }
      delete [] ConvHist_file[iZone];
    }
    delete [] ConvHist_file;
  }
}

void CTurbomachineryDriver::Run() {

  /*--- Run a single iteration of a multi-zone problem by looping over all
   zones and executing the iterations. Note that data transers between zones
   and other intermediate procedures may be required. ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Preprocess(output_container[iZone], integration_container, geometry_container,
                                           solver_container, numerics_container, config_container,
                                           surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  /* --- Update the mixing-plane interface ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    if(mixingplane)SetMixingPlane(iZone);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Iterate(output_container[iZone], integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  for (iZone = 0; iZone < nZone; iZone++) {
    iteration_container[iZone][INST_0]->Postprocess(output_container[iZone], integration_container, geometry_container,
                                      solver_container, numerics_container, config_container,
                                      surface_movement, grid_movement, FFDBox, iZone, INST_0);
  }

  if (rank == MASTER_NODE){
    SetTurboPerformance(ZONE_0);
  }


}

void CTurbomachineryDriver::SetMixingPlane(unsigned short donorZone){

  unsigned short targetZone, nMarkerInt, iMarkerInt ;
  nMarkerInt     = config_container[donorZone]->GetnMarker_MixingPlaneInterface()/2;

  /* --- transfer the average value from the donorZone to the targetZone*/
  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++){
    for (targetZone = 0; targetZone < nZone; targetZone++) {
      if (targetZone != donorZone){
        interface_container[donorZone][targetZone]->AllgatherAverage(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
            geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
            config_container[donorZone], config_container[targetZone], iMarkerInt );
      }
    }
  }
}

void CTurbomachineryDriver::SetTurboPerformance(unsigned short targetZone){

  unsigned short donorZone;
  //IMPORTANT this approach of multi-zone performances rely upon the fact that turbomachinery markers follow the natural (stator-rotor) development of the real machine.
  /* --- transfer the local turboperfomance quantities (for each blade)  from all the donorZones to the targetZone (ZONE_0) ---*/
  for (donorZone = 1; donorZone < nZone; donorZone++) {
    interface_container[donorZone][targetZone]->GatherAverageValues(solver_container[donorZone][INST_0][MESH_0][FLOW_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL], donorZone);
  }

  /* --- compute turboperformance for each stage and the global machine ---*/

 output_legacy->ComputeTurboPerformance(solver_container[targetZone][INST_0][MESH_0][FLOW_SOL], geometry_container[targetZone][INST_0][MESH_0], config_container[targetZone]);

}


bool CTurbomachineryDriver::Monitor(unsigned long ExtIter) {

  su2double rot_z_ini, rot_z_final ,rot_z;
  su2double outPres_ini, outPres_final, outPres;
  unsigned long rampFreq, finalRamp_Iter;
  unsigned short iMarker, KindBC, KindBCOption;
  string Marker_Tag;

  bool print;

  /*--- Synchronization point after a single solver iteration. Compute the
   wall clock time required. ---*/

#ifndef HAVE_MPI
  StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#else
  StopTime = MPI_Wtime();
#endif
  IterCount++;
  UsedTime = (StopTime - StartTime);


  /*--- Check if there is any change in the runtime parameters ---*/
  CConfig *runtime = NULL;
  strcpy(runtime_file_name, "runtime.dat");
  runtime = new CConfig(runtime_file_name, config_container[ZONE_0]);
  runtime->SetInnerIter(ExtIter);
  delete runtime;

  /*--- Update the convergence history file (serial and parallel computations). ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      output_legacy->SetConvHistory_Body(&ConvHist_file[iZone][iInst], geometry_container, solver_container,
          config_container, integration_container, false, UsedTime, iZone, iInst);
  }

  /*--- ROTATING FRAME Ramp: Compute the updated rotational velocity. ---*/
  if (config_container[ZONE_0]->GetGrid_Movement() && config_container[ZONE_0]->GetRampRotatingFrame()) {
    rampFreq       = SU2_TYPE::Int(config_container[ZONE_0]->GetRampRotatingFrame_Coeff(1));
    finalRamp_Iter = SU2_TYPE::Int(config_container[ZONE_0]->GetRampRotatingFrame_Coeff(2));
    rot_z_ini = config_container[ZONE_0]->GetRampRotatingFrame_Coeff(0);
    print = false;
    if(ExtIter % rampFreq == 0 &&  ExtIter <= finalRamp_Iter){

      for (iZone = 0; iZone < nZone; iZone++) {
        rot_z_final = config_container[iZone]->GetFinalRotation_Rate_Z();
        if(abs(rot_z_final) > 0.0){
          rot_z = rot_z_ini + ExtIter*( rot_z_final - rot_z_ini)/finalRamp_Iter;
          config_container[iZone]->SetRotation_Rate(2, rot_z);
          if(rank == MASTER_NODE && print && ExtIter > 0) {
            cout << endl << " Updated rotating frame grid velocities";
            cout << " for zone " << iZone << "." << endl;
          }
          geometry_container[iZone][INST_0][MESH_0]->SetRotationalVelocity(config_container[iZone], print);
          geometry_container[iZone][INST_0][MESH_0]->SetShroudVelocity(config_container[iZone]);
        }
      }

      for (iZone = 0; iZone < nZone; iZone++) {
        geometry_container[iZone][INST_0][MESH_0]->SetAvgTurboValue(config_container[iZone], iZone, INFLOW, false);
        geometry_container[iZone][INST_0][MESH_0]->SetAvgTurboValue(config_container[iZone],iZone, OUTFLOW, false);
        geometry_container[iZone][INST_0][MESH_0]->GatherInOutAverageValues(config_container[iZone], false);

      }

      for (iZone = 1; iZone < nZone; iZone++) {
        interface_container[iZone][ZONE_0]->GatherAverageTurboGeoValues(geometry_container[iZone][INST_0][MESH_0],geometry_container[ZONE_0][INST_0][MESH_0], iZone);
      }

    }
  }


  /*--- Outlet Pressure Ramp: Compute the updated rotational velocity. ---*/
  if (config_container[ZONE_0]->GetRampOutletPressure()) {
    rampFreq       = SU2_TYPE::Int(config_container[ZONE_0]->GetRampOutletPressure_Coeff(1));
    finalRamp_Iter = SU2_TYPE::Int(config_container[ZONE_0]->GetRampOutletPressure_Coeff(2));
    outPres_ini    = config_container[ZONE_0]->GetRampOutletPressure_Coeff(0);
    outPres_final  = config_container[ZONE_0]->GetFinalOutletPressure();

    if(ExtIter % rampFreq == 0 &&  ExtIter <= finalRamp_Iter){
      outPres = outPres_ini + ExtIter*(outPres_final - outPres_ini)/finalRamp_Iter;
      if(rank == MASTER_NODE) config_container[ZONE_0]->SetMonitotOutletPressure(outPres);

      for (iZone = 0; iZone < nZone; iZone++) {
        for (iMarker = 0; iMarker < config_container[iZone]->GetnMarker_All(); iMarker++) {
          KindBC = config_container[iZone]->GetMarker_All_KindBC(iMarker);
          switch (KindBC) {
          case RIEMANN_BOUNDARY:
            Marker_Tag         = config_container[iZone]->GetMarker_All_TagBound(iMarker);
            KindBCOption       = config_container[iZone]->GetKind_Data_Riemann(Marker_Tag);
            if(KindBCOption == STATIC_PRESSURE || KindBCOption == RADIAL_EQUILIBRIUM ){
              SU2_MPI::Error("Outlet pressure ramp only implemented for NRBC", CURRENT_FUNCTION);
            }
            break;
          case GILES_BOUNDARY:
            Marker_Tag         = config_container[iZone]->GetMarker_All_TagBound(iMarker);
            KindBCOption       = config_container[iZone]->GetKind_Data_Giles(Marker_Tag);
            if(KindBCOption == STATIC_PRESSURE || KindBCOption == STATIC_PRESSURE_1D || KindBCOption == RADIAL_EQUILIBRIUM ){
              config_container[iZone]->SetGiles_Var1(outPres, Marker_Tag);
            }
            break;
          }
        }
      }
    }
  }


  /*--- Check whether the current simulation has reached the specified
   convergence criteria, and set StopCalc to true, if so. ---*/

  switch (config_container[ZONE_0]->GetKind_Solver()) {
  case EULER: case NAVIER_STOKES: case RANS:
  case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
    StopCalc = integration_container[ZONE_0][INST_0][FLOW_SOL]->GetConvergence(); break;
  case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
  case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
  case DISC_ADJ_FEM_EULER: case DISC_ADJ_FEM_NS: case DISC_ADJ_FEM_RANS:
    StopCalc = integration_container[ZONE_0][INST_0][ADJFLOW_SOL]->GetConvergence(); break;
  }

  /*--- Set StopCalc to true if max. number of iterations has been reached ---*/

  StopCalc = StopCalc || (ExtIter == Max_Iter - 1);

  return StopCalc;

}

CHBDriver::CHBDriver(char* confFile,
    unsigned short val_nZone,
    SU2_Comm MPICommunicator) : CFluidDriver(confFile,
        val_nZone,
        MPICommunicator) {
  unsigned short kInst;

  nInstHB = nInst[ZONE_0];

  D = NULL;
  /*--- allocate dynamic memory for the Harmonic Balance operator ---*/
  D = new su2double*[nInstHB]; for (kInst = 0; kInst < nInstHB; kInst++) D[kInst] = new su2double[nInstHB];

  output_legacy = new COutputLegacy(config_container[ZONE_0]);

  /*--- Open the convergence history file ---*/
  ConvHist_file = NULL;
  ConvHist_file = new ofstream*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    ConvHist_file[iZone] = NULL;
    if (rank == MASTER_NODE){
      ConvHist_file[iZone] = new ofstream[nInst[iZone]];
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        output_legacy->SetConvHistory_Header(&ConvHist_file[iZone][iInst], config_container[iZone], iZone, iInst);
        config_container[iZone]->SetHistFile(&ConvHist_file[iZone][iInst]);
      }
    }
  }


}

CHBDriver::~CHBDriver(void) {

  unsigned short kInst;

  /*--- delete dynamic memory for the Harmonic Balance operator ---*/
  for (kInst = 0; kInst < nInstHB; kInst++) if (D[kInst] != NULL) delete [] D[kInst];
  if (D[kInst] != NULL) delete [] D;

  if (rank == MASTER_NODE){
  /*--- Close the convergence history file. ---*/
  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInstHB; iInst++) {
      ConvHist_file[iZone][iInst].close();
    }
    delete [] ConvHist_file[iZone];
  }
  delete [] ConvHist_file;
  }
}


void CHBDriver::Run() {

  /*--- Run a single iteration of a Harmonic Balance problem. Preprocess all
   all zones before beginning the iteration. ---*/

  for (iInst = 0; iInst < nInstHB; iInst++)
    iteration_container[ZONE_0][iInst]->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  for (iInst = 0; iInst < nInstHB; iInst++)
    iteration_container[ZONE_0][iInst]->Iterate(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  /*--- Update the convergence history file (serial and parallel computations). ---*/

  for (iZone = 0; iZone < nZone; iZone++) {
    for (iInst = 0; iInst < nInst[iZone]; iInst++)
      output_legacy->SetConvHistory_Body(&ConvHist_file[iZone][iInst], geometry_container, solver_container,
          config_container, integration_container, false, UsedTime, iZone, iInst);
  }

}

void CHBDriver::Update() {

  for (iInst = 0; iInst < nInstHB; iInst++) {
    /*--- Compute the harmonic balance terms across all zones ---*/
    SetHarmonicBalance(iInst);

  }

  /*--- Precondition the harmonic balance source terms ---*/
  if (config_container[ZONE_0]->GetHB_Precondition() == YES) {
    StabilizeHarmonicBalance();

  }

  for (iInst = 0; iInst < nInstHB; iInst++) {

    /*--- Update the harmonic balance terms across all zones ---*/
    iteration_container[ZONE_0][iInst]->Update(output_container[ZONE_0], integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);

  }

}

void CHBDriver::ResetConvergence() {

  for(iInst = 0; iInst < nZone; iInst++) {
    switch (config_container[ZONE_0]->GetKind_Solver()) {

    case EULER: case NAVIER_STOKES: case RANS:
      integration_container[ZONE_0][iInst][FLOW_SOL]->SetConvergence(false);
      if (config_container[ZONE_0]->GetKind_Solver() == RANS) integration_container[ZONE_0][iInst][TURB_SOL]->SetConvergence(false);
      if(config_container[ZONE_0]->GetKind_Trans_Model() == LM) integration_container[ZONE_0][iInst][TRANS_SOL]->SetConvergence(false);
      break;

    case FEM_ELASTICITY:
      integration_container[ZONE_0][iInst][FEA_SOL]->SetConvergence(false);
      break;

    case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      integration_container[ZONE_0][iInst][ADJFLOW_SOL]->SetConvergence(false);
      if( (config_container[ZONE_0]->GetKind_Solver() == ADJ_RANS) || (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS) )
        integration_container[ZONE_0][iInst][ADJTURB_SOL]->SetConvergence(false);
      break;
    }
  }

}

void CHBDriver::SetHarmonicBalance(unsigned short iInst) {

  unsigned short iVar, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());
  if (adjoint) {
    implicit = (config_container[ZONE_0]->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  }

  unsigned long InnerIter = config_container[ZONE_0]->GetInnerIter();

  /*--- Retrieve values from the config file ---*/
  su2double *U = new su2double[nVar];
  su2double *U_old = new su2double[nVar];
  su2double *Psi = new su2double[nVar];
  su2double *Psi_old = new su2double[nVar];
  su2double *Source = new su2double[nVar];
  su2double deltaU, deltaPsi;

  /*--- Compute period of oscillation ---*/
  su2double period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.  */
  period /= config_container[ZONE_0]->GetTime_Ref();

  if (InnerIter == 0)
    ComputeHB_Operator();

  /*--- Compute various source terms for explicit direct, implicit direct, and adjoint problems ---*/
  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][iInst][iMGlevel]->GetnPoint(); iPoint++) {

      for (iVar = 0; iVar < nVar; iVar++) {
        Source[iVar] = 0.0;
      }

      /*--- Step across the columns ---*/
      for (jInst = 0; jInst < nInstHB; jInst++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar; iVar++) {

          if (!adjoint) {
            U[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
            Source[iVar] += U[iVar]*D[iInst][jInst];

            if (implicit) {
              U_old[iVar] = solver_container[ZONE_0][jInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
              deltaU = U[iVar] - U_old[iVar];
              Source[iVar] += deltaU*D[iInst][jInst];
            }

          }

          else {
            Psi[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
            Source[iVar] += Psi[iVar]*D[jInst][iInst];

            if (implicit) {
              Psi_old[iVar] = solver_container[ZONE_0][jInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
              deltaPsi = Psi[iVar] - Psi_old[iVar];
              Source[iVar] += deltaPsi*D[jInst][iInst];
            }
          }
        }

        /*--- Store sources for current row ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (!adjoint) {
            solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iVar]);
          }
          else {
            solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iVar]);
          }
        }

      }
    }
  }

  /*--- Source term for a turbulence model ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == RANS) {

    /*--- Extra variables needed if we have a turbulence model. ---*/
    unsigned short nVar_Turb = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar();
    su2double *U_Turb = new su2double[nVar_Turb];
    su2double *Source_Turb = new su2double[nVar_Turb];

    /*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar_Turb; iVar++) Source_Turb[iVar] = 0.0;
      for (jInst = 0; jInst < nInstHB; jInst++) {

        /*--- Retrieve solution at this node in current zone ---*/
        for (iVar = 0; iVar < nVar_Turb; iVar++) {
          U_Turb[iVar] = solver_container[ZONE_0][jInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint, iVar);
          Source_Turb[iVar] += U_Turb[iVar]*D[iInst][jInst];
        }
      }

      /*--- Store sources for current iZone ---*/
      for (iVar = 0; iVar < nVar_Turb; iVar++)
        solver_container[ZONE_0][iInst][MESH_0][TURB_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source_Turb[iVar]);
    }

    delete [] U_Turb;
    delete [] Source_Turb;
  }

  delete [] Source;
  delete [] U;
  delete [] U_old;
  delete [] Psi;
  delete [] Psi_old;

}

void CHBDriver::StabilizeHarmonicBalance() {

  unsigned short i, j, k, iVar, iInst, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  bool adjoint = (config_container[ZONE_0]->GetContinuous_Adjoint());

  /*--- Retrieve values from the config file ---*/
  su2double *Source     = new su2double[nInstHB];
  su2double *Source_old = new su2double[nInstHB];
  su2double Delta;

  su2double **Pinv     = new su2double*[nInstHB];
  su2double **P        = new su2double*[nInstHB];
  for (iInst = 0; iInst < nInstHB; iInst++) {
    Pinv[iInst]       = new su2double[nInstHB];
    P[iInst]          = new su2double[nInstHB];
  }

  /*--- Loop over all grid levels ---*/
  for (iMGlevel = 0; iMGlevel <= config_container[ZONE_0]->GetnMGLevels(); iMGlevel++) {

    /*--- Loop over each node in the volume mesh ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][iMGlevel]->GetnPoint(); iPoint++) {

      /*--- Get time step for current node ---*/
      Delta = solver_container[ZONE_0][INST_0][iMGlevel][FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint);

      /*--- Setup stabilization matrix for this node ---*/
      for (iInst = 0; iInst < nInstHB; iInst++) {
        for (jInst = 0; jInst < nInstHB; jInst++) {
          if (jInst == iInst ) {
            Pinv[iInst][jInst] = 1.0 + Delta*D[iInst][jInst];
          }
          else {
            Pinv[iInst][jInst] = Delta*D[iInst][jInst];
          }
        }
      }

      /*--- Invert stabilization matrix Pinv with Gauss elimination---*/

      /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
      su2double **temp = new su2double*[nInstHB];
      for (i = 0; i < nInstHB; i++) {
        temp[i] = new su2double[2 * nInstHB];
      }

      /*---  Copy the desired matrix into the temporary matrix ---*/
      for (i = 0; i < nInstHB; i++) {
        for (j = 0; j < nInstHB; j++) {
          temp[i][j] = Pinv[i][j];
          temp[i][nInstHB + j] = 0;
        }
        temp[i][nInstHB + i] = 1;
      }

      su2double max_val;
      unsigned short max_idx;

      /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
      for (k = 0; k < nInstHB - 1; k++) {
        max_idx = k;
        max_val = abs(temp[k][k]);
        /*---  Find the largest value (pivot) in the column  ---*/
        for (j = k; j < nInstHB; j++) {
          if (abs(temp[j][k]) > max_val) {
            max_idx = j;
            max_val = abs(temp[j][k]);
          }
        }

        /*---  Move the row with the highest value up  ---*/
        for (j = 0; j < (nInstHB * 2); j++) {
          su2double d = temp[k][j];
          temp[k][j] = temp[max_idx][j];
          temp[max_idx][j] = d;
        }
        /*---  Subtract the moved row from all other rows ---*/
        for (i = k + 1; i < nInstHB; i++) {
          su2double c = temp[i][k] / temp[k][k];
          for (j = 0; j < (nInstHB * 2); j++) {
            temp[i][j] = temp[i][j] - temp[k][j] * c;
          }
        }
      }

      /*---  Back-substitution  ---*/
      for (k = nInstHB - 1; k > 0; k--) {
        if (temp[k][k] != su2double(0.0)) {
          for (int i = k - 1; i > -1; i--) {
            su2double c = temp[i][k] / temp[k][k];
            for (j = 0; j < (nInstHB * 2); j++) {
              temp[i][j] = temp[i][j] - temp[k][j] * c;
            }
          }
        }
      }

      /*---  Normalize the inverse  ---*/
      for (i = 0; i < nInstHB; i++) {
        su2double c = temp[i][i];
        for (j = 0; j < nInstHB; j++) {
          temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
        }
      }

      /*---  Copy the inverse back to the main program flow ---*/
      for (i = 0; i < nInstHB; i++) {
        for (j = 0; j < nInstHB; j++) {
          P[i][j] = temp[i][j + nInstHB];
        }
      }

      /*---  Delete dynamic template  ---*/
      for (iInst = 0; iInst < nInstHB; iInst++) {
        delete[] temp[iInst];
      }
      delete[] temp;

      /*--- Loop through variables to precondition ---*/
      for (iVar = 0; iVar < nVar; iVar++) {

        /*--- Get current source terms (not yet preconditioned) and zero source array to prepare preconditioning ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          Source_old[iInst] = solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->GetHarmonicBalance_Source(iPoint, iVar);
          Source[iInst] = 0;
        }

        /*--- Step through columns ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          for (jInst = 0; jInst < nInstHB; jInst++) {
            Source[iInst] += P[iInst][jInst]*Source_old[jInst];
          }

          /*--- Store updated source terms for current node ---*/
          if (!adjoint) {
            solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
          }
          else {
            solver_container[ZONE_0][iInst][iMGlevel][ADJFLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
          }
        }

      }
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iInst = 0; iInst < nInstHB; iInst++){
    delete [] P[iInst];
    delete [] Pinv[iInst];
  }
  delete [] P;
  delete [] Pinv;
  delete [] Source;
  delete [] Source_old;

}

void CHBDriver::ComputeHB_Operator() {

  const   complex<su2double> J(0.0,1.0);
  unsigned short i, j, k, iInst;

  su2double *Omega_HB       = new su2double[nInstHB];
  complex<su2double> **E    = new complex<su2double>*[nInstHB];
  complex<su2double> **Einv = new complex<su2double>*[nInstHB];
  complex<su2double> **DD   = new complex<su2double>*[nInstHB];
  for (iInst = 0; iInst < nInstHB; iInst++) {
    E[iInst]    = new complex<su2double>[nInstHB];
    Einv[iInst] = new complex<su2double>[nInstHB];
    DD[iInst]   = new complex<su2double>[nInstHB];
  }

  /*--- Get simualation period from config file ---*/
  su2double Period = config_container[ZONE_0]->GetHarmonicBalance_Period();

  /*--- Non-dimensionalize the input period, if necessary.      */
  Period /= config_container[ZONE_0]->GetTime_Ref();

  /*--- Build the array containing the selected frequencies to solve ---*/
  for (iInst = 0; iInst < nInstHB; iInst++) {
    Omega_HB[iInst]  = config_container[ZONE_0]->GetOmega_HB()[iInst];
    Omega_HB[iInst] /= config_container[ZONE_0]->GetOmega_Ref(); //TODO: check
  }

  /*--- Build the diagonal matrix of the frequencies DD ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      if (k == i ) {
        DD[i][k] = J*Omega_HB[k];
      }
    }
  }


  /*--- Build the harmonic balance inverse matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      Einv[i][k] = complex<su2double>(cos(Omega_HB[k]*(i*Period/nInstHB))) + J*complex<su2double>(sin(Omega_HB[k]*(i*Period/nInstHB)));
    }
  }

  /*---  Invert inverse harmonic balance Einv with Gauss elimination ---*/

  /*--  A temporary matrix to hold the inverse, dynamically allocated ---*/
  complex<su2double> **temp = new complex<su2double>*[nInstHB];
  for (i = 0; i < nInstHB; i++) {
    temp[i] = new complex<su2double>[2 * nInstHB];
  }

  /*---  Copy the desired matrix into the temporary matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (j = 0; j < nInstHB; j++) {
      temp[i][j] = Einv[i][j];
      temp[i][nInstHB + j] = 0;
    }
    temp[i][nInstHB + i] = 1;
  }

  su2double max_val;
  unsigned short max_idx;

  /*---  Pivot each column such that the largest number possible divides the other rows  ---*/
  for (k = 0; k < nInstHB - 1; k++) {
    max_idx = k;
    max_val = abs(temp[k][k]);
    /*---  Find the largest value (pivot) in the column  ---*/
    for (j = k; j < nInstHB; j++) {
      if (abs(temp[j][k]) > max_val) {
        max_idx = j;
        max_val = abs(temp[j][k]);
      }
    }
    /*---  Move the row with the highest value up  ---*/
    for (j = 0; j < (nInstHB * 2); j++) {
      complex<su2double> d = temp[k][j];
      temp[k][j] = temp[max_idx][j];
      temp[max_idx][j] = d;
    }
    /*---  Subtract the moved row from all other rows ---*/
    for (i = k + 1; i < nInstHB; i++) {
      complex<su2double> c = temp[i][k] / temp[k][k];
      for (j = 0; j < (nInstHB * 2); j++) {
        temp[i][j] = temp[i][j] - temp[k][j] * c;
      }
    }
  }
  /*---  Back-substitution  ---*/
  for (k = nInstHB - 1; k > 0; k--) {
    if (temp[k][k] != complex<su2double>(0.0)) {
      for (int i = k - 1; i > -1; i--) {
        complex<su2double> c = temp[i][k] / temp[k][k];
        for (j = 0; j < (nInstHB * 2); j++) {
          temp[i][j] = temp[i][j] - temp[k][j] * c;
        }
      }
    }
  }
  /*---  Normalize the inverse  ---*/
  for (i = 0; i < nInstHB; i++) {
    complex<su2double> c = temp[i][i];
    for (j = 0; j < nInstHB; j++) {
      temp[i][j + nInstHB] = temp[i][j + nInstHB] / c;
    }
  }
  /*---  Copy the inverse back to the main program flow ---*/
  for (i = 0; i < nInstHB; i++) {
    for (j = 0; j < nInstHB; j++) {
      E[i][j] = temp[i][j + nInstHB];
    }
  }
  /*---  Delete dynamic template  ---*/
  for (i = 0; i < nInstHB; i++) {
    delete[] temp[i];
  }
  delete[] temp;


  /*---  Temporary matrix for performing product  ---*/
  complex<su2double> **Temp    = new complex<su2double>*[nInstHB];

  /*---  Temporary complex HB operator  ---*/
  complex<su2double> **Dcpx    = new complex<su2double>*[nInstHB];

  for (iInst = 0; iInst < nInstHB; iInst++){
    Temp[iInst]    = new complex<su2double>[nInstHB];
    Dcpx[iInst]   = new complex<su2double>[nInstHB];
  }


  /*---  Calculation of the HB operator matrix ---*/
  for (int row = 0; row < nInstHB; row++) {
    for (int col = 0; col < nInstHB; col++) {
      for (int inner = 0; inner < nInstHB; inner++) {
        Temp[row][col] += Einv[row][inner] * DD[inner][col];
      }
    }
  }

  unsigned short row, col, inner;

  for (row = 0; row < nInstHB; row++) {
    for (col = 0; col < nInstHB; col++) {
      for (inner = 0; inner < nInstHB; inner++) {
        Dcpx[row][col] += Temp[row][inner] * E[inner][col];
      }
    }
  }

  /*---  Take just the real part of the HB operator matrix ---*/
  for (i = 0; i < nInstHB; i++) {
    for (k = 0; k < nInstHB; k++) {
      D[i][k] = real(Dcpx[i][k]);
    }
  }

  /*--- Deallocate dynamic memory ---*/
  for (iInst = 0; iInst < nInstHB; iInst++){
    delete [] E[iInst];
    delete [] Einv[iInst];
    delete [] DD[iInst];
    delete [] Temp[iInst];
    delete [] Dcpx[iInst];
  }
  delete [] E;
  delete [] Einv;
  delete [] DD;
  delete [] Temp;
  delete [] Dcpx;
  delete [] Omega_HB;

}

CDiscAdjFSIDriver::CDiscAdjFSIDriver(char* confFile,
                                     unsigned short val_nZone,
                                     SU2_Comm MPICommunicator) : CDriver(confFile,
                                                                         val_nZone,
                                                                         MPICommunicator,
                                                                         false) {

  unsigned short iVar;
  unsigned short nVar_Flow = 0, nVar_Struct = 0;
  RecordingState = 0;
  CurrentRecording = 0;

  switch (config_container[ZONE_0]->GetKind_ObjFunc()){
  case DRAG_COEFFICIENT:
  case LIFT_COEFFICIENT:
  case SIDEFORCE_COEFFICIENT:
  case EFFICIENCY:
  case MOMENT_X_COEFFICIENT:
  case MOMENT_Y_COEFFICIENT:
  case MOMENT_Z_COEFFICIENT:
  case EQUIVALENT_AREA:
    Kind_Objective_Function = FLOW_OBJECTIVE_FUNCTION;
    break;
  case REFERENCE_GEOMETRY:
  case REFERENCE_NODE:
  case VOLUME_FRACTION:
  case TOPOL_DISCRETENESS:
  case TOPOL_COMPLIANCE:
    Kind_Objective_Function = FEM_OBJECTIVE_FUNCTION;
    break;
  default:
    Kind_Objective_Function = NO_OBJECTIVE_FUNCTION;
    break;
  }

  direct_iteration = new CIteration*[nZone];

  unsigned short iZone;
  for (iZone = 0; iZone < nZone; iZone++){
    switch (config_container[iZone]->GetKind_Solver()) {
      case DISC_ADJ_INC_RANS: case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES:
       case DISC_ADJ_RANS: case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES:
         direct_iteration[iZone] = new CFluidIteration(config_container[iZone]);
         nVar_Flow = solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL]->GetnVar();
         flow_criteria = -8.0;  // Temporarily hard coded until adapted into new structure
         flow_criteria_rel = 3.0;
         break;
       case DISC_ADJ_FEM:
         direct_iteration[iZone] = new CFEAIteration(config_container[iZone]);
         nVar_Struct = solver_container[iZone][INST_0][MESH_0][ADJFEA_SOL]->GetnVar();
         structure_criteria    = -8.0;
         structure_criteria_rel = 3.0;
         break;
    }
  }

  init_res_flow   = new su2double[nVar_Flow];
  init_res_struct = new su2double[nVar_Struct];

  residual_flow   = new su2double[nVar_Flow];
  residual_struct = new su2double[nVar_Struct];

  residual_flow_rel   = new su2double[nVar_Flow];
  residual_struct_rel = new su2double[nVar_Struct];

  for (iVar = 0; iVar < nVar_Flow; iVar++){
    init_res_flow[iVar] = 0.0;
    residual_flow[iVar] = 0.0;
    residual_flow_rel[iVar] = 0.0;
  }
  for (iVar = 0; iVar < nVar_Struct; iVar++){
    init_res_struct[iVar] = 0.0;
    residual_struct[iVar] = 0.0;
    residual_struct_rel[iVar] = 0.0;
  }


  bool write_history = true;

  /*--- Header of the temporary output file ---*/
  if ((write_history) && (rank == MASTER_NODE)){
    ofstream myfile_res;
    myfile_res.open ("history_adjoint_FSI.csv");

    myfile_res << "BGS_Iter\t";

    for (iVar = 0; iVar < nVar_Flow; iVar++){
      myfile_res << "ResFlow[" << iVar << "]\t";
    }

    for (iVar = 0; iVar < nVar_Struct; iVar++){
      myfile_res << "ResFEA[" << iVar << "]\t";
    }


    bool de_effects = config_container[ZONE_0]->GetDE_Effects();
    for (iVar = 0; iVar < config_container[ZONE_0]->GetnElasticityMod(); iVar++)
        myfile_res << "Sens_E_" << iVar << "\t";

    for (iVar = 0; iVar < config_container[ZONE_0]->GetnPoissonRatio(); iVar++)
      myfile_res << "Sens_Nu_" << iVar << "\t";

    if (de_effects){
        for (iVar = 0; iVar < config_container[ZONE_0]->GetnElectric_Field(); iVar++)
          myfile_res << "Sens_EField_" << iVar << "\t";
    }

    myfile_res << endl;

    myfile_res.close();
  }

  // TEST: for implementation of python framework in standalone structural problems
  if ((config_container[ZONE_1]->GetDV_FEA() != NODV_FEA) && (rank == MASTER_NODE)){

    /*--- Header of the temporary output file ---*/
    ofstream myfile_res;

    switch (config_container[ZONE_1]->GetDV_FEA()) {
      case YOUNG_MODULUS:
        myfile_res.open("grad_young.opt");
        break;
      case POISSON_RATIO:
        myfile_res.open("grad_poisson.opt");
        break;
      case DENSITY_VAL:
      case DEAD_WEIGHT:
        myfile_res.open("grad_density.opt");
        break;
      case ELECTRIC_FIELD:
        myfile_res.open("grad_efield.opt");
        break;
      default:
        myfile_res.open("grad.opt");
        break;
    }

    unsigned short iDV;
    unsigned short nDV = solver_container[ZONE_1][INST_0][MESH_0][ADJFEA_SOL]->GetnDVFEA();

    myfile_res << "INDEX" << "\t" << "GRAD" << endl;

    myfile_res.precision(15);

    for (iDV = 0; iDV < nDV; iDV++){
      myfile_res << iDV;
      myfile_res << "\t";
      myfile_res << scientific << solver_container[ZONE_1][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_DVFEA(iDV);
      myfile_res << endl;
    }

    myfile_res.close();
  }
  output_legacy = new COutputLegacy(config_container[ZONE_0]);

  /*--- TODO: This is a workaround until the TestCases.py script incorporates new classes for nested loops. ---*/
//  config_container[ZONE_0]->SetnExtIter(1);
//  config_container[ZONE_1]->SetnExtIter(1);
  ConvHist_file = NULL;
  ConvHist_file = new ofstream*[nZone];
  for (iZone = 0; iZone < nZone; iZone++) {
    ConvHist_file[iZone] = NULL;
    if (rank == MASTER_NODE){
      ConvHist_file[iZone] = new ofstream[nInst[iZone]];
      for (iInst = 0; iInst < nInst[iZone]; iInst++) {
        output_legacy->SetConvHistory_Header(&ConvHist_file[iZone][iInst], config_container[iZone], iZone, iInst);
        config_container[iZone]->SetHistFile(&ConvHist_file[iZone][INST_0]);
      }
    }
  }
}

CDiscAdjFSIDriver::~CDiscAdjFSIDriver(void) {

  delete [] direct_iteration;
  delete [] init_res_flow;
  delete [] init_res_struct;
  delete [] residual_flow;
  delete [] residual_struct;
  delete [] residual_flow_rel;
  delete [] residual_struct_rel;

}

void CDiscAdjFSIDriver::DynamicMeshUpdate(unsigned long ExtIter){

}

void CDiscAdjFSIDriver::Run( ) {

  /*--- As of now, we are coding it for just 2 zones. ---*/
  /*--- This will become more general, but we need to modify the configuration for that ---*/
  unsigned short ZONE_FLOW = 0, ZONE_STRUCT = 1;
  unsigned short iZone;
  bool BGS_Converged = false;

  unsigned long IntIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetInnerIter(IntIter);
  unsigned long iOuterIter = 0; for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(iOuterIter);
  unsigned long nOuterIter = driver_config->GetnOuter_Iter();

  ofstream myfile_struc, myfile_flow, myfile_geo;

  Preprocess(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);

  for (iOuterIter = 0; iOuterIter < nOuterIter && !BGS_Converged; iOuterIter++){

    if (rank == MASTER_NODE){
      cout << endl << "                    ****** BGS ITERATION ";
      cout << iOuterIter;
      cout << " ******" << endl;
    }

    for (iZone = 0; iZone < nZone; iZone++) config_container[iZone]->SetOuterIter(iOuterIter);

    /*--- Start with structural terms if OF is based on displacements ---*/

    if (Kind_Objective_Function == FEM_OBJECTIVE_FUNCTION)
      Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEA_DISP_VARS);

    /*--- Iterate fluid (including cross term) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FLOW_CONS_VARS);

    /*--- Compute mesh (it is a cross term dF / dMv ) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, MESH_COORDS);

    /*--- Compute mesh cross term (dM / dSv) ---*/

    Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEM_CROSS_TERM_GEOMETRY);

    /*--- End with structural terms if OF is based on fluid variables ---*/

    if (Kind_Objective_Function == FLOW_OBJECTIVE_FUNCTION)
      Iterate_Block(ZONE_FLOW, ZONE_STRUCT, FEA_DISP_VARS);

    /*--- Check convergence of the BGS method ---*/
    BGS_Converged = BGSConvergence(iOuterIter, ZONE_FLOW, ZONE_STRUCT);

  }

  output_container[ZONE_FLOW]->SetResult_Files(geometry_container[ZONE_FLOW][INST_0][MESH_0],
    config_container[ZONE_FLOW], solver_container[ZONE_FLOW][INST_0][MESH_0], 0, true);

  output_container[ZONE_STRUCT]->SetResult_Files(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
    config_container[ZONE_STRUCT], solver_container[ZONE_STRUCT][INST_0][MESH_0], 0, true);

}


void CDiscAdjFSIDriver::Preprocess(unsigned short ZONE_FLOW,
                  unsigned short ZONE_STRUCT,
                  unsigned short kind_recording){

  unsigned long IntIter = 0, iPoint;
  config_container[ZONE_0]->SetInnerIter(IntIter);
  unsigned short ExtIter = config_container[ZONE_FLOW]->GetTimeIter();

  bool dual_time_1st = (config_container[ZONE_FLOW]->GetTime_Marching() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[ZONE_FLOW]->GetTime_Marching() == DT_STEPPING_2ND);
  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  unsigned short iMesh;
  int Direct_Iter_Flow;
  bool update_geo = false;

  /*----------------------------------------------------------------------------*/
  /*------------------------------ FLOW SOLUTION -------------------------------*/
  /*----------------------------------------------------------------------------*/

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config_container[ZONE_FLOW]->GetTime_Marching()) {

    Direct_Iter_Flow = SU2_TYPE::Int(config_container[ZONE_FLOW]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 2;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (dual_time) {
      Direct_Iter_Flow += 1;
    }

    if (ExtIter == 0){

      if (dual_time_2nd) {

        /*--- Load solution at timestep n-2 ---*/

        iteration_container[ZONE_FLOW][INST_0]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, INST_0, Direct_Iter_Flow-2);

        /*--- Push solution back to correct array ---*/

        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
          solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1();
          if (turbulent) {
            solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
            solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1();
          }
        }
      }
      if (dual_time) {

        /*--- Load solution at timestep n-1 ---*/

        iteration_container[ZONE_FLOW][INST_0]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, INST_0, Direct_Iter_Flow-1);

        /*--- Push solution back to correct array ---*/

        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n();
          if (turbulent) {
            solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n();
          }
        }
      }

      /*--- Load solution timestep n ---*/

      iteration_container[ZONE_FLOW][INST_0]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, INST_0, Direct_Iter_Flow);

    }


    if ((ExtIter > 0) && dual_time){

      /*--- Load solution timestep n - 2 ---*/

      iteration_container[ZONE_FLOW][INST_0]->LoadUnsteady_Solution(geometry_container, solver_container,config_container, ZONE_FLOW, INST_0, Direct_Iter_Flow - 2);

      /*--- Temporarily store the loaded solution in the Solution_Old array ---*/

      for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
        solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_OldSolution();
        if (turbulent){
         solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_OldSolution();
        }
      }

      /*--- Set Solution at timestep n to solution at n-1 ---*/

      for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
        for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][INST_0][iMesh]->GetnPoint();iPoint++) {
          solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->SetSolution(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          if (turbulent) {
            solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->SetSolution(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->GetSolution_time_n(iPoint));
          }
        }
      }
      if (dual_time_1st){
      /*--- Set Solution at timestep n-1 to the previously loaded solution ---*/
        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][INST_0][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            if (turbulent) {
              solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
          }
        }
      }
      if (dual_time_2nd){
        /*--- Set Solution at timestep n-1 to solution at n-2 ---*/
        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][INST_0][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            if (turbulent) {
              solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->GetSolution_time_n1(iPoint));
            }
          }
        }
        /*--- Set Solution at timestep n-2 to the previously loaded solution ---*/
        for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
          for(iPoint=0; iPoint<geometry_container[ZONE_FLOW][INST_0][iMesh]->GetnPoint();iPoint++) {
            solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->Set_Solution_time_n1(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint));
            if (turbulent) {
              solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->Set_Solution_time_n1(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][TURB_SOL]->GetNodes()->GetSolution_Old(iPoint));
            }
          }
        }
      }
    }
  }
  else{

    /*--- Load the restart (we need to use the routine in order to get the GEOMETRY, otherwise it's restarted from the base mesh ---*/

    solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[ZONE_FLOW][INST_0], solver_container[ZONE_FLOW][INST_0], config_container[ZONE_FLOW], 0, true);

  if (ExtIter == 0 || dual_time) {

    for (iMesh=0; iMesh<=config_container[ZONE_FLOW]->GetnMGLevels();iMesh++) {
      for (iPoint = 0; iPoint < geometry_container[ZONE_FLOW][INST_0][iMesh]->GetnPoint(); iPoint++) {
        solver_container[ZONE_FLOW][INST_0][iMesh][ADJFLOW_SOL]->GetNodes()->SetSolution_Direct(iPoint, solver_container[ZONE_FLOW][INST_0][iMesh][FLOW_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
    if (turbulent && !config_container[ZONE_FLOW]->GetFrozen_Visc_Disc()) {
      for (iPoint = 0; iPoint < geometry_container[ZONE_FLOW][INST_0][MESH_0]->GetnPoint(); iPoint++) {
        solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->GetNodes()->SetSolution_Direct(iPoint, solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint));
      }
    }
  }

    /*--- Store geometry of the converged solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[ZONE_FLOW][INST_0][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetNodes()->SetGeometry_Direct(iPoint, geometry_container[ZONE_FLOW][INST_0][MESH_0]->node[iPoint]->GetCoord());
    }

  }

  /*----------------------------------------------------------------------------*/
  /*-------------------------- STRUCTURAL SOLUTION -----------------------------*/
  /*----------------------------------------------------------------------------*/

  IntIter = 0;
  config_container[ZONE_STRUCT]->SetInnerIter(IntIter);
  ExtIter = config_container[ZONE_STRUCT]->GetTimeIter();
  bool dynamic = (config_container[ZONE_STRUCT]->GetTime_Domain());

  int Direct_Iter_FEA;

  /*--- For the dynamic adjoint, load direct solutions from restart files. ---*/

  if (dynamic) {

    Direct_Iter_FEA = SU2_TYPE::Int(config_container[ZONE_STRUCT]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;

    /*--- We want to load the already converged solution at timesteps n and n-1 ---*/

    /*--- Load solution at timestep n-1 ---*/

    iteration_container[ZONE_STRUCT][INST_0]->LoadDynamic_Solution(geometry_container, solver_container,config_container, ZONE_STRUCT, INST_0, Direct_Iter_FEA-1);

    /*--- Push solution back to correct array ---*/

    solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->Set_Solution_time_n();

    /*--- Push solution back to correct array ---*/

    solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->SetSolution_Accel_time_n();

    /*--- Push solution back to correct array ---*/

    solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->SetSolution_Vel_time_n();

    /*--- Load solution timestep n ---*/

    iteration_container[ZONE_STRUCT][INST_0]->LoadDynamic_Solution(geometry_container, solver_container,config_container, ZONE_STRUCT, INST_0, Direct_Iter_FEA);

    /*--- Store FEA solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][INST_0][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Direct(iPoint, solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->GetSolution(iPoint));
    }

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][INST_0][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Accel_Direct(iPoint, solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->GetSolution_Accel(iPoint));
    }

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][INST_0][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Vel_Direct(iPoint, solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->GetSolution_Vel(iPoint));
    }

  }
  else {

    solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->LoadRestart(geometry_container[ZONE_STRUCT][INST_0], solver_container[ZONE_STRUCT][INST_0], config_container[ZONE_STRUCT], 0, update_geo);

    /*--- Store FEA solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[ZONE_STRUCT][INST_0][MESH_0]->GetnPoint(); iPoint++){
      solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetNodes()->SetSolution_Direct(iPoint, solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetNodes()->GetSolution(iPoint));
    }

  }

  /*----------------------------------------------------------------------------*/
  /*--------------------- ADJOINT SOLVER PREPROCESSING -------------------------*/
  /*----------------------------------------------------------------------------*/

  solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0], solver_container[ZONE_FLOW][INST_0][MESH_0],  config_container[ZONE_FLOW] , MESH_0, 0, RUNTIME_ADJFLOW_SYS, false);

  if (turbulent){
    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->Preprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0], solver_container[ZONE_FLOW][INST_0][MESH_0],  config_container[ZONE_FLOW] , MESH_0, 0, RUNTIME_ADJTURB_SYS, false);
  }

  solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->Preprocessing(geometry_container[ZONE_STRUCT][INST_0][MESH_0], solver_container[ZONE_STRUCT][INST_0][MESH_0],  config_container[ZONE_STRUCT] , MESH_0, 0, RUNTIME_ADJFEA_SYS, false);



}

void CDiscAdjFSIDriver::PrintDirect_Residuals(unsigned short ZONE_FLOW,
                                                          unsigned short ZONE_STRUCT,
                                                          unsigned short kind_recording){

  unsigned short ExtIter = config_container[ZONE_FLOW]->GetTimeIter();
  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool nonlinear_analysis = (config_container[ZONE_STRUCT]->GetGeometricConditions() == LARGE_DEFORMATIONS);   // Nonlinear analysis.
  bool unsteady = config_container[ZONE_FLOW]->GetTime_Marching() != NONE;
  bool dynamic = (config_container[ZONE_STRUCT]->GetTime_Domain());

  su2double val_OFunction = 0.0;
  string kind_OFunction;

  cout.precision(6);
  cout.setf(ios::scientific, ios::floatfield);

  if ((kind_recording == FLOW_CONS_VARS) || (kind_recording == MESH_COORDS)) {

    /*--- Print residuals in the first iteration ---*/

    if (rank == MASTER_NODE && ((ExtIter == 0) || unsteady )){
      cout << "log10[RMS Density]: "<< log10(solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))
                     <<", Drag: " <<solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CD()
                     <<", Lift: " << solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CL() << "." << endl;

      if (turbulent){
        cout << "log10[RMS k]: " << log10(solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
      }
      if (Kind_Objective_Function == FLOW_OBJECTIVE_FUNCTION){
        switch (config_container[ZONE_FLOW]->GetKind_ObjFunc()){
        case DRAG_COEFFICIENT:
          kind_OFunction = "(Drag coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CD();
          break;
        case LIFT_COEFFICIENT:
          kind_OFunction = "(Lift coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CL();
          break;
        case SIDEFORCE_COEFFICIENT:
          kind_OFunction = "(Sideforce coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CSF();
          break;
        case EFFICIENCY:
          kind_OFunction = "(Efficiency): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CEff();
          break;
        case MOMENT_X_COEFFICIENT:
          kind_OFunction = "(Moment X coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CMx();
          break;
        case MOMENT_Y_COEFFICIENT:
          kind_OFunction = "(Moment Y coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CMy();
          break;
        case MOMENT_Z_COEFFICIENT:
          kind_OFunction = "(Moment Z coefficient): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CMz();
          break;
        case EQUIVALENT_AREA:
          kind_OFunction = "(Equivalent area): ";
          val_OFunction = solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->GetTotal_CEquivArea();
          break;
        default:
          val_OFunction = 0.0;  // If the objective function is computed in a different physical problem
          break;
        }
        cout << "Objective function " << kind_OFunction << val_OFunction << endl;
      }
    }

  }

  if ((kind_recording == FEA_DISP_VARS) || (kind_recording == FLOW_CROSS_TERM) || (kind_recording == GEOMETRY_CROSS_TERM)) {

    if (rank == MASTER_NODE && ((ExtIter == 0) || dynamic )){
      if (nonlinear_analysis){
        cout << "UTOL-A: "   << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(0))
             << ", RTOL-A: " << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(1))
             << ", ETOL-A: " << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_FEM(2)) << "." << endl;
      }
      else{
        if (solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetnVar() == 2){
          cout << "log10[RMS Ux]: "   << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(0))
               << ", log10[RMS Uy]: " << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(1)) << "." << endl;

        }
        else{
          cout << "log10[RMS Ux]: "   << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(0))
               << ", log10[RMS Uy]: " << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(1))
               << ", log10[RMS Uz]: " << log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetRes_RMS(2))<< "." << endl;
        }

      }
      if (Kind_Objective_Function == FEM_OBJECTIVE_FUNCTION){
        switch (config_container[ZONE_STRUCT]->GetKind_ObjFunc()){
        case REFERENCE_GEOMETRY:
          kind_OFunction = "(Reference Geometry): ";
          val_OFunction = solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetTotal_OFRefGeom();
          break;
        case REFERENCE_NODE:
          kind_OFunction = "(Reference Node): ";
          val_OFunction = solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetTotal_OFRefNode();
          break;
        case VOLUME_FRACTION:
          kind_OFunction = "(Volume Fraction): ";
          val_OFunction = solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetTotal_OFVolFrac();
          break;
        case TOPOL_DISCRETENESS:
          kind_OFunction = "(Topology discreteness): ";
          val_OFunction = solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetTotal_OFVolFrac();
          break;
        case TOPOL_COMPLIANCE:
          kind_OFunction = "(Topology compliance): ";
          val_OFunction = solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->GetTotal_OFCompliance();
          break;
        default:
          val_OFunction = 0.0;  // If the objective function is computed in a different physical problem
          break;
        }
        cout << "Objective function " << kind_OFunction << val_OFunction << endl;
      }
    }

  }

}

void CDiscAdjFSIDriver::Iterate_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT, unsigned short kind_recording){

  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {

    Fluid_Iteration_Direct(ZONE_FLOW, ZONE_STRUCT);


  }

  if ((kind_recording == FEA_DISP_VARS) ||
      (kind_recording == FLOW_CROSS_TERM) ||
      (kind_recording == GEOMETRY_CROSS_TERM)) {

    Structural_Iteration_Direct(ZONE_FLOW, ZONE_STRUCT);

  }


  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {

    Mesh_Deformation_Direct(ZONE_FLOW, ZONE_STRUCT);

  }


}

void CDiscAdjFSIDriver::Fluid_Iteration_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool frozen_visc = config_container[ZONE_FLOW]->GetFrozen_Visc_Disc();

  /*-----------------------------------------------------------------*/
  /*------------------- Set Dependency on Geometry ------------------*/
  /*-----------------------------------------------------------------*/

  geometry_container[ZONE_FLOW][INST_0][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW][INST_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->InitiateComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);
  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->CompleteComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0],solver_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  if (turbulent && !frozen_visc) {
    solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->Postprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0], solver_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], MESH_0);

    solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->InitiateComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION_EDDY);
    solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->CompleteComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION_EDDY);

  }

  /*-----------------------------------------------------------------*/
  /*----------------- Iterate the flow solver -----------------------*/
  /*---- Sets all the cross dependencies for the flow variables -----*/
  /*-----------------------------------------------------------------*/

  config_container[ZONE_FLOW]->SetInnerIter(0);

  direct_iteration[ZONE_FLOW]->Iterate(output_container[ZONE_FLOW], integration_container, geometry_container,
      solver_container, numerics_container, config_container,
      surface_movement, grid_movement, FFDBox, ZONE_FLOW, INST_0);

  /*-----------------------------------------------------------------*/
  /*--------------------- Set MPI Solution --------------------------*/
  /*-----------------------------------------------------------------*/

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->InitiateComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);
  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->CompleteComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);

}

void CDiscAdjFSIDriver::Structural_Iteration_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool frozen_visc = config_container[ZONE_FLOW]->GetFrozen_Visc_Disc();

  /*-----------------------------------------------------------------*/
  /*---------- Set Dependencies on Geometry and Flow ----------------*/
  /*-----------------------------------------------------------------*/

  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->InitiateComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);
  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->CompleteComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);

  geometry_container[ZONE_FLOW][INST_0][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW][INST_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->InitiateComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);
  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->CompleteComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0],solver_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  if (turbulent && !frozen_visc) {
    solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->Postprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0], solver_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], MESH_0);

    solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->InitiateComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION_EDDY);
    solver_container[ZONE_FLOW][INST_0][MESH_0][TURB_SOL]->CompleteComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION_EDDY);

  }

  /*-----------------------------------------------------------------*/
  /*-------------------- Transfer Tractions -------------------------*/
  /*-----------------------------------------------------------------*/

  Transfer_Tractions(ZONE_FLOW, ZONE_STRUCT);

  /*-----------------------------------------------------------------*/
  /*--------------- Iterate the structural solver -------------------*/
  /*-----------------------------------------------------------------*/

  direct_iteration[ZONE_STRUCT]->Iterate(output_container[ZONE_STRUCT], integration_container, geometry_container,
                                        solver_container, numerics_container, config_container,
                                        surface_movement, grid_movement, FFDBox, ZONE_STRUCT, INST_0);

  /*-----------------------------------------------------------------*/
  /*--------------------- Set MPI Solution --------------------------*/
  /*-----------------------------------------------------------------*/

  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->InitiateComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);
  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->CompleteComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);

}

void CDiscAdjFSIDriver::Mesh_Deformation_Direct(unsigned short ZONE_FLOW, unsigned short ZONE_STRUCT) {

  unsigned long ExtIter = config_container[ZONE_STRUCT]->GetTimeIter();

  /*-----------------------------------------------------------------*/
  /*--------------------- Set MPI Solution --------------------------*/
  /*-----------------------------------------------------------------*/

  geometry_container[ZONE_FLOW][INST_0][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW][INST_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->InitiateComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);
  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->CompleteComms(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], SOLUTION);

  solver_container[ZONE_FLOW][INST_0][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[ZONE_FLOW][INST_0][MESH_0],solver_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);

  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->InitiateComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);
  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->CompleteComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);

  /*-----------------------------------------------------------------*/
  /*------------------- Transfer Displacements ----------------------*/
  /*-----------------------------------------------------------------*/

  Transfer_Displacements(ZONE_STRUCT, ZONE_FLOW);

  /*-----------------------------------------------------------------*/
  /*------------------- Set the Grid movement -----------------------*/
  /*---- No longer done in the preprocess of the flow iteration -----*/
  /*---- as the flag Grid_Movement is set to false in this case -----*/
  /*-----------------------------------------------------------------*/

  direct_iteration[ZONE_FLOW]->SetGrid_Movement(geometry_container[ZONE_FLOW][INST_0],
                                                               surface_movement[ZONE_FLOW], grid_movement[ZONE_FLOW][INST_0],
                                                               solver_container[ZONE_FLOW][INST_0], config_container[ZONE_FLOW], 0, ExtIter );

  geometry_container[ZONE_FLOW][INST_0][MESH_0]->UpdateGeometry(geometry_container[ZONE_FLOW][INST_0], config_container[ZONE_FLOW]);

  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->InitiateComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);
  solver_container[ZONE_STRUCT][INST_0][MESH_0][FEA_SOL]->CompleteComms(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT], SOLUTION_FEA);

}

void CDiscAdjFSIDriver::SetRecording(unsigned short ZONE_FLOW,
                                              unsigned short ZONE_STRUCT,
                                              unsigned short kind_recording){

  unsigned long IntIter = config_container[ZONE_0]->GetInnerIter();
  bool unsteady = (config_container[ZONE_FLOW]->GetTime_Marching() != NONE);
  bool dynamic = (config_container[ZONE_STRUCT]->GetTime_Domain());

  string kind_DirectIteration = " ";
  string kind_AdjointIteration = " ";

  if (unsteady || dynamic){
    SU2_MPI::Error("DYNAMIC ADJOINT SOLVER NOT IMPLEMENTED FOR FSI APPLICATIONS", CURRENT_FUNCTION);
  }


  if (rank == MASTER_NODE){
    cout << endl;
    switch (kind_recording){
    case FLOW_CONS_VARS:
      kind_AdjointIteration = "Flow iteration: flow input -> flow output";
      kind_DirectIteration = "flow ";
      break;
    case MESH_COORDS:
      kind_AdjointIteration = "Geometry cross term from flow: geometry input -> flow output";
      kind_DirectIteration = "flow ";
      break;
    case FEA_DISP_VARS:
      kind_AdjointIteration = "Structural iteration: structural input -> structural output";
      kind_DirectIteration = "structural ";
      break;
    case FLOW_CROSS_TERM:
      kind_AdjointIteration = "Flow cross term: flow input -> structural output";
      kind_DirectIteration = "structural ";
      break;
    case GEOMETRY_CROSS_TERM:
      kind_AdjointIteration = "Geometry cross term from structure: geometry input -> structural output";
      kind_DirectIteration = "structural ";
      break;
    case FEM_CROSS_TERM_GEOMETRY:
      kind_AdjointIteration = "Structural cross term from geometry: structural input -> geometry output";
      kind_DirectIteration = "mesh deformation ";
      break;
    }
    cout << kind_AdjointIteration << endl;
    cout << "Direct " << kind_DirectIteration << "iteration to store computational graph." << endl;
    switch (kind_recording){
    case FLOW_CONS_VARS: case MESH_COORDS:
    case FEA_DISP_VARS: case FLOW_CROSS_TERM: case GEOMETRY_CROSS_TERM:
      cout << "Compute residuals to check the convergence of the direct problem." << endl; break;
    case FEM_CROSS_TERM_GEOMETRY:
      cout << "Deform the grid using the converged solution of the direct problem." << endl; break;
    }
  }


  AD::Reset();

  if (CurrentRecording != kind_recording && (CurrentRecording != NONE) ){

    /*--- Clear indices ---*/

    PrepareRecording(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);

    /*--- Clear indices of coupling variables ---*/

    SetDependencies(ZONE_FLOW, ZONE_STRUCT, ALL_VARIABLES);

    /*--- Run one iteration while tape is passive - this clears all indices ---*/
    Iterate_Direct(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  }

  /*--- Prepare for recording ---*/

  PrepareRecording(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Start the recording of all operations ---*/

  AD::StartRecording();

  /*--- Register input variables ---*/

  RegisterInput(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Set dependencies for flow, geometry and structural solvers ---*/

  SetDependencies(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Run a direct iteration ---*/
  Iterate_Direct(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Register objective function and output variables ---*/

  RegisterOutput(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Stop the recording ---*/
  AD::StopRecording();

  /*--- Set the recording status ---*/

  CurrentRecording = kind_recording;

  /* --- Reset the number of the internal iterations---*/

  config_container[ZONE_0]->SetInnerIter(IntIter);


}

void CDiscAdjFSIDriver::PrepareRecording(unsigned short ZONE_FLOW,
                                                   unsigned short ZONE_STRUCT,
                                                   unsigned short kind_recording){

  unsigned short iMesh;
  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);

  /*--- Set fluid variables to direct solver values ---*/
  for (iMesh = 0; iMesh <= config_container[ZONE_FLOW]->GetnMGLevels(); iMesh++){
    solver_container[ZONE_FLOW][INST_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW]);
  }
  if (turbulent){
    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[ZONE_FLOW][INST_0][MESH_0], config_container[ZONE_FLOW]);
  }

  /*--- Set geometry to the converged values ---*/

  solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->SetMesh_Recording(geometry_container[ZONE_FLOW][INST_0], grid_movement[ZONE_FLOW][INST_0], config_container[ZONE_FLOW]);

  /*--- Set structural variables to direct solver values ---*/

  solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->SetRecording(geometry_container[ZONE_STRUCT][INST_0][MESH_0], config_container[ZONE_STRUCT]);

}

void CDiscAdjFSIDriver::RegisterInput(unsigned short ZONE_FLOW,
                                               unsigned short ZONE_STRUCT,
                                               unsigned short kind_recording){

  /*--- Register flow variables ---*/
  if (kind_recording == FLOW_CONS_VARS) {
    iteration_container[ZONE_FLOW][INST_0]->RegisterInput(solver_container, geometry_container, config_container, ZONE_FLOW, INST_0, kind_recording);
  }

  /*--- Register geometry variables ---*/
  if (kind_recording == MESH_COORDS) {
    iteration_container[ZONE_FLOW][INST_0]->RegisterInput(solver_container, geometry_container, config_container, ZONE_FLOW, INST_0, kind_recording);
  }

  /*--- Register structural variables ---*/
  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {
    iteration_container[ZONE_STRUCT][INST_0]->RegisterInput(solver_container, geometry_container, config_container, ZONE_STRUCT, INST_0, kind_recording);
  }

  /*--- Register all variables ---*/
  if (kind_recording == FEA_DISP_VARS) {
    iteration_container[ZONE_STRUCT][INST_0]->RegisterInput(solver_container, geometry_container, config_container, ZONE_STRUCT, INST_0, FEA_DISP_VARS);
    iteration_container[ZONE_FLOW][INST_0]->RegisterInput(solver_container, geometry_container, config_container, ZONE_FLOW, INST_0, FLOW_CROSS_TERM);
    iteration_container[ZONE_FLOW][INST_0]->RegisterInput(solver_container, geometry_container, config_container, ZONE_FLOW, INST_0, GEOMETRY_CROSS_TERM);
  }

}

void CDiscAdjFSIDriver::SetDependencies(unsigned short ZONE_FLOW,
                                                  unsigned short ZONE_STRUCT,
                                                  unsigned short kind_recording){

  /*--- Add dependencies for geometrical and turbulent variables ---*/

  iteration_container[ZONE_FLOW][INST_0]->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_FLOW, INST_0, kind_recording);

  /*--- Add dependencies for E, Nu, Rho, and Rho_DL variables ---*/

  iteration_container[ZONE_STRUCT][INST_0]->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_STRUCT, INST_0, kind_recording);


}

void CDiscAdjFSIDriver::RegisterOutput(unsigned short ZONE_FLOW,
                                                 unsigned short ZONE_STRUCT,
                                                 unsigned short kind_recording){

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool frozen_visc = config_container[ZONE_FLOW]->GetFrozen_Visc_Disc();


  /*--- Register a flow-type objective function and the conservative variables of the flow as output of the iteration. ---*/
  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {
    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[ZONE_FLOW]);

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[ZONE_FLOW][INST_0][MESH_0],config_container[ZONE_FLOW]);

    if (turbulent && !frozen_visc) {
      solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->RegisterOutput(geometry_container[ZONE_FLOW][INST_0][MESH_0],config_container[ZONE_FLOW]);
    }
  }


  /*--- Register a structural-type objective function and the displacements of the structure as output of the iteration. ---*/
  if (kind_recording == FEA_DISP_VARS) {
    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->RegisterObj_Func(config_container[ZONE_STRUCT]);

    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->RegisterOutput(geometry_container[ZONE_STRUCT][INST_0][MESH_0],config_container[ZONE_STRUCT]);
  }


  /*--- The FEM_CROSS_TERM_GEOMETRY evaluates the mesh routines, they do not throw any dependency on the objective function. ---*/
  /*--- Register the displacements of the fluid nodes as output of the iteration. ---*/
  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {
    geometry_container[ZONE_FLOW][INST_0][MESH_0]->RegisterOutput_Coordinates(config_container[ZONE_FLOW]);
  }

}


void CDiscAdjFSIDriver::Iterate_Block(unsigned short ZONE_FLOW,
                                                unsigned short ZONE_STRUCT,
                                                unsigned short kind_recording){

  unsigned long IntIter=0, nIntIter = 1;
  bool dual_time_1st = (config_container[ZONE_0]->GetTime_Marching() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[ZONE_0]->GetTime_Marching() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  bool dynamic = (config_container[ZONE_STRUCT]->GetTime_Domain());

  bool adjoint_convergence = false;

  /*--- Record one direct iteration with kind_recording as input ---*/

  SetRecording(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Print the residuals of the direct subiteration ---*/

  PrintDirect_Residuals(ZONE_FLOW, ZONE_STRUCT, kind_recording);

  /*--- Run the iteration ---*/

  switch (kind_recording){
  case FLOW_CONS_VARS:
    nIntIter = config_container[ZONE_FLOW]->GetnInner_Iter();
    break;
  case FEA_DISP_VARS:
    nIntIter = config_container[ZONE_STRUCT]->GetnInner_Iter();
    break;
  case MESH_COORDS:
  case FEM_CROSS_TERM_GEOMETRY:
  case FLOW_CROSS_TERM:
  case GEOMETRY_CROSS_TERM:
    nIntIter = 1;
    break;
  }

  for (unsigned short iZone = 0; iZone < config_container[ZONE_FLOW]->GetnZone(); iZone++)
    config_container[iZone]->SetInnerIter(IntIter);

  for(IntIter = 0; IntIter < nIntIter; IntIter++){

    /*--- Set the internal iteration ---*/

    for (unsigned short iZone = 0; iZone < config_container[ZONE_FLOW]->GetnZone(); iZone++)
      config_container[iZone]->SetInnerIter(IntIter);

    /*--- Set the adjoint values of the flow and objective function ---*/

    InitializeAdjoint(ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Run the adjoint computation ---*/

    AD::ComputeAdjoint();

    /*--- Extract the adjoints of the input variables and store them for the next iteration ---*/

    ExtractAdjoint(ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/
    AD::ClearAdjoints();

    /*--- Check the convergence of the adjoint block ---*/

    adjoint_convergence = CheckConvergence(IntIter, ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Write the convergence history (only screen output) ---*/

    ConvergenceHistory(IntIter, nIntIter, ZONE_FLOW, ZONE_STRUCT, kind_recording);

    /*--- Break the loop if converged ---*/

    if (adjoint_convergence) break;


  }

  if (dual_time){
    integration_container[ZONE_FLOW][INST_0][ADJFLOW_SOL]->SetConvergence(false);
  }
  if (dynamic){
    integration_container[ZONE_FLOW][INST_0][ADJFLOW_SOL]->SetConvergence(false);
  }

}


void CDiscAdjFSIDriver::InitializeAdjoint(unsigned short ZONE_FLOW,
                                                     unsigned short ZONE_STRUCT,
                                                     unsigned short kind_recording){

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool frozen_visc = config_container[ZONE_FLOW]->GetFrozen_Visc_Disc();

  /*--- Seed a fluid-type objective function and initialize the adjoints of fluid conservative variables. ---*/
  if ((kind_recording == FLOW_CONS_VARS) ||
      (kind_recording == MESH_COORDS)) {
    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                             config_container[ZONE_FLOW]);

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                                config_container[ZONE_FLOW]);
    if (turbulent && !frozen_visc) {
      solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->SetAdjoint_Output(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                                  config_container[ZONE_FLOW]);
    }
  }

  /*--- Seed a structural-type objective function and initialize the adjoints of structural displacements. ---*/
  if (kind_recording == FEA_DISP_VARS) {
    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->SetAdj_ObjFunc(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                              config_container[ZONE_STRUCT]);

    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->SetAdjoint_Output(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                                 config_container[ZONE_STRUCT]);
  }

  /*--- Initialize the adjoints of fluid grid nodes. ---*/
  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {
    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->SetAdjoint_OutputMesh(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                                    config_container[ZONE_FLOW]);
  }
}

void CDiscAdjFSIDriver::ExtractAdjoint(unsigned short ZONE_FLOW,
                                                  unsigned short ZONE_STRUCT,
                                                  unsigned short kind_recording){

  bool turbulent = (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_RANS) ||
                   (config_container[ZONE_FLOW]->GetKind_Solver() == DISC_ADJ_INC_RANS);
  bool frozen_visc = config_container[ZONE_FLOW]->GetFrozen_Visc_Disc();

  /*--- Extract the adjoint of the fluid conservative variables ---*/

  if (kind_recording == FLOW_CONS_VARS) {

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                      config_container[ZONE_FLOW]);

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                      config_container[ZONE_FLOW]);

    if (turbulent && !frozen_visc) {
      solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                        config_container[ZONE_FLOW]);
    }

  }

  /*--- Extract the adjoint of the mesh coordinates ---*/

  if (kind_recording == MESH_COORDS) {

    /*--- Extract the adjoints of the flow geometry and store them for the next iteration ---*/

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_CrossTerm_Geometry_Flow(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                      config_container[ZONE_FLOW]);

  }

  /*--- Extract the adjoint of the structural displacements ---*/

  if (kind_recording == FEA_DISP_VARS) {

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Solution(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                               config_container[ZONE_STRUCT]);

    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Variables(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                                config_container[ZONE_STRUCT]);

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_CrossTerm(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                      config_container[ZONE_FLOW]);

    if (turbulent && !frozen_visc)
      solver_container[ZONE_FLOW][INST_0][MESH_0][ADJTURB_SOL]->ExtractAdjoint_CrossTerm(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                                         config_container[ZONE_FLOW]);

    solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_CrossTerm_Geometry(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                                                config_container[ZONE_FLOW]);
  }


  if (kind_recording == FEM_CROSS_TERM_GEOMETRY) {

    /*--- Extract the adjoints of the displacements (input variables) and store them for the next iteration ---*/

    solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->ExtractAdjoint_CrossTerm_Geometry(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                                config_container[ZONE_STRUCT]);
  }

}


bool CDiscAdjFSIDriver::CheckConvergence(unsigned long IntIter,
                                                   unsigned short ZONE_FLOW,
                                                   unsigned short ZONE_STRUCT,
                                                   unsigned short kind_recording){

  bool flow_convergence    = false,
       struct_convergence  = false,
       adjoint_convergence = false;

//  su2double residual_1, residual_2;

//  if (kind_recording == FLOW_CONS_VARS) {

//      /*--- Set the convergence criteria (only residual possible as of now) ---*/

//      residual_1 = log10(solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0));
//      residual_2 = log10(solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(1));

//      flow_convergence = ((residual_1 < config_container[ZONE_FLOW]->GetMinLogResidual()) &&
//                          (residual_2 < config_container[ZONE_FLOW]->GetMinLogResidual()));

//  }

//  if (kind_recording == FEA_DISP_VARS) {

//    /*--- Set the convergence criteria (only residual possible as of now) ---*/

//    residual_1 = log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetRes_RMS(0));
//    residual_2 = log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetRes_RMS(1));

//    // Temporary, until function is added
//    struct_convergence = ((residual_1 < config_container[ZONE_STRUCT]->GetResidual_FEM_UTOL()) &&
//                          (residual_2 < config_container[ZONE_STRUCT]->GetResidual_FEM_UTOL()));

//  }

  switch (kind_recording){
  case FLOW_CONS_VARS:      adjoint_convergence = flow_convergence; break;
  case MESH_COORDS:  adjoint_convergence = true; break;
  case FEA_DISP_VARS:       adjoint_convergence = struct_convergence; break;
  case FLOW_CROSS_TERM:     adjoint_convergence = true; break;
  case FEM_CROSS_TERM_GEOMETRY:      adjoint_convergence = true; break;
  case GEOMETRY_CROSS_TERM: adjoint_convergence = true; break;
  default:                  adjoint_convergence = false; break;
  }

  /*--- Apply the same convergence criteria to all the processors ---*/

#ifdef HAVE_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv[0] = adjoint_convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }

  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv[0] == 1) { adjoint_convergence = true;}
  else { adjoint_convergence = false;}

  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif

  return adjoint_convergence;

}

void CDiscAdjFSIDriver::ConvergenceHistory(unsigned long IntIter,
                                                      unsigned long nIntIter,
                                                      unsigned short ZONE_FLOW,
                                                      unsigned short ZONE_STRUCT,
                                                      unsigned short kind_recording){

  unsigned long BGS_Iter = config_container[ZONE_FLOW]->GetOuterIter();


  if (rank == MASTER_NODE)
     output_legacy->SetConvHistory_Header(&ConvHist_file[ZONE_0][INST_0], config_container[ZONE_0], ZONE_0, INST_0);

  if (kind_recording == FLOW_CONS_VARS) {

    if (rank == MASTER_NODE){
      if (IntIter == 0){
        cout << endl;
        cout << " IntIter" << "    BGSIter" << "   Res[Psi_Rho]" << "     Res[Psi_E]" << endl;
      }

      if (IntIter % config_container[ZONE_FLOW]->GetWrt_Con_Freq() == 0){
        /*--- Output the flow convergence ---*/
        /*--- This is temporary as it requires several changes in the output structure ---*/
        unsigned short nVar_Flow = solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetnVar();
        cout.width(8);     cout << IntIter;
        cout.width(11);    cout << BGS_Iter + 1;
        cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
        cout.width(15);    cout << log10(solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0));
        cout.width(15);    cout << log10(solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(nVar_Flow-1));
        cout << endl;
      }

    }
  }

  if (kind_recording == FEA_DISP_VARS) {
    /*--- Set the convergence criteria (only residual possible) ---*/
       output_legacy->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_STRUCT, INST_0);

  }


}


bool CDiscAdjFSIDriver::BGSConvergence(unsigned long IntIter,
                                                 unsigned short ZONE_FLOW,
                                                 unsigned short ZONE_STRUCT){

  unsigned short nVar_Flow = solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetnVar(),
                   nVar_Struct = solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetnVar();
  unsigned short iRes;

  bool flow_converged_absolute = false,
        flow_converged_relative = false,
        struct_converged_absolute = false,
        struct_converged_relative = false;

  bool Convergence = false;

  /*--- Compute the residual for the flow and structural zones ---*/

  /*--- Flow ---*/

  solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->ComputeResidual_Multizone(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                        config_container[ZONE_FLOW]);

  /*--- Structure ---*/

  solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->ComputeResidual_Multizone(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                         config_container[ZONE_STRUCT]);


  /*--- Retrieve residuals ---*/

  /*--- Flow residuals ---*/

  for (iRes = 0; iRes < nVar_Flow; iRes++){
    residual_flow[iRes] = log10(solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->GetRes_BGS(iRes));
    if (IntIter == 0) init_res_flow[iRes] = residual_flow[iRes];
    residual_flow_rel[iRes] = fabs(residual_flow[iRes] - init_res_flow[iRes]);
  }

  /*--- Structure residuals ---*/

  for (iRes = 0; iRes < nVar_Struct; iRes++){
    residual_struct[iRes] = log10(solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetRes_BGS(iRes));
    if (IntIter == 0) init_res_struct[iRes] = residual_struct[iRes];
    residual_struct_rel[iRes] = fabs(residual_struct[iRes] - init_res_struct[iRes]);
  }

  /*--- Check convergence ---*/
  flow_converged_absolute = ((residual_flow[0] < flow_criteria) && (residual_flow[nVar_Flow-1] < flow_criteria));
  flow_converged_relative = ((residual_flow_rel[0] > flow_criteria_rel) && (residual_flow_rel[nVar_Flow-1] > flow_criteria_rel));

  struct_converged_absolute = ((residual_struct[0] < structure_criteria) && (residual_struct[nVar_Flow-1] < structure_criteria));
  struct_converged_relative = ((residual_struct_rel[0] > structure_criteria_rel) && (residual_struct_rel[nVar_Flow-1] > structure_criteria_rel));

  Convergence = ((flow_converged_absolute && struct_converged_absolute) ||
                 (flow_converged_absolute && struct_converged_relative) ||
                 (flow_converged_relative && struct_converged_relative) ||
                 (flow_converged_relative && struct_converged_absolute));

  if (rank == MASTER_NODE){

    cout << "\n-------------------------------------------------------------------------\n\n";
    cout << "Convergence summary for BGS iteration " << IntIter << "\n\n";
    /*--- TODO: This is a workaround until the TestCases.py script incorporates new classes for nested loops. ---*/
    cout << "Iter[ID]" << "  BGSRes[Psi_Rho]" << "  BGSRes[Psi_E]" << "  BGSRes[Psi_Ux]" << "  BGSRes[Psi_Uy]\n";
    cout.precision(6); cout.setf(ios::fixed, ios::floatfield);
    cout << "|"; cout.width(8);  cout << IntIter*1000;
    cout << "|"; cout.width(17); cout << residual_flow[0];
    cout << "|"; cout.width(15); cout << residual_flow[nVar_Flow-1];
    cout << "|"; cout.width(16); cout << residual_struct[0];
    cout << "|"; cout.width(16); cout << residual_struct[1];
    cout << "|"; cout << "\n\n-------------------------------------------------------------------------" << endl;


    bool write_history = true;
    unsigned short iVar;

    /*--- Header of the temporary output file ---*/
    if ((write_history) && (rank == MASTER_NODE)){
      ofstream myfile_res;
      bool de_effects = config_container[ZONE_STRUCT]->GetDE_Effects();

      myfile_res.open ("history_adjoint_FSI.csv", ios::app);

      myfile_res << IntIter << "\t";

      myfile_res.precision(15);

      for (iVar = 0; iVar < nVar_Flow; iVar++){
        myfile_res << fixed << residual_flow[iVar] << "\t";
      }

      for (iVar = 0; iVar < nVar_Struct; iVar++){
        myfile_res << fixed << residual_struct[iVar] << "\t";
      }

      for (iVar = 0; iVar < config_container[ZONE_STRUCT]->GetnElasticityMod(); iVar++)
         myfile_res << scientific << solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_E(iVar) << "\t";
      for (iVar = 0; iVar < config_container[ZONE_STRUCT]->GetnPoissonRatio(); iVar++)
         myfile_res << scientific << solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_Nu(iVar) << "\t";
      if (de_effects){
        for (iVar = 0; iVar < config_container[ZONE_STRUCT]->GetnElectric_Field(); iVar++)
          myfile_res << scientific << solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_EField(0) << "\t";
      }

      myfile_res << endl;

      myfile_res.close();
    }

    // TEST: for implementation of python framework in coupled FSI problems
    if ((config_container[ZONE_1]->GetDV_FEA() != NODV_FEA) && (rank == MASTER_NODE)){

      /*--- Header of the temporary output file ---*/
      ofstream myfile_res;

      switch (config_container[ZONE_1]->GetDV_FEA()) {
        case YOUNG_MODULUS:
          myfile_res.open("grad_young.opt");
          break;
        case POISSON_RATIO:
          myfile_res.open("grad_poisson.opt");
          break;
        case DENSITY_VAL:
        case DEAD_WEIGHT:
          myfile_res.open("grad_density.opt");
          break;
        case ELECTRIC_FIELD:
          myfile_res.open("grad_efield.opt");
          break;
        default:
          myfile_res.open("grad.opt");
          break;
      }

      unsigned short iDV;
      unsigned short nDV = solver_container[ZONE_1][INST_0][MESH_0][ADJFEA_SOL]->GetnDVFEA();

      myfile_res << "INDEX" << "\t" << "GRAD" << endl;

      myfile_res.precision(15);

      for (iDV = 0; iDV < nDV; iDV++){
        myfile_res << iDV;
        myfile_res << "\t";
        myfile_res << scientific << solver_container[ZONE_1][INST_0][MESH_0][ADJFEA_SOL]->GetGlobal_Sens_DVFEA(iDV);
        myfile_res << endl;
      }

      myfile_res.close();
    }


  }

  /*--- Apply the same convergence criteria to all the processors ---*/

#ifdef HAVE_MPI

  unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
  sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
  rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv[0] = Convergence;
  SU2_MPI::Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv[0] = 0;
  if (rank == MASTER_NODE) {
    if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
    else sbuf_conv[0] = 0;
  }

  SU2_MPI::Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv[0] == 1) { Convergence = true;}
  else { Convergence = false;}

  delete [] sbuf_conv;
  delete [] rbuf_conv;

#endif

  /*--- Update the solution for the flow and structural zones ---*/

  /*--- Flow ---*/

  solver_container[ZONE_FLOW][INST_0][MESH_0][ADJFLOW_SOL]->UpdateSolution_BGS(geometry_container[ZONE_FLOW][INST_0][MESH_0],
                                                                       config_container[ZONE_FLOW]);

  /*--- Structure ---*/

  solver_container[ZONE_STRUCT][INST_0][MESH_0][ADJFEA_SOL]->UpdateSolution_BGS(geometry_container[ZONE_STRUCT][INST_0][MESH_0],
                                                                       config_container[ZONE_STRUCT]);

  return Convergence;
}

void CDiscAdjFSIDriver::Transfer_Displacements(unsigned short donorZone, unsigned short targetZone) {


  interface_container[donorZone][targetZone]->BroadcastData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                     geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                     config_container[donorZone], config_container[targetZone]);

}

void CDiscAdjFSIDriver::Transfer_Tractions(unsigned short donorZone, unsigned short targetZone) {

  interface_container[donorZone][targetZone]->BroadcastData(solver_container[donorZone][INST_0][MESH_0][FEA_SOL],solver_container[targetZone][INST_0][MESH_0][FLOW_SOL],
                                                                     geometry_container[donorZone][INST_0][MESH_0],geometry_container[targetZone][INST_0][MESH_0],
                                                                     config_container[donorZone], config_container[targetZone]);
}
