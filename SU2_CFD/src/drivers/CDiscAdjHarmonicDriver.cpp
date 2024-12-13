/*!
 * \file CDiscAdjHarmonicDriver.cpp
 * \brief The main subroutines for driving adjoint harmonic balance single-zone problems.
 * \author 
 * \version 7.3.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/drivers/CDiscAdjHarmonicDriver.hpp"
#include "../../include/output/COutputLegacy.hpp"
#include "../../include/output/COutputFactory.hpp"
#include "../../include/output/COutput.hpp"
#include "../../include/iteration/CIterationFactory.hpp"
#include "../../include/iteration/CTurboIteration.hpp"


CDiscAdjHarmonicDriver::CDiscAdjHarmonicDriver(char* confFile,
                                                   unsigned short val_nZone,
                                                   SU2_Comm MPICommunicator) : CSinglezoneDriver(confFile,
                                                                                                 val_nZone,
                                                                                                 MPICommunicator) {


  /*--- Store the number of internal iterations that will be run by the adjoint solver ---*/
  nAdjoint_Iter = config_container[ZONE_0]->GetnInner_Iter();

  /*--- Store the pointers ---*/
  config      = config_container[ZONE_0];
  iteration   = iteration_container[ZONE_0][INST_0];
  /*solver      = solver_container[ZONE_0][INST_0][MESH_0];
  numerics    = numerics_container[ZONE_0][INST_0][MESH_0];
  geometry    = geometry_container[ZONE_0][INST_0][MESH_0];
  integration = integration_container[ZONE_0][INST_0];*/


  nInstHB = config->GetnTimeInstances();
  direct_iteration   = new CIteration*    [nInstHB] ();

  D = nullptr;
  /*--- allocate dynamic memory for the Harmonic Balance operator ---*/
  D = new su2double*[nInstHB]; for (auto kInst = 0; kInst < nInstHB; kInst++) D[kInst] = new su2double[nInstHB];

  /*--- Store the recording state ---*/
  RecordingState = RECORDING::CLEAR_INDICES;

  /*--- Initialize the direct iteration ---*/

  switch (config->GetKind_Solver()) {

  case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
  case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;

    if (config->GetBoolTurbomachinery()) {
      for (auto iInst = 0; iInst < nInstHB; iInst++) {
        direct_iteration[iInst] = new CTurboIteration(config);
      }
    }
    else {
      for (auto iInst = 0; iInst < nInstHB; iInst++) {
        direct_iteration[iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::EULER, config);
      }
    }
    output_legacy = COutputFactory::CreateLegacyOutput(config_container[ZONE_0]);

    if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {
      direct_output = COutputFactory::CreateOutput(MAIN_SOLVER::EULER, config, nDim);
    }
    else {
      direct_output =  COutputFactory::CreateOutput(MAIN_SOLVER::INC_EULER, config, nDim);
    }

    MainVariables = RECORDING::SOLUTION_VARIABLES;
    if (config->GetDeform_Mesh()) {
      SecondaryVariables = RECORDING::MESH_DEFORM;
    }
    else { SecondaryVariables = RECORDING::MESH_COORDS; }
    MainSolver = ADJFLOW_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM_EULER : case MAIN_SOLVER::DISC_ADJ_FEM_NS : case MAIN_SOLVER::DISC_ADJ_FEM_RANS :
    if (rank == MASTER_NODE)
      cout << "Direct iteration: Euler/Navier-Stokes/RANS equation." << endl;
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
        direct_iteration[iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_EULER, config);
      }
    output_legacy = COutputFactory::CreateLegacyOutput(config_container[ZONE_0]);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJFLOW_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_FEM:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: elasticity equation." << endl;
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
        direct_iteration[iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::FEM_ELASTICITY, config);
      }
    output_legacy = COutputFactory::CreateLegacyOutput(config_container[ZONE_0]);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJFEA_SOL;
    break;

  case MAIN_SOLVER::DISC_ADJ_HEAT:
    if (rank == MASTER_NODE)
      cout << "Direct iteration: heat equation." << endl;
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
        direct_iteration[iInst] = CIterationFactory::CreateIteration(MAIN_SOLVER::HEAT_EQUATION, config);
      }
    output_legacy = COutputFactory::CreateLegacyOutput(config_container[ZONE_0]);
    MainVariables = RECORDING::SOLUTION_VARIABLES;
    SecondaryVariables = RECORDING::MESH_COORDS;
    MainSolver = ADJHEAT_SOL;
    break;

  default:
    break;

  }
  direct_output->PreprocessHistoryOutput(config, false);

}

CDiscAdjHarmonicDriver::~CDiscAdjHarmonicDriver(void) {

  /*--- delete dynamic memory for the Harmonic Balance operator ---*/
  /*for (auto kInst = 0; kInst < nInstHB; kInst++) delete [] D[kInst];
  delete [] D;

  for (auto kInst = 0; kInst < nInstHB; kInst++){
    delete direct_iteration[iInst];
  }

  delete [] direct_iteration;*/
  delete direct_output;

}

void CDiscAdjHarmonicDriver::Preprocess(unsigned long TimeIter) {
  config_container[ZONE_0]->SetTimeIter(TimeIter);

  /*--- Preprocess the adjoint iteration ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    iteration_container[ZONE_0][iInst]->Preprocess(output_container[ZONE_0], integration_container, geometry_container,
                          solver_container, numerics_container, config_container,
                          surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
  }
  /*--- For the adjoint iteration we need the derivatives of the iteration function with
   *--- respect to the conservative variables. Since these derivatives do not change in the steady state case
   *--- we only have to record if the current recording is different from the main variables. ---*/

  if (RecordingState != MainVariables){
    MainRecording();
  }
}

void CDiscAdjHarmonicDriver::Run() {

  for (auto Adjoint_Iter = 0ul; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) {

    /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *--- of the previous iteration. The values are passed to the AD tool.
     *--- Issues with iteration number should be dealt with once the output structure is in place. ---*/

    config->SetInnerIter(Adjoint_Iter);
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
      iteration_container[ZONE_0][iInst]->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, iInst);
    }

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
      iteration_container[ZONE_0][iInst]->IterateDiscAdj(geometry_container, solver_container,
                                config_container, ZONE_0, iInst, false);
    }
    /*--- Monitor the pseudo-time ---*/
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
      StopCalc = iteration->Monitor(output_container[ZONE_0], integration_container, geometry_container,
                                    solver_container, numerics_container, config_container,
                                    surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
    }

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    /*--- Output files for steady state simulations. ---*/

    if (!config->GetTime_Domain()) {
      for (auto iInst = 0; iInst < nInstHB; iInst++) {
        config_container[ZONE_0]->SetiInst(iInst);
        iteration_container[ZONE_0][iInst]->Output(output_container[ZONE_0], geometry_container, solver_container,
                        config_container, Adjoint_Iter, false, ZONE_0, iInst);
        /*output_container[ZONE_0]->SetResult_Files(geometry_container[ZONE_0][iInst][MESH_0],
                                               config_container[ZONE_0],
                                               solver_container[ZONE_0][iInst][MESH_0],
                                               Adjoint_Iter, StopCalc);*/
      }
    }
    if (StopCalc) break;

  }

}

void CDiscAdjHarmonicDriver::Postprocess() {

  switch(config->GetKind_Solver())
  {
    case MAIN_SOLVER::DISC_ADJ_EULER :     case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES :     case MAIN_SOLVER::DISC_ADJ_RANS :
    case MAIN_SOLVER::DISC_ADJ_INC_EULER : case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES : case MAIN_SOLVER::DISC_ADJ_INC_RANS :
    case MAIN_SOLVER::DISC_ADJ_HEAT :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();
      break;

    case MAIN_SOLVER::DISC_ADJ_FEM :

      /*--- Compute the geometrical sensitivities ---*/
      SecondaryRecording();

      iteration->Postprocess(output_container[ZONE_0], integration_container, geometry_container,
                             solver_container, numerics_container, config_container,
                             surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);
      break;

    default:
      break;

  }//switch

}

void CDiscAdjHarmonicDriver::SetRecording(RECORDING kind_recording){
  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution. ---*/

  for (unsigned short iSol=0; iSol < MAX_SOLS; iSol++) {
    for (unsigned short iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++) {
      for (unsigned short iInst=INST_0; iInst < nInstHB; iInst++) {
        auto solver = solver_container[ZONE_0][iInst][iMesh][iSol];
        if (solver && solver->GetAdjoint()) {
          solver->SetRecording(geometry_container[ZONE_0][iInst][iMesh], config_container[ZONE_0]);
        }
      }
    }
  }

  if (rank == MASTER_NODE) {
    cout << "\n-------------------------------------------------------------------------\n";
    switch(kind_recording) {
    case RECORDING::CLEAR_INDICES: cout << "Clearing the computational graph." << endl; break;
    case RECORDING::MESH_COORDS:   cout << "Storing computational graph wrt MESH COORDINATES." << endl; break;
    case RECORDING::SOLUTION_VARIABLES:
      cout << "Direct iteration to store the primal computational graph." << endl;
      cout << "Computing residuals to check the convergence of the direct problem." << endl; break;
    default: break;
    }
  }

  /*---Enable recording and register input of the iteration --- */

  if (kind_recording != RECORDING::CLEAR_INDICES){

    AD::StartRecording();
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
      iteration_container[ZONE_0][iInst]->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, iInst, kind_recording);
    }
  }

  /*--- Set the dependencies of the iteration ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    iteration_container[ZONE_0][iInst]->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0,
                              iInst, kind_recording);
  }
  /*--- Do one iteration of the direct solver ---*/

  DirectRun(kind_recording);

  /*--- Store the recording state ---*/

  RecordingState = kind_recording;

  /*--- Register Output of the iteration ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    iteration_container[ZONE_0][iInst]->RegisterOutput(solver_container, geometry_container, config_container, ZONE_0, iInst);
  }
  /*--- Extract the objective function and store it --- */

  SetObjFunction();

  if (kind_recording != RECORDING::CLEAR_INDICES && config_container[ZONE_0]->GetWrt_AD_Statistics()) {
    if (rank == MASTER_NODE) AD::PrintStatistics();
#ifdef CODI_REVERSE_TYPE
    if (size > SINGLE_NODE) {
      su2double myMem = AD::getGlobalTape().getTapeValues().getUsedMemorySize(), totMem = 0.0;
      SU2_MPI::Allreduce(&myMem, &totMem, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      if (rank == MASTER_NODE) {
        cout << "MPI\n";
        cout << "-------------------------------------\n";
        cout << "  Total memory used      :  " << totMem << " MB\n";
        cout << "-------------------------------------\n" << endl;
      }
    }
#endif
  }

  AD::StopRecording();

}

void CDiscAdjHarmonicDriver::SetAdj_ObjFunction(){

  su2double seeding = 1.0; // No need for windowing

  if (rank == MASTER_NODE){
    SU2_TYPE::SetDerivative(ObjFunc, SU2_TYPE::GetValue(seeding));
  } else {
    SU2_TYPE::SetDerivative(ObjFunc, 0.0);
  }
}

void CDiscAdjHarmonicDriver::SetObjFunction(){

  ObjFunc = 0.0;

  /*--- Specific scalar objective functions ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    switch (config->GetKind_Solver()) {
    case MAIN_SOLVER::DISC_ADJ_INC_EULER:       case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:      case MAIN_SOLVER::DISC_ADJ_INC_RANS:
    case MAIN_SOLVER::DISC_ADJ_EULER:           case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:          case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_FEM_EULER:       case MAIN_SOLVER::DISC_ADJ_FEM_NS:                 case MAIN_SOLVER::DISC_ADJ_FEM_RANS:
      direct_output->SetHistory_Output(geometry_container[ZONE_0][iInst][MESH_0], solver_container[ZONE_0][iInst][MESH_0], config, config->GetTimeIter(),
                                     config->GetOuterIter(), config->GetInnerIter());
      /*--- Surface based obj. function ---*/
      ObjFunc += solver_container[ZONE_0][iInst][MESH_0][FLOW_SOL]->GetTotal_ComboObj();

      /*--- These calls to be moved to a generic framework at a next stage        ---*/
      /*--- Some things that are currently hacked into output must be reorganized ---*/
      if (config->GetBoolTurbomachinery()) {
        output_legacy->ComputeTurboPerformance(solver_container[ZONE_0][iInst][MESH_0][FLOW_SOL], geometry_container[ZONE_0][iInst][MESH_0], config);

        unsigned short nMarkerTurboPerf = config->GetnMarker_TurboPerformance();
        unsigned short nSpanSections = config->GetnSpanWiseSections();

        switch (config_container[ZONE_0]->GetKind_ObjFunc()){
        case ENTROPY_GENERATION:
          ObjFunc += output_legacy->GetEntropyGen(nMarkerTurboPerf-1, nSpanSections);
          break;
        case FLOW_ANGLE_OUT:
          ObjFunc += output_legacy->GetFlowAngleOut(nMarkerTurboPerf-1, nSpanSections);
          break;
        case MASS_FLOW_IN:
          ObjFunc += output_legacy->GetMassFlowIn(nMarkerTurboPerf-1, nSpanSections);
          break;
        default:
          break;
        }
      }
      break;

    case MAIN_SOLVER::DISC_ADJ_HEAT:
      ObjFunc = solver_container[ZONE_0][iInst][MESH_0][HEAT_SOL]->GetTotal_ComboObj();
      break;

    case MAIN_SOLVER::DISC_ADJ_FEM:
      solver_container[ZONE_0][iInst][MESH_0][FEA_SOL]->Postprocessing(geometry_container[ZONE_0][iInst][MESH_0], config, numerics_container[ZONE_0][iInst][MESH_0][FEA_SOL], true);

      ObjFunc = solver_container[ZONE_0][iInst][MESH_0][FEA_SOL]->GetTotal_ComboObj();
      break;

    default:
      break;
    }
  }

  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc);
  }

}

void CDiscAdjHarmonicDriver::DirectRun(RECORDING kind_recording){

  ComputeHB_Operator();

  /*--- Mesh movement ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    direct_iteration[iInst]->SetMesh_Deformation(geometry_container[ZONE_0][iInst], solver_container[ZONE_0][iInst][MESH_0], numerics_container[ZONE_0][iInst][MESH_0], config, kind_recording);
  }

  for(auto iInst = 0; iInst < nInstHB; iInst++) {
    for (unsigned long jPoint = 0; jPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); jPoint++) {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        su2double GridVel = 0.0;
        for (auto kInst = 0; kInst < nInstHB; kInst++) {
          const su2double Disp = geometry_container[ZONE_0][kInst][MESH_0]->nodes->GetCoord(jPoint, iDim);
          GridVel += D[iInst][kInst]*Disp;
        }
        geometry_container[ZONE_0][iInst][MESH_0]->nodes->SetGridVel(jPoint, iDim, GridVel);
      }
    }
    /*--- The velocity was computed for nPointDomain, now we communicate it. ---*/
    geometry_container[ZONE_0][iInst][MESH_0]->InitiateComms(geometry_container[ZONE_0][iInst][MESH_0], config_container[ZONE_0], GRID_VELOCITY);
    geometry_container[ZONE_0][iInst][MESH_0]->CompleteComms(geometry_container[ZONE_0][iInst][MESH_0], config_container[ZONE_0], GRID_VELOCITY);
  }

  if(rank == MASTER_NODE) cout << " Updating multigrid structure." << endl;
  for(unsigned int jInst = 0; jInst < nInstHB; jInst++) {
    grid_movement[ZONE_0][jInst]->UpdateMultiGrid(geometry_container[ZONE_0][jInst], config_container[ZONE_0]);
  }

  /*--- Zone preprocessing ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    direct_iteration[iInst]->Preprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
  }
  /*--- Iterate the direct solver ---*/
  SetHarmonicBalance(false);
  /*--- Precondition the harmonic balance source terms ---*/
  if (config_container[ZONE_0]->GetHB_Precondition() == YES) {
    StabilizeHarmonicBalance();
  }
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    direct_iteration[iInst]->Iterate(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
  }

  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    direct_iteration[iInst]->Update(direct_output, integration_container, geometry_container,
        solver_container, numerics_container, config_container,
        surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
  }
  /*--- Postprocess the direct solver ---*/
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    direct_iteration[iInst]->Postprocess(direct_output, integration_container, geometry_container, solver_container, numerics_container, config_container, surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
  }

  /*--- Print the direct residual to screen ---*/

  Print_DirectResidual(kind_recording);

}

void CDiscAdjHarmonicDriver::MainRecording(){
  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with
   *    RECORDING::CLEAR_INDICES as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES);

  /*--- Store the computational graph of one direct iteration with the solution variables as input. ---*/

  SetRecording(MainVariables);

}

void CDiscAdjHarmonicDriver::SecondaryRecording(){
  /*--- SetRecording stores the computational graph on one iteration of the direct problem. Calling it with
   *    RECORDING::CLEAR_INDICES as argument ensures that all information from a previous recording is removed. ---*/

  SetRecording(RECORDING::CLEAR_INDICES);

  /*--- Store the computational graph of one direct iteration with the secondary variables as input. ---*/

  SetRecording(SecondaryVariables);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the current iteration. The values are passed to the AD tool. ---*/

  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    iteration_container[ZONE_0][iInst]->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, iInst);
  }

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed sensitivity values. ---*/

  if (SecondaryVariables == RECORDING::MESH_COORDS) {
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
      solver_container[ZONE_0][iInst][MESH_0][MainSolver]->SetSensitivity(geometry_container[ZONE_0][iInst][MESH_0], config);
    }
  }
  else { // MESH_DEFORM
    for (auto iInst = 0; iInst < nInstHB; iInst++) {
      solver_container[ZONE_0][iInst][MESH_0][ADJMESH_SOL]->SetSensitivity(geometry_container[ZONE_0][iInst][MESH_0], config, solver_container[ZONE_0][iInst][MESH_0][MainSolver]);
    }
  }

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();

}

void CDiscAdjHarmonicDriver::ComputeHB_Operator() {

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
    //Omega_HB[iInst] /= config_container[ZONE_0]->GetOmega_Ref(); //TODO: check
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

void CDiscAdjHarmonicDriver::SetHarmonicBalance(bool implicit) {

  unsigned short iVar, iInst, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;
  unsigned long InnerIter = config_container[ZONE_0]->GetInnerIter();

  /*--- Retrieve values from the config file ---*/
  su2double *U = new su2double[nInstHB];
  su2double *UV = new su2double[nInstHB];
  su2double *U_old = new su2double[nInstHB];
  su2double *UV_old = new su2double[nInstHB];
  su2double *Psi = new su2double[nInstHB];
  su2double *PsiV = new su2double[nInstHB];
  su2double *Psi_old = new su2double[nInstHB];
  su2double *PsiV_old = new su2double[nInstHB];
  su2double *Source = new su2double[nInstHB];
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
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][iMGlevel]->GetnPoint(); iPoint++) {
      /*--- Loop over all variables ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        /*--- Calculate volume source term ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          su2double Volume = geometry_container[ZONE_0][iInst][iMGlevel]->nodes->GetVolume(iPoint);
          Source[iInst] = 0.0;
          U[iInst] = solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution(iPoint, iVar);
          UV[iInst] = Volume*U[iInst];
          if (implicit) {
            U_old[iInst] = solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
            UV_old[iInst] = Volume*U_old[iInst];
          }
        }
        /*--- Loop over the time instances ---*/
        for (iInst = 0; iInst < nInstHB; iInst++) {
          /*--- Step across the columns ---*/
          for (jInst = 0; jInst < nInstHB; jInst++) {
            /*--- Retrieve solution at this node in current zone ---*/
            Source[iInst] += UV[jInst]*D[iInst][jInst];
            if (implicit) {
              deltaU = UV[jInst] - UV_old[jInst];
              Source[iInst] += deltaU*D[iInst][jInst];
            }
          }
          solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
        }
      }
    }
  }

  /*--- Source term for a turbulence model ---*/
  if (config_container[ZONE_0]->GetKind_Solver() == MAIN_SOLVER::RANS) {

    /*--- Extra variables needed if we have a turbulence model. ---*/
    unsigned short nVar_Turb = solver_container[ZONE_0][INST_0][MESH_0][TURB_SOL]->GetnVar();

    /*--- Loop over only the finest mesh level (turbulence is always solved
     on the original grid only). ---*/
    for (iPoint = 0; iPoint < geometry_container[ZONE_0][INST_0][MESH_0]->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar_Turb; iVar++){
        for (iInst = 0; iInst < nInstHB; iInst++) {
          su2double Volume = geometry_container[ZONE_0][iInst][MESH_0]->nodes->GetVolume(iPoint);
          Source[iInst] = 0.0;
          U[iInst] = solver_container[ZONE_0][iInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution(iPoint, iVar);
          UV[iInst] = Volume*U[iInst];
          if (implicit) {
            U_old[iInst] = solver_container[ZONE_0][iInst][MESH_0][TURB_SOL]->GetNodes()->GetSolution_Old(iPoint, iVar);
            UV_old[iInst] = Volume*U_old[iInst];
          }
        }
        for (iInst = 0; iInst < nInstHB; iInst++) {
          /*--- Step across the columns ---*/
          for (jInst = 0; jInst < nInstHB; jInst++) {
            /*--- Retrieve solution at this node in current zone ---*/
            Source[iInst] += UV[jInst]*D[iInst][jInst];
            if (implicit) {
              deltaU = UV[jInst] - UV_old[jInst];
              Source[iInst] += deltaU*D[iInst][jInst];
            }
          }
          /*--- Store sources for current iInst ---*/
          solver_container[ZONE_0][iInst][MESH_0][TURB_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
        }
      }
    }
  }

  delete [] Source;
  delete [] U;
  delete [] UV;
  delete [] U_old;
  delete [] UV_old;
  delete [] Psi;
  delete [] PsiV;
  delete [] Psi_old;
  delete [] PsiV_old;

}

void CDiscAdjHarmonicDriver::StabilizeHarmonicBalance() {

  unsigned short i, j, k, iVar, iInst, jInst, iMGlevel;
  unsigned short nVar = solver_container[ZONE_0][INST_0][MESH_0][FLOW_SOL]->GetnVar();
  unsigned long iPoint;

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
            su2double Volume = geometry_container[ZONE_0][jInst][iMGlevel]->nodes->GetVolume(iPoint);
            Source[iInst] += P[iInst][jInst]*Source_old[jInst];
          }

          /*--- Store updated source terms for current node ---*/
          solver_container[ZONE_0][iInst][iMGlevel][FLOW_SOL]->GetNodes()->SetHarmonicBalance_Source(iPoint, iVar, Source[iInst]);
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

void CDiscAdjHarmonicDriver::Output(unsigned long TimeIter) {

  /*--- Time the output for performance benchmarking. ---*/
  bool wrote_files = false;
  StopTime = SU2_MPI::Wtime();

  UsedTimeCompute += StopTime-StartTime;

  StartTime = SU2_MPI::Wtime();
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    config_container[ZONE_0]->SetiInst(iInst);
    wrote_files = output_container[ZONE_0]->SetResult_Files(geometry_container[ZONE_0][iInst][MESH_0],
                                                                config_container[ZONE_0],
                                                                solver_container[ZONE_0][iInst][MESH_0],
                                                                TimeIter, StopCalc);
  }

  if (wrote_files){

    StopTime = SU2_MPI::Wtime();

    UsedTimeOutput += StopTime-StartTime;
    OutputCount++;
    BandwidthSum = config_container[ZONE_0]->GetRestart_Bandwidth_Agg();

    StartTime = SU2_MPI::Wtime();

    config_container[ZONE_0]->Set_StartTime(StartTime);
  }
}

void CDiscAdjHarmonicDriver::Update() {
  for (auto iInst = 0; iInst < nInstHB; iInst++) {
    iteration_container[ZONE_0][iInst]->Update(output_container[ZONE_0], integration_container, geometry_container,
          solver_container, numerics_container, config_container,
          surface_movement, grid_movement, FFDBox, ZONE_0, iInst);
  }

}
