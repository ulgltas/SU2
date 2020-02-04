/*!
 * \file integration_time.cpp
 * \brief Time dependent numerical methods
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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


#include "../include/integration_structure.hpp"

CMultiGridIntegration::CMultiGridIntegration(CConfig *config) : CIntegration(config) {}

CMultiGridIntegration::~CMultiGridIntegration(void) { }

void CMultiGridIntegration::MultiGrid_Iteration(CGeometry ****geometry,
                                                CSolver *****solver_container,
                                                CNumerics ******numerics_container,
                                                CConfig **config,
                                                unsigned short RunTime_EqSystem,
                                                unsigned short iZone,
                                                unsigned short iInst) {
  unsigned short FinestMesh;
  su2double monitor = 1.0;
  bool FullMG = false;
  
  const bool direct = ((config[iZone]->GetKind_Solver() == EULER)                         ||
                       (config[iZone]->GetKind_Solver() == NAVIER_STOKES)                 ||
                       (config[iZone]->GetKind_Solver() == RANS)                          ||
                       (config[iZone]->GetKind_Solver() == FEM_EULER)                     ||
                       (config[iZone]->GetKind_Solver() == FEM_NAVIER_STOKES)             ||
                       (config[iZone]->GetKind_Solver() == FEM_RANS)                      ||
                       (config[iZone]->GetKind_Solver() == FEM_LES)                       ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_EULER)                ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)        ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_FEM_EULER)            ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_FEM_NS)               ||
                       (config[iZone]->GetKind_Solver() == DISC_ADJ_RANS));
  const unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  unsigned short RecursiveParam = config[iZone]->GetMGCycle();
  
  if (config[iZone]->GetMGCycle() == FULLMG_CYCLE) {
    RecursiveParam = V_CYCLE;
    FullMG = true;
  }
	
  /*--- Full multigrid strategy and start up with fine grid only works with the direct problem ---*/

  if (!config[iZone]->GetRestart() && FullMG && direct && ( Convergence_FullMG && (config[iZone]->GetFinestMesh() != MESH_0 ))) {
    SetProlongated_Solution(RunTime_EqSystem, solver_container[iZone][iInst][config[iZone]->GetFinestMesh()-1][SolContainer_Position],
                            solver_container[iZone][iInst][config[iZone]->GetFinestMesh()][SolContainer_Position],
                            geometry[iZone][iInst][config[iZone]->GetFinestMesh()-1], geometry[iZone][iInst][config[iZone]->GetFinestMesh()],
                            config[iZone]);
    config[iZone]->SubtractFinestMesh();
  }

  /*--- Set the current finest grid (full multigrid strategy) ---*/
  
  FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Perform the Full Approximation Scheme multigrid ---*/
  
  MultiGrid_Cycle(geometry, solver_container, numerics_container, config,
                  FinestMesh, RecursiveParam, RunTime_EqSystem,
                  iZone, iInst);

  /*--- Computes primitive variables and gradients in the finest mesh (useful for the next solver (turbulence) and output ---*/

   solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][MESH_0],
                                                                         solver_container[iZone][iInst][MESH_0], config[iZone],
                                                                         MESH_0, NO_RK_ITER, RunTime_EqSystem, true);
  
  /*--- Compute non-dimensional parameters and the convergence monitor ---*/
  
  NonDimensional_Parameters(geometry[iZone][iInst], solver_container[iZone][iInst],
                            numerics_container[iZone][iInst], config[iZone],
                            FinestMesh, RunTime_EqSystem, &monitor);
  
}

void CMultiGridIntegration::MultiGrid_Cycle(CGeometry ****geometry,
                                            CSolver *****solver_container,
                                            CNumerics ******numerics_container,
                                            CConfig **config,
                                            unsigned short iMesh,
                                            unsigned short RecursiveParam,
                                            unsigned short RunTime_EqSystem,
                                            unsigned short iZone,
                                            unsigned short iInst) {
  
  unsigned short iPreSmooth, iPostSmooth, iRKStep, iRKLimit = 1;
    unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  
  /*--- Do a presmoothing on the grid iMesh to be restricted to the grid iMesh+1 ---*/
  
  for (iPreSmooth = 0; iPreSmooth < config[iZone]->GetMG_PreSmooth(iMesh); iPreSmooth++) {
    
    switch (config[iZone]->GetKind_TimeIntScheme()) {
      case RUNGE_KUTTA_EXPLICIT: iRKLimit = config[iZone]->GetnRKStep(); break;
      case CLASSICAL_RK4_EXPLICIT: iRKLimit = 4; break;
      case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; }

    /*--- Time and space integration ---*/
    
    for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
      
      /*--- Send-Receive boundary conditions, and preprocessing ---*/
      
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, iRKStep, RunTime_EqSystem, false);
      
      if (iRKStep == 0) {
        
        /*--- Set the old solution ---*/
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][iMesh]);

        if (config[iZone]->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT)
          solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_NewSolution(geometry[iZone][iInst][iMesh]);

        /*--- Compute time step, max eigenvalue, and integration scheme (steady and unsteady problems) ---*/
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh,  config[iZone]->GetTimeIter());
        
        /*--- Restrict the solution and gradient for the adjoint problem ---*/
        
        Adjoint_Setup(geometry, solver_container, config, RunTime_EqSystem, config[iZone]->GetTimeIter(), iZone);
        
      }
      
      /*--- Space integration ---*/
      
      Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, iRKStep, RunTime_EqSystem);
      
      /*--- Time integration, update solution using the old solution plus the solution increment ---*/
      
      Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iRKStep, RunTime_EqSystem);
      
      /*--- Send-Receive boundary conditions, and postprocessing ---*/
      
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
      
    }
    
  }
  
  /*--- Compute Forcing Term $P_(k+1) = I^(k+1)_k(P_k+F_k(u_k))-F_(k+1)(I^(k+1)_k u_k)$ and update solution for multigrid ---*/
  if ( iMesh < config[iZone]->GetnMGLevels() ) {
    /*--- Compute $r_k = P_k + F_k(u_k)$ ---*/
    
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, NO_RK_ITER, RunTime_EqSystem, false);
    Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, NO_RK_ITER, RunTime_EqSystem);
    SetResidual_Term(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh][SolContainer_Position]);
    
    /*--- Compute $r_(k+1) = F_(k+1)(I^(k+1)_k u_k)$ ---*/
    
    SetRestricted_Solution(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    solver_container[iZone][iInst][iMesh+1][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], config[iZone], iMesh+1, NO_RK_ITER, RunTime_EqSystem, false);
    Space_Integration(geometry[iZone][iInst][iMesh+1], solver_container[iZone][iInst][iMesh+1], numerics_container[iZone][iInst][iMesh+1][SolContainer_Position], config[iZone], iMesh+1, NO_RK_ITER, RunTime_EqSystem);
    
    /*--- Compute $P_(k+1) = I^(k+1)_k(r_k) - r_(k+1) ---*/
    
    SetForcing_Term(solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone], iMesh+1);
    
    /*--- Recursive call to MultiGrid_Cycle ---*/
    
    for (unsigned short imu = 0; imu <= RecursiveParam; imu++) {
      if (iMesh == config[iZone]->GetnMGLevels()-2) MultiGrid_Cycle(geometry, solver_container, numerics_container, config, iMesh+1, 0, RunTime_EqSystem, iZone, iInst);
      else MultiGrid_Cycle(geometry, solver_container, numerics_container, config, iMesh+1, RecursiveParam, RunTime_EqSystem, iZone, iInst);
    }
    
    /*--- Compute prolongated solution, and smooth the correction $u^(new)_k = u_k +  Smooth(I^k_(k+1)(u_(k+1)-I^(k+1)_k u_k))$ ---*/
    
    GetProlongated_Correction(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position],
                              geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    SmoothProlongated_Correction(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], geometry[iZone][iInst][iMesh],
                                 config[iZone]->GetMG_CorrecSmooth(iMesh), 1.25, config[iZone]);
    SetProlongated_Correction(solver_container[iZone][iInst][iMesh][SolContainer_Position], geometry[iZone][iInst][iMesh], config[iZone], iMesh);
    
    /*--- Solution postsmoothing in the prolongated grid ---*/
    
    for (iPostSmooth = 0; iPostSmooth < config[iZone]->GetMG_PostSmooth(iMesh); iPostSmooth++) {
      
      switch (config[iZone]->GetKind_TimeIntScheme()) {
        case RUNGE_KUTTA_EXPLICIT: iRKLimit = config[iZone]->GetnRKStep(); break;
        case CLASSICAL_RK4_EXPLICIT: iRKLimit = 4; break;
        case EULER_EXPLICIT: case EULER_IMPLICIT: iRKLimit = 1; break; }

      for (iRKStep = 0; iRKStep < iRKLimit; iRKStep++) {
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh, iRKStep, RunTime_EqSystem, false);
        
        if (iRKStep == 0) {
          solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][iMesh]);
          if (config[iZone]->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT)
            solver_container[iZone][iInst][iMesh][SolContainer_Position]->Set_NewSolution(geometry[iZone][iInst][iMesh]);
          solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh,  config[iZone]->GetTimeIter());
        }
        
        Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], numerics_container[iZone][iInst][iMesh][SolContainer_Position], config[iZone], iMesh, iRKStep, RunTime_EqSystem);
        Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iRKStep, RunTime_EqSystem);
        
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh], config[iZone], iMesh);
        
      }
    }
  }
  
}

void CMultiGridIntegration::GetProlongated_Correction(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                                      CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short Boundary, iMarker, iChildren, iVar;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Solution_Coarse;
  
  const unsigned short nVar = sol_coarse->GetnVar();
  
  su2double *Solution = new su2double[nVar];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] -= Solution_Fine[iVar]*Area_Children/Area_Parent;
    }
    
    Solution_Coarse = sol_coarse->GetNodes()->GetSolution(Point_Coarse);
    
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] += Solution_Coarse[iVar];
    
    for (iVar = 0; iVar < nVar; iVar++)
      sol_coarse->GetNodes()->SetSolution_Old(Point_Coarse,Solution);
    
  }
  
  /*--- Remove any contributions from no-slip walls. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    if ((Boundary == HEAT_FLUX             ) ||
        (Boundary == ISOTHERMAL            ) ||
        (Boundary == CHT_WALL_INTERFACE    )) {
      
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        
        /*--- For dirichlet boundary condtions, set the correction to zero.
         Note that Solution_Old stores the correction not the actual value ---*/
        
        sol_coarse->GetNodes()->SetVelSolutionOldZero(Point_Coarse);
        
      }
      
    }
  }
  
  /*--- MPI the set solution old ---*/
  
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION_OLD);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION_OLD);

  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      sol_fine->LinSysRes.SetBlock(Point_Fine, sol_coarse->GetNodes()->GetSolution_Old(Point_Coarse));
    }
  }
  
  delete [] Solution;
  
}

void CMultiGridIntegration::SmoothProlongated_Correction (unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                                          unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  su2double *Residual_Old, *Residual_Sum, *Residual, *Residual_i, *Residual_j;
  unsigned short iVar, iSmooth, iMarker, nneigh;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  
  const unsigned short nVar = solver->GetnVar();
  
  if (val_nSmooth > 0) {
    
    Residual = new su2double [nVar];
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Residual_Old = solver->LinSysRes.GetBlock(iPoint);
      solver->GetNodes()->SetResidual_Old(iPoint,Residual_Old);
    }
    
    /*--- Jacobi iterations ---*/
    
    for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
      solver->GetNodes()->SetResidualSumZero();
      
      /*--- Loop over Interior edges ---*/
      
      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        
        Residual_i = solver->LinSysRes.GetBlock(iPoint);
        Residual_j = solver->LinSysRes.GetBlock(jPoint);
        
        /*--- Accumulate nearest neighbor Residual to Res_sum for each variable ---*/
        
        solver->GetNodes()->AddResidual_Sum(iPoint,Residual_j);
        solver->GetNodes()->AddResidual_Sum(jPoint,Residual_i);
      }
      
      /*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        nneigh = geometry->node[iPoint]->GetnPoint();
        Residual_Sum = solver->GetNodes()->GetResidual_Sum(iPoint);
        Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);
        for (iVar = 0; iVar < nVar; iVar++) {
          Residual[iVar] =(Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar])
          /(1.0 + val_smooth_coeff*su2double(nneigh));
        }
        solver->LinSysRes.SetBlock(iPoint, Residual);
      }
      
      /*--- Copy boundary values ---*/
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
        if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
            (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Residual_Old = solver->GetNodes()->GetResidual_Old(iPoint);
          solver->LinSysRes.SetBlock(iPoint, Residual_Old);
        }
        }
    }
    
    delete [] Residual;
    
  }
}

void CMultiGridIntegration::Smooth_Solution(unsigned short RunTime_EqSystem, CSolver *solver, CGeometry *geometry,
                                            unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  su2double *Solution_Old, *Solution_Sum, *Solution, *Solution_i, *Solution_j;
  unsigned short iVar, iSmooth, iMarker, nneigh;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  
  const unsigned short nVar = solver->GetnVar();
  
  if (val_nSmooth > 0) {
    
    Solution = new su2double [nVar];
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Solution_Old = solver->GetNodes()->GetSolution(iPoint);
      solver->GetNodes()->SetResidual_Old(iPoint,Solution_Old);
    }
    
    /*--- Jacobi iterations ---*/
    
    for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
      solver->GetNodes()->SetResidualSumZero();
      
      /*--- Loop over Interior edges ---*/
      
      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        iPoint = geometry->edge[iEdge]->GetNode(0);
        jPoint = geometry->edge[iEdge]->GetNode(1);
        
        Solution_i = solver->GetNodes()->GetSolution(iPoint);
        Solution_j = solver->GetNodes()->GetSolution(jPoint);
        
        /*--- Accumulate nearest neighbor Residual to Res_sum for each variable ---*/
        
        solver->GetNodes()->AddResidual_Sum(iPoint,Solution_j);
        solver->GetNodes()->AddResidual_Sum(jPoint,Solution_i);
      }
      
      /*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        nneigh = geometry->node[iPoint]->GetnPoint();
        Solution_Sum = solver->GetNodes()->GetResidual_Sum(iPoint);
        Solution_Old = solver->GetNodes()->GetResidual_Old(iPoint);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] =(Solution_Old[iVar] + val_smooth_coeff*Solution_Sum[iVar])
          /(1.0 + val_smooth_coeff*su2double(nneigh));
        }
        solver->GetNodes()->SetSolution(iPoint,Solution);
      }
      
      /*--- Copy boundary values ---*/
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
        if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
            (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
          for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Solution_Old = solver->GetNodes()->GetResidual_Old(iPoint);
          solver->GetNodes()->SetSolution(iPoint,Solution_Old);
        }
        }
    }
    
    delete [] Solution;
    
  }
  
}

void CMultiGridIntegration::SetProlongated_Correction(CSolver *sol_fine, CGeometry *geo_fine, CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine;
  unsigned short iVar;
  su2double *Solution_Fine, *Residual_Fine;
  
  const unsigned short nVar = sol_fine->GetnVar();
  su2double factor = config->GetDamp_Correc_Prolong(); //pow(config->GetDamp_Correc_Prolong(), iMesh+1);
  
  su2double *Solution = new su2double [nVar];
  
  for (Point_Fine = 0; Point_Fine < geo_fine->GetnPointDomain(); Point_Fine++) {
    Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
    Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
    for (iVar = 0; iVar < nVar; iVar++) {
      /*--- Prevent a fine grid divergence due to a coarse grid divergence ---*/
      if (Residual_Fine[iVar] != Residual_Fine[iVar]) Residual_Fine[iVar] = 0.0;
      Solution[iVar] = Solution_Fine[iVar]+factor*Residual_Fine[iVar];
    }
    sol_fine->GetNodes()->SetSolution(Point_Fine,Solution);
  }
  
  /*--- MPI the new interpolated solution ---*/
  
  sol_fine->InitiateComms(geo_fine, config, SOLUTION);
  sol_fine->CompleteComms(geo_fine, config, SOLUTION);
  
  delete [] Solution;
}


void CMultiGridIntegration::SetProlongated_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iChildren;
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      sol_fine->GetNodes()->SetSolution(Point_Fine, sol_coarse->GetNodes()->GetSolution(Point_Coarse));
    }
  }
}

void CMultiGridIntegration::SetForcing_Term(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config, unsigned short iMesh) {
  unsigned long Point_Fine, Point_Coarse, iVertex;
  unsigned short iMarker, iVar, iChildren;
  su2double *Residual_Fine;
  
  const unsigned short nVar = sol_coarse->GetnVar();
  su2double factor = config->GetDamp_Res_Restric(); //pow(config->GetDamp_Res_Restric(), iMesh);
  
  su2double *Residual = new su2double[nVar];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);
    
    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += factor*Residual_Fine[iVar];
    }
    sol_coarse->GetNodes()->AddRes_TruncError(Point_Coarse,Residual);
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE    )) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetVel_ResTruncError_Zero(Point_Coarse);
      }
    }
  }
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->GetNodes()->SubtractRes_TruncError(Point_Coarse, sol_coarse->LinSysRes.GetBlock(Point_Coarse));
  }
  
  delete [] Residual;
}

void CMultiGridIntegration::SetResidual_Term(CGeometry *geometry, CSolver *solver) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    solver->LinSysRes.AddBlock(iPoint, solver->GetNodes()->GetResTruncError(iPoint));
  
}

void CMultiGridIntegration::SetRestricted_Residual(CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iVar, iChildren;
  su2double *Residual_Fine;
  
  const unsigned short nVar = sol_coarse->GetnVar();
  
  su2double *Residual = new su2double[nVar];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    sol_coarse->GetNodes()->SetRes_TruncErrorZero(Point_Coarse);
    
    for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Residual_Fine = sol_fine->LinSysRes.GetBlock(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] += Residual_Fine[iVar];
    }
    sol_coarse->GetNodes()->AddRes_TruncError(Point_Coarse,Residual);
  }
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE    )) {
      for (iVertex = 0; iVertex<geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetVel_ResTruncError_Zero(Point_Coarse);
      }
    }
  }
  
  delete [] Residual;
}

void CMultiGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iVar, iChildren, iDim;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Grid_Vel, Vector[3];
  
  const unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
  const unsigned short nVar = sol_coarse->GetnVar();
  const unsigned short nDim = geo_fine->GetnDim();
  const bool grid_movement  = config->GetGrid_Movement();
  
  su2double *Solution = new su2double[nVar];
  
  /*--- Compute coarse solution from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
      }
    }
    
    sol_coarse->GetNodes()->SetSolution(Point_Coarse,Solution);
    
  }
  
  /*--- Update the solution at the no-slip walls ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE    )) {
      
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        
        if (SolContainer_Position == FLOW_SOL) {
          
          /*--- At moving walls, set the solution based on the new density and wall velocity ---*/
          
          if (grid_movement) {
            Grid_Vel = geo_coarse->node[Point_Coarse]->GetGridVel();
            for (iDim = 0; iDim < nDim; iDim++)
              Vector[iDim] = sol_coarse->GetNodes()->GetSolution(Point_Coarse,0)*Grid_Vel[iDim];
            sol_coarse->GetNodes()->SetVelSolutionVector(Point_Coarse,Vector);
          } else {
            
            /*--- For stationary no-slip walls, set the velocity to zero. ---*/
            
            sol_coarse->GetNodes()->SetVelSolutionZero(Point_Coarse);
          }
          
        }
        
        if (SolContainer_Position == ADJFLOW_SOL) {
          sol_coarse->GetNodes()->SetVelSolutionDVector(Point_Coarse);
        }
        
      }
    }
  }
  
  /*--- MPI the new interpolated solution ---*/
  
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION);
  
  delete [] Solution;
  
}

void CMultiGridIntegration::SetRestricted_Gradient(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine,
                                                   CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iVar, iDim, iChildren;
  su2double Area_Parent, Area_Children, **Gradient_fine;
  
  const unsigned short nDim = geo_coarse->GetnDim();
  const unsigned short nVar = sol_coarse->GetnVar();
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double [nDim];
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Gradient[iVar][iDim] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Gradient_fine = sol_fine->GetNodes()->GetGradient(Point_Fine);
      
      for (iVar = 0; iVar < nVar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
          Gradient[iVar][iDim] += Gradient_fine[iVar][iDim]*Area_Children/Area_Parent;
    }
    sol_coarse->GetNodes()->SetGradient(Point_Coarse,Gradient);
  }
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;
  
}

void CMultiGridIntegration::NonDimensional_Parameters(CGeometry **geometry, CSolver ***solver_container, CNumerics ****numerics_container,
                                                      CConfig *config, unsigned short FinestMesh, unsigned short RunTime_EqSystem, 
                                                      su2double *monitor) {
    
  switch (RunTime_EqSystem) {
      
    case RUNTIME_FLOW_SYS:
      
      /*--- Calculate the inviscid and viscous forces ---*/
      
      solver_container[FinestMesh][FLOW_SOL]->Pressure_Forces(geometry[FinestMesh], config);
      solver_container[FinestMesh][FLOW_SOL]->Momentum_Forces(geometry[FinestMesh], config);
      solver_container[FinestMesh][FLOW_SOL]->Friction_Forces(geometry[FinestMesh], config);
          
      /*--- Evaluate the buffet metric if requested ---*/
      
      if(config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
          solver_container[FinestMesh][FLOW_SOL]->Buffet_Monitoring(geometry[FinestMesh], config);
      }
      
      break;
      
    case RUNTIME_ADJFLOW_SYS:
      
      /*--- Calculate the inviscid and viscous sensitivities ---*/
      
      solver_container[FinestMesh][ADJFLOW_SOL]->Inviscid_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh], numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      solver_container[FinestMesh][ADJFLOW_SOL]->Viscous_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh], numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      
      /*--- Smooth the inviscid and viscous sensitivities ---*/
      
      if (config->GetKind_SensSmooth() != NONE) solver_container[FinestMesh][ADJFLOW_SOL]->Smooth_Sensitivity(geometry[FinestMesh], solver_container[FinestMesh], numerics_container[FinestMesh][ADJFLOW_SOL][CONV_BOUND_TERM], config);
      
      break;
          
  }
  
}

CSingleGridIntegration::CSingleGridIntegration(CConfig *config) : CIntegration(config) { }

CSingleGridIntegration::~CSingleGridIntegration(void) { }

void CSingleGridIntegration::SingleGrid_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) {
  unsigned short iMesh;
  
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- Preprocessing ---*/
  
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh, 0, RunTime_EqSystem, false);
  
  /*--- Set the old solution ---*/
  
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Set_OldSolution(geometry[iZone][iInst][FinestMesh]);
  
  /*--- Time step evaluation ---*/
  
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh, config[iZone]->GetTimeIter());
  
  /*--- Space integration ---*/
  
  Space_Integration(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], numerics_container[iZone][iInst][FinestMesh][SolContainer_Position],
                    config[iZone], FinestMesh, NO_RK_ITER, RunTime_EqSystem);
  
  /*--- Time integration ---*/
  
  Time_Integration(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], NO_RK_ITER,
                   RunTime_EqSystem);
  
  /*--- Postprocessing ---*/
  
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone], FinestMesh);
  
  if (RunTime_EqSystem == RUNTIME_HEAT_SYS) {
    solver_container[iZone][iInst][FinestMesh][HEAT_SOL]->Heat_Fluxes(geometry[iZone][iInst][FinestMesh], solver_container[iZone][iInst][FinestMesh], config[iZone]);
  }
  
  /*--- If turbulence model, copy the turbulence variables to the coarse levels ---*/
  
  if (RunTime_EqSystem == RUNTIME_TURB_SYS) {
    for (iMesh = FinestMesh; iMesh < config[iZone]->GetnMGLevels(); iMesh++) {
      SetRestricted_Solution(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
      SetRestricted_EddyVisc(RunTime_EqSystem, solver_container[iZone][iInst][iMesh][SolContainer_Position], solver_container[iZone][iInst][iMesh+1][SolContainer_Position], geometry[iZone][iInst][iMesh], geometry[iZone][iInst][iMesh+1], config[iZone]);
    }
  }
}

void CSingleGridIntegration::SetRestricted_Solution(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iVar, iChildren;
  su2double Area_Parent, Area_Children, *Solution_Fine, *Solution;
  
  unsigned short nVar = sol_coarse->GetnVar();
  
  Solution = new su2double[nVar];
  
  /*--- Compute coarse solution from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      Solution_Fine = sol_fine->GetNodes()->GetSolution(Point_Fine);
      for (iVar = 0; iVar < nVar; iVar++)
        Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
    }
    
    sol_coarse->GetNodes()->SetSolution(Point_Coarse,Solution);
    
  }
  
  /*--- MPI the new interpolated solution ---*/
  
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION);
  
  delete [] Solution;
  
}

void CSingleGridIntegration::SetRestricted_EddyVisc(unsigned short RunTime_EqSystem, CSolver *sol_fine, CSolver *sol_coarse, CGeometry *geo_fine, CGeometry *geo_coarse, CConfig *config) {
  
  unsigned long iVertex, Point_Fine, Point_Coarse;
  unsigned short iMarker, iChildren;
  su2double Area_Parent, Area_Children, EddyVisc_Fine, EddyVisc;
  
  /*--- Compute coarse Eddy Viscosity from fine solution ---*/
  
  for (Point_Coarse = 0; Point_Coarse < geo_coarse->GetnPointDomain(); Point_Coarse++) {
    Area_Parent = geo_coarse->node[Point_Coarse]->GetVolume();
    
    EddyVisc = 0.0;
    
    for (iChildren = 0; iChildren < geo_coarse->node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = geo_coarse->node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geo_fine->node[Point_Fine]->GetVolume();
      EddyVisc_Fine = sol_fine->GetNodes()->GetmuT(Point_Fine);
      EddyVisc += EddyVisc_Fine*Area_Children/Area_Parent;
    }
    
    sol_coarse->GetNodes()->SetmuT(Point_Coarse,EddyVisc);
    
  }
  
  /*--- Update solution at the no slip wall boundary, only the first
   variable (nu_tilde -in SA and SA_NEG- and k -in SST-), to guarantee that the eddy viscoisty
   is zero on the surface ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX              ) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL             ) ||
        (config->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE     )) {
      for (iVertex = 0; iVertex < geo_coarse->nVertex[iMarker]; iVertex++) {
        Point_Coarse = geo_coarse->vertex[iMarker][iVertex]->GetNode();
        sol_coarse->GetNodes()->SetmuT(Point_Coarse,0.0);
      }
    }
  }

  /*--- MPI the new interpolated solution (this also includes the eddy viscosity) ---*/
    
  sol_coarse->InitiateComms(geo_coarse, config, SOLUTION_EDDY);
  sol_coarse->CompleteComms(geo_coarse, config, SOLUTION_EDDY);
  
}


CStructuralIntegration::CStructuralIntegration(CConfig *config) : CIntegration(config) { }

CStructuralIntegration::~CStructuralIntegration(void) { }

void CStructuralIntegration::Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) {

  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  /*--- Preprocessing ---*/

  solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0],
      config[iZone], numerics_container[iZone][iInst][MESH_0][SolContainer_Position], MESH_0, NO_RK_ITER, RunTime_EqSystem, false);


  /*--- Space integration ---*/

  Space_Integration_FEM(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                    config[iZone], RunTime_EqSystem);

  /*--- Time integration ---*/

  Time_Integration_FEM(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                config[iZone], RunTime_EqSystem);

  /*--- Postprocessing ---*/

  solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0],
      config[iZone], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],  MESH_0);

}

CFEM_DG_Integration::CFEM_DG_Integration(CConfig *config) : CIntegration(config) { }

CFEM_DG_Integration::~CFEM_DG_Integration(void) { }

void CFEM_DG_Integration::SingleGrid_Iteration(CGeometry ****geometry,
                                               CSolver *****solver_container,
                                               CNumerics ******numerics_container,
                                               CConfig **config,
                                               unsigned short RunTime_EqSystem,
                                               unsigned short iZone,
                                               unsigned short iInst) {

  unsigned short iMesh, iStep, iLimit = 1;
  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);
  unsigned short FinestMesh = config[iZone]->GetFinestMesh();

  /*--- For now, we assume no geometric multigrid. ---*/
  iMesh = FinestMesh;

  /*--- Check if only the Jacobian of the spatial discretization must
        be computed. If so, call the appropriate function and return. ---*/
  if (config[iZone]->GetJacobian_Spatial_Discretization_Only()) {
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->ComputeSpatialJacobian(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                         numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                                                                                         config[iZone], iMesh, RunTime_EqSystem);
    return;
  }

  /*--- Determine the number of stages in the time stepping algorithm.
        For the Runge-Kutta schemes this is the number of RK stages,
        while for ADER-DG this information is not used, because a more
        complicated algorithm must be used to facilitate time accurate
        local time stepping.  Note that we are currently hard-coding
        the classical RK4 scheme. ---*/
  bool useADER = false;
  switch (config[iZone]->GetKind_TimeIntScheme()) {
    case RUNGE_KUTTA_EXPLICIT: iLimit = config[iZone]->GetnRKStep(); break;
    case CLASSICAL_RK4_EXPLICIT: iLimit = 4; break;
    case ADER_DG: iLimit = 1; useADER = true; break;
    case EULER_EXPLICIT: case EULER_IMPLICIT: iLimit = 1; break; }

  /*--- In case an unsteady simulation is carried out, it is possible that a
        synchronization time step is specified. If so, set the boolean
        TimeSynSpecified to true, which leads to an outer loop in the
        algorithm below. ---*/
  bool TimeSyncSpecified   = false;
  const su2double TimeSync = config[iZone]->GetTime_Step()/config[iZone]->GetTime_Ref();
  if(config[iZone]->GetTime_Marching() == TIME_STEPPING &&
     config[iZone]->GetUnst_CFL()            != 0.0           &&
     TimeSync                                != 0.0) TimeSyncSpecified = true;

  /*--- Outer loop, which is only active when a synchronization time has been
        specified for an unsteady simulation. ---*/
  bool syncTimeReached = false;
  su2double timeEvolved = 0.0;
  while( !syncTimeReached ) {

    /* Compute the time step for stability. */
    solver_container[iZone][iInst][iMesh][SolContainer_Position]->SetTime_Step(geometry[iZone][iInst][iMesh],
                                                                               solver_container[iZone][iInst][iMesh],
                                                                               config[iZone], iMesh, config[iZone]->GetTimeIter());
    /* Possibly overrule the specified time step when a synchronization time was
       specified and determine whether or not the time loop must be continued.
       When TimeSyncSpecified is false, the loop is always terminated. */
    if( TimeSyncSpecified )
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->CheckTimeSynchronization(config[iZone],
                                                                                             TimeSync, timeEvolved,
                                                                                             syncTimeReached);
    else
      syncTimeReached = true;

    /*--- For ADER in combination with time accurate local time stepping, the
          space and time integration are tightly coupled and cannot be treated
          segregatedly. Therefore a different function is called for ADER to
          carry out the space and time integration. ---*/
    if( useADER ) {
      solver_container[iZone][iInst][iMesh][SolContainer_Position]->ADER_SpaceTimeIntegration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                              numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                                                                                              config[iZone], iMesh, RunTime_EqSystem);
    }
    else {

      /*--- Time and space integration can be decoupled. ---*/
      for (iStep = 0; iStep < iLimit; iStep++) {

        /*--- Preprocessing ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                    config[iZone], iMesh, iStep, RunTime_EqSystem, false);

        /*--- Space integration ---*/
        Space_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                          numerics_container[iZone][iInst][iMesh][SolContainer_Position],
                          config[iZone], iMesh, iStep, RunTime_EqSystem);

        /*--- Time integration, update solution using the old solution plus the solution increment ---*/
        Time_Integration(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                         config[iZone], iStep, RunTime_EqSystem);

        /*--- Postprocessing ---*/
        solver_container[iZone][iInst][iMesh][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][iMesh], solver_container[iZone][iInst][iMesh],
                                                                                     config[iZone], iMesh);
      }
    }
  }

  /*--- Calculate the inviscid and viscous forces ---*/
  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Pressure_Forces(geometry[iZone][iInst][iMesh], config[iZone]);

  solver_container[iZone][iInst][FinestMesh][SolContainer_Position]->Friction_Forces(geometry[iZone][iInst][iMesh], config[iZone]);

  /*--- Convergence strategy ---*/

  //Convergence_Monitoring(geometry[iZone][iInst][FinestMesh], config[iZone], Iteration, monitor, FinestMesh);
}

void CFEM_DG_Integration::Space_Integration(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics **numerics,
                                            CConfig *config, unsigned short iMesh,
                                            unsigned short iStep,
                                            unsigned short RunTime_EqSystem) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Runge-Kutta type of time integration schemes. In the first step, i.e.
        if iStep == 0, set the old solution (working solution for the DG part),
        and if needed, the new solution. ---*/
  if (iStep == 0) {
    solver_container[MainSolver]->Set_OldSolution(geometry);

    if (config->GetKind_TimeIntScheme() == CLASSICAL_RK4_EXPLICIT) {
      solver_container[MainSolver]->Set_NewSolution(geometry);
    }
  }

  /*--- Compute the spatial residual by processing the task list. ---*/
  solver_container[MainSolver]->ProcessTaskList_DG(geometry, solver_container, numerics, config, iMesh);
}

void CFEM_DG_Integration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iStep,
                                    unsigned short RunTime_EqSystem) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  /*--- Perform the time integration ---*/
  switch (config->GetKind_TimeIntScheme()) {
    case (RUNGE_KUTTA_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iStep);
      break;
    case (CLASSICAL_RK4_EXPLICIT):
      solver_container[MainSolver]->ClassicalRK4_Iteration(geometry, solver_container, config, iStep);
      break;
    default:
      SU2_MPI::Error("Time integration scheme not implemented.", CURRENT_FUNCTION);
  }
}
