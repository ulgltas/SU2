/*!
 * \file CMeshSolver.cpp
 * \brief Main subroutines to solve moving meshes using a pseudo-linear elastic approach.
 * \author Ruben Sanchez
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


#include "../../../Common/include/adt_structure.hpp"
#include "../../../Common/include/omp_structure.hpp"
#include "../../include/solvers/CMeshSolver.hpp"
#include "../../include/variables/CMeshBoundVariable.hpp"


CMeshSolver::CMeshSolver(CGeometry *geometry, CConfig *config) : CFEASolver(true) {

  /*--- Initialize some booleans that determine the kind of problem at hand. ---*/

  time_domain = config->GetTime_Domain();
  multizone = config->GetMultizone_Problem();

  /*--- Determine if the stiffness per-element is set ---*/
  switch (config->GetDeform_Stiffness_Type()) {
  case INVERSE_VOLUME:
  case SOLID_WALL_DISTANCE:
    stiffness_set = false;
    break;
  case CONSTANT_STIFFNESS:
    stiffness_set = true;
    break;
  }

  /*--- Initialize the number of spatial dimensions, length of the state
   vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/

  unsigned short iDim;
  unsigned long iPoint, iElem;

  nDim         = geometry->GetnDim();
  nVar         = geometry->GetnDim();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  nElement     = geometry->GetnElem();

  MinVolume_Ref = 0.0;
  MinVolume_Curr = 0.0;

  MaxVolume_Ref = 0.0;
  MaxVolume_Curr = 0.0;

  /*--- Initialize the node structure ---*/

  nodes = new CMeshBoundVariable(nPoint, nDim, config);
  SetBaseClassPointerToNodes();

  /*--- Set which points are vertices and allocate boundary data. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    for (iDim = 0; iDim < nDim; ++iDim)
      nodes->SetMesh_Coord(iPoint, iDim, geometry->node[iPoint]->GetCoord(iDim));

    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      long iVertex = geometry->node[iPoint]->GetVertex(iMarker);
      if (iVertex >= 0) {
        nodes->Set_isVertex(iPoint,true);
        break;
      }
    }
  }
  static_cast<CMeshBoundVariable*>(nodes)->AllocateBoundaryVariables(config);

  /*--- Initialize the element structure ---*/

  element.resize(nElement);

  /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

#ifdef HAVE_OMP
  /*--- Get the element coloring. ---*/

  const auto& coloring = geometry->GetElementColoring();

  auto nColor = coloring.getOuterSize();
  ElemColoring.resize(nColor);

  for(auto iColor = 0ul; iColor < nColor; ++iColor) {
    ElemColoring[iColor].size = coloring.getNumNonZeros(iColor);
    ElemColoring[iColor].indices = coloring.innerIdx(iColor);
  }

  ColorGroupSize = geometry->GetElementColorGroupSize();

  omp_chunk_size = computeStaticChunkSize(nPointDomain, omp_get_max_threads(), OMP_MAX_SIZE);
#endif

  /*--- Structural parameters ---*/

  E      = config->GetDeform_ElasticityMod();
  Nu     = config->GetDeform_PoissonRatio();

  Mu     = E / (2.0*(1.0 + Nu));
  Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));

  /*--- Element container structure ---*/

  if (nDim == 2) {
    for(int thread = 0; thread < omp_get_max_threads(); ++thread) {

      const int offset = thread*MAX_FE_KINDS;
      element_container[FEA_TERM][EL_TRIA+offset] = new CTRIA1();
      element_container[FEA_TERM][EL_QUAD+offset] = new CQUAD4();
    }
  }
  else {
    for(int thread = 0; thread < omp_get_max_threads(); ++thread) {

      const int offset = thread*MAX_FE_KINDS;
      element_container[FEA_TERM][EL_TETRA+offset] = new CTETRA1();
      element_container[FEA_TERM][EL_HEXA +offset] = new CHEXA8 ();
      element_container[FEA_TERM][EL_PYRAM+offset] = new CPYRAM5();
      element_container[FEA_TERM][EL_PRISM+offset] = new CPRISM6();
    }
  }

  unsigned short iVar;

  /*--- Initialize the BGS residuals in multizone problems. ---*/
  if (config->GetMultizone_Residual()){

    Residual_BGS      = new su2double[nVar]();
    Residual_Max_BGS  = new su2double[nVar]();

    /*--- Define some structures for locating max residuals ---*/

    Point_Max_BGS       = new unsigned long[nVar]();
    Point_Max_Coord_BGS = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord_BGS[iVar] = new su2double[nDim]();
    }
  }


  /*--- Allocate element properties - only the index, to allow further integration with CFEASolver on a later stage ---*/
  element_properties = new CProperty*[nElement];
  for (iElem = 0; iElem < nElement; iElem++){
    element_properties[iElem] = new CProperty(iElem);
  }

  /*--- Compute the element volumes using the reference coordinates ---*/
  SetMinMaxVolume(geometry, config, false);

  /*--- Compute the wall distance using the reference coordinates ---*/
  SetWallDistance(geometry, config);

}

void CMeshSolver::SetMinMaxVolume(CGeometry *geometry, CConfig *config, bool updated) {

  /*--- Shared reduction variables. ---*/

  unsigned long ElemCounter = 0;
  su2double MaxVolume = -1E22, MinVolume = -1E22;

  SU2_OMP_PARALLEL
  {
    /*--- Loop over the elements in the domain. ---*/

    SU2_OMP(for schedule(dynamic,omp_chunk_size) reduction(max:MaxVolume,MinVolume))
    for (unsigned long iElem = 0; iElem < nElement; iElem++) {

      int thread = omp_get_thread_num();

      int EL_KIND;
      unsigned short iNode, nNodes, iDim;

      GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

      CElement* fea_elem = element_container[FEA_TERM][EL_KIND + thread*MAX_FE_KINDS];

      /*--- For the number of nodes, we get the coordinates from
       *    the connectivity matrix and the geometry structure. ---*/

      for (iNode = 0; iNode < nNodes; iNode++) {

        auto indexNode = geometry->elem[iElem]->GetNode(iNode);

        /*--- Compute the volume with the reference or current coordinates. ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          su2double val_Coord = nodes->GetMesh_Coord(indexNode,iDim);
          if (updated)
            val_Coord += nodes->GetSolution(indexNode,iDim);

          fea_elem->SetRef_Coord(iNode, iDim, val_Coord);
        }
      }

      /*--- Compute the volume of the element (or the area in 2D cases ). ---*/

      su2double ElemVolume;
      if (nDim == 2) ElemVolume = fea_elem->ComputeArea();
      else           ElemVolume = fea_elem->ComputeVolume();

      MaxVolume = max(MaxVolume, ElemVolume);
      MinVolume = max(MinVolume, -1.0*ElemVolume);

      if (updated) element[iElem].SetCurr_Volume(ElemVolume);
      else element[iElem].SetRef_Volume(ElemVolume);

      /*--- Count distorted elements. ---*/
      if (ElemVolume <= 0.0) {
        SU2_OMP(atomic)
        ElemCounter++;
      }
    }
    MinVolume *= -1.0;

    SU2_OMP_MASTER
    {
      unsigned long ElemCounter_Local = ElemCounter;
      su2double MaxVolume_Local = MaxVolume;
      su2double MinVolume_Local = MinVolume;
      SU2_MPI::Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    SU2_OMP_BARRIER

    /*--- Volume from  0 to 1 ---*/

    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iElem = 0; iElem < nElement; iElem++) {
      if (updated) {
        su2double ElemVolume = element[iElem].GetCurr_Volume()/MaxVolume;
        element[iElem].SetCurr_Volume(ElemVolume);
      }
      else {
        su2double ElemVolume = element[iElem].GetRef_Volume()/MaxVolume;
        element[iElem].SetRef_Volume(ElemVolume);
      }
    }

  } // end SU2_OMP_PARALLEL

  /*--- Store the maximum and minimum volume. ---*/
  if (updated) {
    MaxVolume_Curr = MaxVolume;
    MinVolume_Curr = MinVolume;
  }
  else {
    MaxVolume_Ref = MaxVolume;
    MinVolume_Ref = MinVolume;
  }

  if ((ElemCounter != 0) && (rank == MASTER_NODE))
    cout <<"There are " << ElemCounter << " elements with negative volume.\n" << endl;

}

void CMeshSolver::SetWallDistance(CGeometry *geometry, CConfig *config) {

  unsigned long nVertex_SolidWall, ii, jj, iVertex, iPoint, pointID, iElem;
  unsigned short iNodes, nNodes = 0;
  unsigned short iMarker, iDim;
  su2double dist, MaxDistance_Local, MinDistance_Local;
  su2double nodeDist, ElemDist;
  int rankID;

  /*--- Initialize min and max distance ---*/

  MaxDistance = -1E22; MinDistance = 1E22;

  /*--- Compute the total number of nodes on no-slip boundaries ---*/

  nVertex_SolidWall = 0;
  for(iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if( (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {
      nVertex_SolidWall += geometry->GetnVertex(iMarker);
    }
  }

  /*--- Allocate the vectors to hold boundary node coordinates
   and its local ID. ---*/

  vector<su2double>     Coord_bound(nDim*nVertex_SolidWall);
  vector<unsigned long> PointIDs(nVertex_SolidWall);

  /*--- Retrieve and store the coordinates of the no-slip boundary nodes
   and their local point IDs. ---*/

  ii = 0; jj = 0;
  for (iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {
    if ( (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) ||
         (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)  ||
         (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {
      for (iVertex=0; iVertex<geometry->GetnVertex(iMarker); ++iVertex) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        PointIDs[jj++] = iPoint;
        for (iDim=0; iDim<nDim; ++iDim){
          Coord_bound[ii++] = nodes->GetMesh_Coord(iPoint,iDim);
        }
      }
    }
  }

  /*--- Build the ADT of the boundary nodes. ---*/

  CADTPointsOnlyClass WallADT(nDim, nVertex_SolidWall, Coord_bound.data(),
                              PointIDs.data(), true);


  /*--- Loop over all interior mesh nodes and compute the distances to each
   of the no-slip boundary nodes. Store the minimum distance to the wall
   for each interior mesh node. ---*/

  if( WallADT.IsEmpty() ) {

    /*--- No solid wall boundary nodes in the entire mesh.
     Set the wall distance to zero for all nodes. ---*/

    for (iPoint=0; iPoint<geometry->GetnPoint(); ++iPoint)
      geometry->node[iPoint]->SetWall_Distance(0.0);
  }
  else {

    /*--- Solid wall boundary nodes are present. Compute the wall
     distance for all nodes. ---*/

    for(iPoint=0; iPoint< nPoint; ++iPoint) {

      WallADT.DetermineNearestNode(nodes->GetMesh_Coord(iPoint), dist,
                                   pointID, rankID);
      nodes->SetWallDistance(iPoint,dist);

      MaxDistance = max(MaxDistance, dist);

      /*--- To discard points on the surface we use > EPS ---*/

      if (sqrt(dist) > EPS)  MinDistance = min(MinDistance, dist);

    }

    MaxDistance_Local = MaxDistance; MaxDistance = 0.0;
    MinDistance_Local = MinDistance; MinDistance = 0.0;

    SU2_MPI::Allreduce(&MaxDistance_Local, &MaxDistance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MinDistance_Local, &MinDistance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  }

  /*--- Normalize distance from 0 to 1 ---*/
  for (iPoint=0; iPoint < nPoint; ++iPoint) {
    nodeDist = nodes->GetWallDistance(iPoint)/MaxDistance;
    nodes->SetWallDistance(iPoint,nodeDist);
  }

  /*--- Compute the element distances ---*/
  for (iElem = 0; iElem < nElement; iElem++) {

    int EL_KIND;
    GetElemKindAndNumNodes(geometry->elem[iElem]->GetVTK_Type(), EL_KIND, nNodes);

    /*--- Average the distance of the nodes in the element ---*/

    ElemDist = 0.0;
    for (iNodes = 0; iNodes < nNodes; iNodes++){
      iPoint = geometry->elem[iElem]->GetNode(iNodes);
      ElemDist += nodes->GetWallDistance(iPoint);
    }
    ElemDist = ElemDist/su2double(nNodes);

    element[iElem].SetWallDistance(ElemDist);

  }

}

void CMeshSolver::SetMesh_Stiffness(CGeometry **geometry, CNumerics **numerics, CConfig *config){

  unsigned long iElem;

  if (!stiffness_set) {
    for (iElem = 0; iElem < nElement; iElem++) {

      switch (config->GetDeform_Stiffness_Type()) {
      /*--- Stiffness inverse of the volume of the element ---*/
      case INVERSE_VOLUME: E = 1.0 / element[iElem].GetRef_Volume();  break;
      /*--- Stiffness inverse of the distance of the element to the closest wall ---*/
      case SOLID_WALL_DISTANCE: E = 1.0 / element[iElem].GetWallDistance(); break;
      }

      /*--- Set the element elastic properties in the numerics container ---*/
      numerics[FEA_TERM]->SetMeshElasticProperties(iElem, E);

    }

    stiffness_set = true;
  }

}

void CMeshSolver::DeformMesh(CGeometry **geometry, CNumerics **numerics, CConfig *config){

  if (multizone) nodes->Set_BGSSolution_k();

  /*--- Initialize sparse matrix ---*/
  Jacobian.SetValZero();

  /*--- Compute the stiffness matrix. ---*/
  Compute_StiffMatrix(geometry[MESH_0], numerics, config);

  /*--- Initialize vectors and clean residual ---*/
  LinSysSol.SetValZero();
  LinSysRes.SetValZero();

  /*--- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   Hence we still need a communication of the transformed coordinates, otherwise periodicity
   is not maintained. ---*/
  geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
  geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);

  /*--- In the same way, communicate the displacements in the solver to make sure the halo
   nodes receive the correct value of the displacement. ---*/
  InitiateComms(geometry[MESH_0], config, SOLUTION);
  CompleteComms(geometry[MESH_0], config, SOLUTION);

  InitiateComms(geometry[MESH_0], config, MESH_DISPLACEMENTS);
  CompleteComms(geometry[MESH_0], config, MESH_DISPLACEMENTS);

  /*--- Impose boundary conditions (all of them are ESSENTIAL BC's - displacements). ---*/
  SetBoundaryDisplacements(geometry[MESH_0], numerics[FEA_TERM], config);

  /*--- Solve the linear system. ---*/
  Solve_System(geometry[MESH_0], config);

  /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/
  UpdateGridCoord(geometry[MESH_0], config);

  /*--- Update the dual grid. ---*/
  UpdateDualGrid(geometry[MESH_0], config);

  /*--- Check for failed deformation (negative volumes). ---*/
  /*--- In order to do this, we recompute the minimum and maximum area/volume for the mesh using the current coordinates. ---*/
  SetMinMaxVolume(geometry[MESH_0], config, true);

  /*--- The Grid Velocity is only computed if the problem is time domain ---*/
  if (time_domain) ComputeGridVelocity(geometry[MESH_0], config);

  /*--- Update the multigrid structure. ---*/
  UpdateMultiGrid(geometry, config);

}

void CMeshSolver::UpdateGridCoord(CGeometry *geometry, CConfig *config){

  /*--- Update the grid coordinates using the solution of the linear system ---*/

  /*--- LinSysSol contains the absolute x, y, z displacements. ---*/
  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      auto total_index = iPoint*nDim + iDim;
      /*--- Retrieve the displacement from the solution of the linear system ---*/
      su2double val_disp = LinSysSol[total_index];
      /*--- Store the displacement of the mesh node ---*/
      nodes->SetSolution(iPoint, iDim, val_disp);
      /*--- Compute the current coordinate as Mesh_Coord + Displacement ---*/
      su2double val_coord = nodes->GetMesh_Coord(iPoint,iDim) + val_disp;
      /*--- Update the geometry container ---*/
      geometry->node[iPoint]->SetCoord(iDim, val_coord);
    }
  }

  /*--- LinSysSol contains the non-transformed displacements in the periodic halo cells.
   Hence we still need a communication of the transformed coordinates, otherwise periodicity
   is not maintained. ---*/
  geometry->InitiateComms(geometry, config, COORDINATES);
  geometry->CompleteComms(geometry, config, COORDINATES);

  /*--- In the same way, communicate the displacements in the solver to make sure the halo
   nodes receive the correct value of the displacement. ---*/
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

}

void CMeshSolver::UpdateDualGrid(CGeometry *geometry, CConfig *config){

  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/

  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, UPDATE);
  geometry->SetBoundControlVolume(config, UPDATE);
  geometry->SetMaxLength(config);

}

void CMeshSolver::ComputeGridVelocity(CGeometry *geometry, CConfig *config){

  /*--- Compute the velocity of each node in the domain of the current rank
   (halo nodes are not computed as the grid velocity is later communicated). ---*/

  SU2_OMP(parallel for schedule(static,omp_chunk_size))
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Coordinates of the current point at n+1, n, & n-1 time levels. ---*/

    const su2double* Disp_nM1 = nodes->GetSolution_time_n1(iPoint);
    const su2double* Disp_n   = nodes->GetSolution_time_n(iPoint);
    const su2double* Disp_nP1 = nodes->GetSolution(iPoint);

    /*--- Unsteady time step ---*/

    su2double TimeStep = config->GetDelta_UnstTimeND();

    /*--- Compute mesh velocity with 1st or 2nd-order approximation. ---*/

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {

      su2double GridVel = 0.0;

      if (config->GetTime_Marching() == DT_STEPPING_1ST)
        GridVel = ( Disp_nP1[iDim] - Disp_n[iDim] ) / TimeStep;
      if (config->GetTime_Marching() == DT_STEPPING_2ND)
        GridVel = ( 3.0*Disp_nP1[iDim] - 4.0*Disp_n[iDim]
                    +  1.0*Disp_nM1[iDim] ) / (2.0*TimeStep);

      /*--- Store grid velocity for this point ---*/

      geometry->node[iPoint]->SetGridVel(iDim, GridVel);

    }
  }

  /*--- The velocity was computed for nPointDomain, now we communicate it. ---*/
  geometry->InitiateComms(geometry, config, GRID_VELOCITY);
  geometry->CompleteComms(geometry, config, GRID_VELOCITY);

}

void CMeshSolver::UpdateMultiGrid(CGeometry **geometry, CConfig *config){

  unsigned short iMGfine, iMGlevel, nMGlevel = config->GetnMGLevels();

  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels
   when the problem is solved in unsteady conditions. ---*/

  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel-1;
    geometry[iMGlevel]->SetControlVolume(config, geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetBoundControlVolume(config, geometry[iMGfine],UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (time_domain)
      geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine], config);
  }

}

void CMeshSolver::SetBoundaryDisplacements(CGeometry *geometry, CNumerics *numerics, CConfig *config){

  unsigned short iMarker;

  /*--- Impose zero displacements of all non-moving surfaces (also at nodes in multiple moving/non-moving boundaries). ---*/
  /*--- Exceptions: symmetry plane, the receive boundaries and periodic boundaries should get a different treatment. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_Deform_Mesh(iMarker) == NO) &&
        (config->GetMarker_All_KindBC(iMarker) != SYMMETRY_PLANE) &&
        (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {

      BC_Clamped(geometry, numerics, config, iMarker);
    }
  }

  /*--- Symmetry plane is, for now, clamped. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_Deform_Mesh(iMarker) == NO) &&
        (config->GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE)) {

      BC_Clamped(geometry, numerics, config, iMarker);
    }
  }

  /*--- Impose displacement boundary conditions. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Deform_Mesh(iMarker) == YES) {

      BC_Deforming(geometry, numerics, config, iMarker);
    }
  }

}

void CMeshSolver::SetDualTime_Mesh(void){

  nodes->Set_Solution_time_n1();
  nodes->Set_Solution_time_n();
}

void CMeshSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  string filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  unsigned long iPoint_Global, counter = 0;

  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local >= 0) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      auto index = counter*Restart_Vars[1];

      for (unsigned short iDim = 0; iDim < nDim; iDim++){
        /*--- Update the coordinates of the mesh ---*/
        su2double curr_coord = Restart_Data[index+iDim];
        /// TODO: "Double deformation" in multizone adjoint if this is set here?
        ///       In any case it should not be needed as deformation is called before other solvers
        ///geometry[MESH_0]->node[iPoint_Local]->SetCoord(iDim, curr_coord);

        /*--- Store the displacements computed as the current coordinates
         minus the coordinates of the reference mesh file ---*/
        su2double displ = curr_coord - nodes->GetMesh_Coord(iPoint_Local, iDim);
        nodes->SetSolution(iPoint_Local, iDim, displ);
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/

  if (counter != nPointDomain) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Communicate the loaded displacements. ---*/
  solver[MESH_0][MESH_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][MESH_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  /*--- Communicate the new coordinates at the halos ---*/
  geometry[MESH_0]->InitiateComms(geometry[MESH_0], config, COORDINATES);
  geometry[MESH_0]->CompleteComms(geometry[MESH_0], config, COORDINATES);

  /*--- Recompute the edges and dual mesh control volumes in the
   domain and on the boundaries. ---*/
  UpdateDualGrid(geometry[MESH_0], config);

  /*--- For time-domain problems, we need to compute the grid velocities ---*/
  if (time_domain){
    /*--- Update the old geometry (coordinates n and n-1) ---*/
    Restart_OldGeometry(geometry[MESH_0], config);
    /*--- Once Displacement_n and Displacement_n1 are filled,
     we can compute the Grid Velocity ---*/
    ComputeGridVelocity(geometry[MESH_0], config);
  }

  /*--- Update the multigrid structure after setting up the finest grid,
   including computing the grid velocities on the coarser levels
   when the problem is unsteady. ---*/
  UpdateMultiGrid(geometry, config);

  /*--- Store the boundary displacements at the Bound_Disp variable. ---*/

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_Deform_Mesh(iMarker) == YES) {

      for (unsigned long iVertex = 0; iVertex < geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

        /*--- Get node index. ---*/
        auto iNode = geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();

        /*--- Set boundary solution. ---*/
        nodes->SetBound_Disp(iNode, nodes->GetSolution(iNode));
      }
    }
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  delete [] Restart_Vars; Restart_Vars = nullptr;
  delete [] Restart_Data; Restart_Data = nullptr;

}

void CMeshSolver::Restart_OldGeometry(CGeometry *geometry, CConfig *config) {

  /*--- This function is intended for dual time simulations ---*/

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  string filename = config->GetSolution_FileName();

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone, "");

  /*--- Determine how many files need to be read. ---*/

  unsigned short nSteps = (config->GetTime_Marching() == DT_STEPPING_2ND) ? 2 : 1;

  for(unsigned short iStep = 1; iStep <= nSteps; ++iStep) {

    unsigned short CommType = (iStep == 1) ? SOLUTION_TIME_N : SOLUTION_TIME_N1;

    /*--- Modify file name for an unsteady restart ---*/
    int Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-iStep;
    string filename_n = config->GetUnsteady_FileName(filename, Unst_RestartIter, "");

    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry, config, filename_n);
    } else {
      Read_SU2_Restart_ASCII(geometry, config, filename_n);
    }

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long iPoint_Global, counter = 0;

    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++) {

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      auto iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local >= 0) {

        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        auto index = counter*Restart_Vars[1];

        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          su2double curr_coord = Restart_Data[index+iDim];
          su2double displ = curr_coord - nodes->GetMesh_Coord(iPoint_Local,iDim);

          if(iStep==1)
            nodes->Set_Solution_time_n(iPoint_Local, iDim, displ);
          else
            nodes->Set_Solution_time_n1(iPoint_Local, iDim, displ);
        }

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }

    }

    /*--- Detect a wrong solution file. ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + filename_n + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Delete the class memory that is used to load the restart. ---*/

    delete [] Restart_Vars; Restart_Vars = nullptr;
    delete [] Restart_Data; Restart_Data = nullptr;

    InitiateComms(geometry, config, CommType);
    CompleteComms(geometry, config, CommType);

  } // iStep

}
