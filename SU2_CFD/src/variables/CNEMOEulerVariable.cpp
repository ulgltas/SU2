/*!
 * \file CNEMOEulerVariable.cpp
 * \brief Definition of the solution fields.
 * \author C. Garbacz, W. Maier, S.R. Copeland
 * \version 7.1.1 "Blackbird"
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

#include "../../include/variables/CNEMOEulerVariable.hpp"
#include <math.h>

CNEMOEulerVariable::CNEMOEulerVariable(su2double val_pressure,
                                       const su2double *val_massfrac,
                                       su2double *val_mach,
                                       su2double val_temperature,
                                       su2double val_temperature_ve,
                                       unsigned long npoint,
                                       unsigned long ndim,
                                       unsigned long nvar,
                                       unsigned long nvarprim,
                                       unsigned long nvarprimgrad,
                                       CConfig *config,
                                       CNEMOGas *fluidmodel) : CVariable(npoint,
                                                                         ndim,
                                                                         nvar,
                                                                         config ),
                                       Gradient_Reconstruction(config->GetReconstructionGradientRequired() ? Gradient_Aux : Gradient_Primitive) {

  vector<su2double> energies;
  unsigned short iDim, iSpecies;
  su2double soundspeed, sqvel, rho;

  /*--- Setting variable amounts ---*/
  nDim            = ndim;
  nPrimVar        = nvarprim;
  nPrimVarGrad    = nvarprimgrad;

  nSpecies        = config->GetnSpecies();
  RHOS_INDEX      = 0;
  T_INDEX         = nSpecies;
  TVE_INDEX       = nSpecies+1;
  VEL_INDEX       = nSpecies+2;
  P_INDEX         = nSpecies+nDim+2;
  RHO_INDEX       = nSpecies+nDim+3;
  H_INDEX         = nSpecies+nDim+4;
  A_INDEX         = nSpecies+nDim+5;
  RHOCVTR_INDEX   = nSpecies+nDim+6;
  RHOCVVE_INDEX   = nSpecies+nDim+7;
  LAM_VISC_INDEX  = nSpecies+nDim+8;
  EDDY_VISC_INDEX = nSpecies+nDim+9;

  /*--- Set monoatomic flag ---*/
  if (config->GetMonoatomic()) {
    monoatomic = true;
    Tve_Freestream = config->GetTemperature_ve_FreeStream();
  }

  /*--- Allocate & initialize residual vectors ---*/
  Res_TruncError.resize(nPoint,nVar) = su2double(0.0);

  /*--- Size Grad_AuxVar for axiysmmetric ---*/
  if (config->GetAxisymmetric()){
    nAuxVar = 3;
    Grad_AuxVar.resize(nPoint,nAuxVar,nDim,0.0);
    AuxVar.resize(nPoint,nAuxVar) = su2double(0.0);
  }

  /*--- Only for residual smoothing (multigrid) ---*/
  for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (config->GetMG_CorrecSmooth(iMesh) > 0) {
      Residual_Sum.resize(nPoint,nVar);
      Residual_Old.resize(nPoint,nVar);
      break;
    }
  }

  /*--- Allocate undivided laplacian (centered) and limiter (upwind)---*/
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED)
    Undivided_Laplacian.resize(nPoint,nVar);

  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/
  Limiter.resize(nPoint,nVar) = su2double(0.0);
  Limiter_Primitive.resize(nPoint,nPrimVarGrad) = su2double(0.0);

  Solution_Max.resize(nPoint,nPrimVarGrad) = su2double(0.0);
  Solution_Min.resize(nPoint,nPrimVarGrad) = su2double(0.0);

  /*--- Primitive and secondary variables ---*/
  Primitive.resize(nPoint,nPrimVar) = su2double(0.0);
  Primitive_Aux.resize(nPoint,nPrimVar) = su2double(0.0);
  Secondary.resize(nPoint,nPrimVar) = su2double(0.0);

  dPdU.resize(nPoint, nVar)      = su2double(0.0);
  dTdU.resize(nPoint, nVar)      = su2double(0.0);
  dTvedU.resize(nPoint, nVar)    = su2double(0.0);
  Cvves.resize(nPoint, nSpecies) = su2double(0.0);
  eves.resize(nPoint, nSpecies)  = su2double(0.0);
  Gamma.resize(nPoint)           = su2double(0.0);

  /*--- Compressible flow, gradients primitive variables ---*/
  Gradient_Primitive.resize(nPoint,nPrimVarGrad,nDim,0.0);
  Gradient.resize(nPoint,nVar,nDim,0.0);

  if (config->GetReconstructionGradientRequired()) {
    Gradient_Aux.resize(nPoint,nPrimVarGrad,nDim,0.0);
  }

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  Velocity2.resize(nPoint) = su2double(0.0);
  Max_Lambda_Inv.resize(nPoint) = su2double(0.0);
  Delta_Time.resize(nPoint) = su2double(0.0);
  Lambda.resize(nPoint) = su2double(0.0);
  Sensor.resize(nPoint) = su2double(0.0);

  /* Non-physical point (first-order) initialization. */
  Non_Physical.resize(nPoint) = false;
  Non_Physical_Counter.resize(nPoint) = 0;

  /* Under-relaxation parameter. */
  UnderRelaxation.resize(nPoint) = su2double(1.0);
  LocalCFL.resize(nPoint) = su2double(0.0);

  /*--- Loop over all points --*/
  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint){

    /*--- Reset velocity^2 [m2/s2] to zero ---*/
    sqvel = 0.0;

    /*--- Set mixture state ---*/
    fluidmodel->SetTDStatePTTv(val_pressure, val_massfrac, val_temperature, val_temperature_ve);

    /*--- Compute necessary quantities ---*/
    rho = fluidmodel->GetDensity();
    soundspeed = fluidmodel->ComputeSoundSpeed();
    for (iDim = 0; iDim < nDim; iDim++){
      sqvel += val_mach[iDim]*soundspeed * val_mach[iDim]*soundspeed;
    }
    energies = fluidmodel->ComputeMixtureEnergies();

    /*--- Initialize Solution & Solution_Old vectors ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Solution(iPoint,iSpecies)     = rho*val_massfrac[iSpecies];
    for (iDim = 0; iDim < nDim; iDim++)
      Solution(iPoint,nSpecies+iDim)     = rho*val_mach[iDim]*soundspeed;

    Solution(iPoint,nSpecies+nDim)       = rho*(energies[0]+0.5*sqvel);
    Solution(iPoint,nSpecies+nDim+1)     = rho*(energies[1]);

    Solution_Old = Solution;

    /*--- Assign primitive variables ---*/
    Primitive(iPoint,T_INDEX)   = val_temperature;
    Primitive(iPoint,TVE_INDEX) = val_temperature_ve;
    Primitive(iPoint,P_INDEX)   = val_pressure;
  }
}

void CNEMOEulerVariable::SetVelocity2(unsigned long iPoint) {

  unsigned short iDim;

  Velocity2(iPoint) = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Primitive(iPoint,VEL_INDEX+iDim) = Solution(iPoint,nSpecies+iDim) / Primitive(iPoint,RHO_INDEX);
    Velocity2(iPoint) +=  Solution(iPoint,nSpecies+iDim)*Solution(iPoint,nSpecies+iDim)
        / (Primitive(iPoint,RHO_INDEX)*Primitive(iPoint,RHO_INDEX));
  }
}

bool CNEMOEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  bool nonPhys;
  unsigned short iVar;

  fluidmodel = static_cast<CNEMOGas*>(FluidModel);

  /*--- Convert conserved to primitive variables ---*/
  nonPhys = Cons2PrimVar(Solution[iPoint], Primitive[iPoint],
                         dPdU[iPoint], dTdU[iPoint], dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);

  /*--- Reset solution to previous one, if nonphys ---*/
  if (nonPhys) {
    for (iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);
  }

  /*--- Set additional point quantaties ---*/
  Gamma(iPoint) = fluidmodel->ComputeGamma();

  SetVelocity2(iPoint);

  return nonPhys;
}

bool CNEMOEulerVariable::Cons2PrimVar(su2double *U, su2double *V,
                                      su2double *val_dPdU, su2double *val_dTdU,
                                      su2double *val_dTvedU, su2double *val_eves,
                                      su2double *val_Cvves) {

  unsigned short iDim, iSpecies;
  su2double Tmin, Tmax, Tvemin, Tvemax;
  vector<su2double> rhos;

  rhos.resize(nSpecies,0.0);

  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  // V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T

  /*--- Set booleans ---*/
  bool nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0; Tmax   = 8E4;
  Tvemin = 50.0; Tvemax = 8E4;

  /*--- Rename variables for convenience ---*/
  su2double rhoE   = U[nSpecies+nDim];     // Density * energy [J/m3]
  su2double rhoEve = U[nSpecies+nDim+1];   // Density * energy_ve [J/m3]

  /*--- Assign species & mixture density ---*/
  // Note: if any species densities are < 0, these values are re-assigned
  //       in the primitive AND conserved vectors to ensure positive density
  V[RHO_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (U[iSpecies] < 0.0) {
      U[iSpecies]            = 1E-20;
      V[RHOS_INDEX+iSpecies] = 1E-20;
      rhos[iSpecies]         = 1E-20;
      //nonPhys                = true;
    } else {
      V[RHOS_INDEX+iSpecies] = U[iSpecies];
      rhos[iSpecies]         = U[iSpecies];
    }
    V[RHO_INDEX]            += U[iSpecies];
  }

  // Rename for convenience
  su2double rho = V[RHO_INDEX];

  /*--- Assign velocity^2 ---*/
  su2double sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    V[VEL_INDEX+iDim] = U[nSpecies+iDim]/V[RHO_INDEX];
    sqvel            += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  }

  /*--- Assign temperatures ---*/
  const auto& T = fluidmodel->ComputeTemperatures(rhos, rhoE, rhoEve, 0.5*rho*sqvel);

  /*--- Temperatures ---*/
  V[T_INDEX]   = T[0];
  V[TVE_INDEX] = T[1];

  // Determine if the temperature lies within the acceptable range
  //TODO: fIX THIS
  if (V[T_INDEX] == Tmin) {
    nonPhys = true;
  } else if (V[T_INDEX] == Tmax){
    nonPhys = true;
  }

  /*--- Vibrational-Electronic Temperature ---*/
  vector<su2double> eves_min = fluidmodel->ComputeSpeciesEve(Tvemin);
  vector<su2double> eves_max = fluidmodel->ComputeSpeciesEve(Tvemax);

  // Check for non-physical solutions
  if (!monoatomic){
    su2double rhoEve_min = 0.0;
    su2double rhoEve_max = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      rhoEve_min += U[iSpecies] * eves_min[iSpecies];
      rhoEve_max += U[iSpecies] * eves_max[iSpecies];
    }

    if (rhoEve < rhoEve_min) {

      nonPhys      = true;
      V[TVE_INDEX] = Tvemin;
      U[nSpecies+nDim+1] = rhoEve_min;
    } else if (rhoEve > rhoEve_max) {
      nonPhys      = true;
      V[TVE_INDEX] = Tvemax;
      U[nSpecies+nDim+1] = rhoEve_max;
    }
  } else {
    //TODO: can e-modes/vibe modes be active?
    V[TVE_INDEX] = Tve_Freestream;
  }

  // Determine other properties of the mixture at the current state
  fluidmodel->SetTDStateRhosTTv(rhos, V[T_INDEX], V[TVE_INDEX]);
  const auto& cvves = fluidmodel->ComputeSpeciesCvVibEle();
  vector<su2double> eves  = fluidmodel->ComputeSpeciesEve(V[TVE_INDEX]);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    val_eves[iSpecies]  = eves[iSpecies];
    val_Cvves[iSpecies] = cvves[iSpecies];
  }

  su2double rhoCvtr = fluidmodel->ComputerhoCvtr();
  su2double rhoCvve = fluidmodel->ComputerhoCvve();

  V[RHOCVTR_INDEX] = rhoCvtr;
  V[RHOCVVE_INDEX] = rhoCvve;

  /*--- Pressure ---*/
  V[P_INDEX] = fluidmodel->ComputePressure();

  if (V[P_INDEX] < 0.0) {
    V[P_INDEX] = 1E-20;
    nonPhys = true;
  }

  /*--- Partial derivatives of pressure and temperature ---*/
  fluidmodel->ComputedPdU  (V, eves, val_dPdU  );
  fluidmodel->ComputedTdU  (V, val_dTdU );
  fluidmodel->ComputedTvedU(V, eves, val_dTvedU);

  /*--- Sound speed ---*/
  V[A_INDEX] = fluidmodel->ComputeSoundSpeed();

  /*--- Enthalpy ---*/
  V[H_INDEX] = (U[nSpecies+nDim] + V[P_INDEX])/V[RHO_INDEX];

  return nonPhys;
}

void CNEMOEulerVariable::SetSolution_New() { Solution_New = Solution; }
