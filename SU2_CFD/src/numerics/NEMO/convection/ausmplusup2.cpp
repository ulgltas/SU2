/*!
 * \file ausmplusup2.cpp
 * \brief Implementations of the AUSM-family of schemes - AUSM+UP2.
 * \author W. Maier, A. Sachedeva, C. Garbacz
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

#include "../../../../include/numerics/NEMO/convection/ausmplusup2.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSMPLUSUP2_NEMO::CUpwAUSMPLUSUP2_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                                           unsigned short val_nPrimVar,
                                           unsigned short val_nPrimVarGrad,
                                           CConfig *config): CNEMONumerics (val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config){

  /*--- Define useful constants ---*/
  Kp       = 0.25;
  sigma    = 1.0;

  /*--- Allocate data structures ---*/
  FcL    = new su2double [nVar];
  FcR    = new su2double [nVar];
  rhos_i = new su2double [nSpecies];
  rhos_j = new su2double [nSpecies];
  u_i    = new su2double [nDim];
  u_j    = new su2double [nDim];

  Flux   = new su2double[nVar];
}

CUpwAUSMPLUSUP2_NEMO::~CUpwAUSMPLUSUP2_NEMO(void) {

  delete [] FcL;
  delete [] FcR;
  delete [] rhos_i;
  delete [] rhos_j;
  delete [] u_i;
  delete [] u_j;

  delete [] Flux;
}

CNumerics::ResidualType<> CUpwAUSMPLUSUP2_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iDim, iVar, iSpecies;
  su2double rho_i, rho_j,
  e_ve_i, e_ve_j, sq_veli, sq_velj;

  /*--- Face area ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  Minf  = config->GetMach();

  /*--- Extracting primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }

  sq_veli = 0.0; sq_velj = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim]  = V_i[VEL_INDEX+iDim];
    u_j[iDim]  = V_j[VEL_INDEX+iDim];
    sq_veli   += u_i[iDim]*u_i[iDim];
    sq_velj   += u_j[iDim]*u_j[iDim];
  }

  P_i   = V_i[P_INDEX];   P_j   = V_j[P_INDEX];
  h_i   = V_i[H_INDEX];   h_j   = V_j[H_INDEX];
  a_i   = V_i[A_INDEX];   a_j   = V_j[A_INDEX];
  rho_i = V_i[RHO_INDEX]; rho_j = V_j[RHO_INDEX];
  
  rhoCvtr_i = V_i[RHOCVTR_INDEX]; rhoCvtr_j = V_j[RHOCVTR_INDEX];
  rhoCvve_i = V_i[RHOCVVE_INDEX]; rhoCvve_j = V_j[RHOCVVE_INDEX];

  e_ve_i = 0.0; e_ve_j = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    e_ve_i += (rhos_i[iSpecies]*eve_i[iSpecies])/rho_i;
    e_ve_j += (rhos_j[iSpecies]*eve_j[iSpecies])/rho_j;
  }

  /*--- Projected velocities ---*/
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }

  /*--- Compute C*  ---*/
  CstarL = sqrt(2.0*(Gamma_i-1.0)/(Gamma_i+1.0)*h_i);
  CstarR = sqrt(2.0*(Gamma_j-1.0)/(Gamma_j+1.0)*h_j);

  /*--- Compute C^ ---*/
  ChatL = CstarL*CstarL/max(CstarL,ProjVel_i);
  ChatR = CstarR*CstarR/max(CstarR,-ProjVel_j);

  /*--- Interface speed of sound ---*/
  aF = min(ChatL,ChatR);

  mL  = ProjVel_i/aF;
  mR  = ProjVel_j/aF;

  rhoF = 0.5*(rho_i+rho_j);
  MFsq = 0.5*(mL*mL+mR*mR);

  param1 = max(MFsq, Minf*Minf);
  Mrefsq = (min(1.0, param1));
  fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  beta = 1.0/8.0;

  /*--- Pressure diffusion term ---*/
  Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(P_j-P_i)/(rhoF*aF*aF);

  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0)+beta*(mL*mL-1.0)*(mL*mL-1.0);
  else                 mLP = 0.5*(mL+fabs(mL));

  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0)-beta*(mR*mR-1.0)*(mR*mR-1.0);
  else                 mRM = 0.5*(mR-fabs(mR));

  mF = mLP + mRM + Mp;

  if (fabs(mL) <= 1.0) pLP = (0.25*(mL+1.0)*(mL+1.0)*(2.0-mL)+alpha*mL*(mL*mL-1.0)*(mL*mL-1.0));
  else                 pLP = 0.5*(mL+fabs(mL))/mL;

  if (fabs(mR) <= 1.0) pRM = (0.25*(mR-1.0)*(mR-1.0)*(2.0+mR)-alpha*mR*(mR*mR-1.0)*(mR*mR-1.0));
  else                 pRM = 0.5*(mR-fabs(mR))/mR;

  /*... Modified pressure flux ...*/
  //Use this definition
  pFi = sqrt(0.5*(sq_veli+sq_velj))*(pLP+pRM-1.0)*0.5*(rho_j+rho_i)*aF;
  pF  = 0.5*(P_j+P_i)+0.5*(pLP-pRM)*(P_i-P_j)+pFi;

  Phi = fabs(mF);

  mfP=0.5*(mF+Phi);
  mfM=0.5*(mF-Phi);

  /*--- Assign left & right covective fluxes ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    FcL[iSpecies] = rhos_i[iSpecies]*aF;
    FcR[iSpecies] = rhos_j[iSpecies]*aF;
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    FcL[nSpecies+iDim] = rho_i*aF*u_i[iDim];
    FcR[nSpecies+iDim] = rho_j*aF*u_j[iDim];
  }
  FcL[nSpecies+nDim]   = rho_i*aF*h_i;
  FcR[nSpecies+nDim]   = rho_j*aF*h_j;
  FcL[nSpecies+nDim+1] = rho_i*aF*e_ve_i;
  FcR[nSpecies+nDim+1] = rho_j*aF*e_ve_j;

  /*--- Compute numerical flux ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = (mfP*FcL[iVar]+mfM*FcR[iVar])*Area;

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[nSpecies+iDim] += pF*UnitNormal[iDim]*Area;

  /*--- Roe's Jacobian -> checking if there is an improvement over NEMO AUSM Jacobian---*/
  //if (implicit){
  //
  //  /*--- Compute Roe Variables ---*/
  //  R    = sqrt(abs(V_j[RHO_INDEX]/V_i[RHO_INDEX]));
  //  for (iVar = 0; iVar < nVar; iVar++)
  //    RoeU[iVar] = (R*U_j[iVar] + U_i[iVar])/(R+1);
  //  for (iVar = 0; iVar < nPrimVar; iVar++)
  //    RoeV[iVar] = (R*V_j[iVar] + V_i[iVar])/(R+1);

  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    RoeEve[iSpecies] = var->CalcEve(config, RoeV[TVE_INDEX], iSpecies);

  /*--- Calculate derivatives of pressure ---*/
  //    var->CalcdPdU(RoeV, RoeEve, config, RoedPdU);

  /*--- Calculate dual grid tangent vectors for P & invP ---*/
  //    CreateBasis(UnitNormal);

  /*--- Compute projected P, invP, and Lambda ---*/
  //  GetPMatrix(RoeU, RoeV, RoedPdU, UnitNormal, l, m, P_Tensor);
  //  GetPMatrix_inv(RoeU, RoeV, RoedPdU, UnitNormal, l, m, invP_Tensor);

  //  RoeSoundSpeed = sqrt((1.0+RoedPdU[nSpecies+nDim])*
  //      RoeV[P_INDEX]/RoeV[RHO_INDEX]);

  /*--- Compute projected velocities ---*/
  //  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  //  for (iDim = 0; iDim < nDim; iDim++) {
  //    ProjVelocity   += RoeV[VEL_INDEX+iDim] * UnitNormal[iDim];
  //    ProjVelocity_i += V_i[VEL_INDEX+iDim]  * UnitNormal[iDim];
  //    ProjVelocity_j += V_j[VEL_INDEX+iDim]  * UnitNormal[iDim];
  //  }

  /*--- Calculate eigenvalues ---*/
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    Lambda[iSpecies] = ProjVelocity;
  //  for (iDim = 0; iDim < nDim-1; iDim++)
  //    Lambda[nSpecies+iDim] = ProjVelocity;
  //  Lambda[nSpecies+nDim-1] = ProjVelocity + RoeSoundSpeed;
  //  Lambda[nSpecies+nDim]   = ProjVelocity - RoeSoundSpeed;
  //  Lambda[nSpecies+nDim+1] = ProjVelocity;

  /*--- Harten and Hyman (1983) entropy correction ---*/
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    Epsilon[iSpecies] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
  //                                         ProjVelocity_j-Lambda[iDim] ));
  //  for (iDim = 0; iDim < nDim-1; iDim++)
  //    Epsilon[nSpecies+iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
  //                                              ProjVelocity_j-Lambda[iDim] ));
  //  Epsilon[nSpecies+nDim-1] = 4.0*max(0.0, max(Lambda[nSpecies+nDim-1]-(ProjVelocity_i+V_i[A_INDEX]),
  //                                     (ProjVelocity_j+V_j[A_INDEX])-Lambda[nSpecies+nDim-1]));
  //  Epsilon[nSpecies+nDim]   = 4.0*max(0.0, max(Lambda[nSpecies+nDim]-(ProjVelocity_i-V_i[A_INDEX]),
  //                                     (ProjVelocity_j-V_j[A_INDEX])-Lambda[nSpecies+nDim]));
  //  Epsilon[nSpecies+nDim+1] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i,
  //                                              ProjVelocity_j-Lambda[iDim] ));
  //  for (iVar = 0; iVar < nVar; iVar++)
  //    if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
  //      Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
  //    else
  //      Lambda[iVar] = fabs(Lambda[iVar]);

  //  for (iVar = 0; iVar < nVar; iVar++)
  //    Lambda[iVar] = fabs(Lambda[iVar]);

  /*--- Calculate inviscid projected Jacobians ---*/
  // Note: Scaling value is 0.5 because inviscid flux is based on 0.5*(Fc_i+Fc_j)
  //  GetInviscidProjJac(U_i, V_i, dPdU_i, Normal, 0.5, val_Jacobian_i);
  //  GetInviscidProjJac(U_j, V_j, dPdU_j, Normal, 0.5, val_Jacobian_j);

  /*--- Roe's Flux approximation ---*/
  //  for (iVar = 0; iVar < nVar; iVar++) {
  //    for (jVar = 0; jVar < nVar; jVar++) {

  /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
  //      Proj_ModJac_Tensor_ij = 0.0;
  //      for (kVar = 0; kVar < nVar; kVar++)
  //        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
  //      val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
  //      val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
  //    }
  //  }
  //}

  /*--- AUSM's Jacobian....requires tiny CFL's (this must be fixed) ---*/
//  if (implicit)  {
//
//    /*--- Initialize the Jacobians ---*/
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_i[iVar][jVar] = 0.0;
//        val_Jacobian_j[iVar][jVar] = 0.0;
//      }
//    }
//
//    if (mF >= 0.0) FcLR = FcL;
//    else           FcLR = FcR;
//
//    /*--- Sound speed derivatives: Species density ---*/
//    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
//      Cvtrs = (3.0/2.0+xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
//      daL[iSpecies] = 1.0/(2.0*aF) * (1/rhoCvtr_i*(Ru/Ms[iSpecies] - Cvtrs*dPdU_i[nSpecies+nDim])*P_i/rho_i
//          + 1.0/rho_i*(1.0+dPdU_i[nSpecies+nDim])*(dPdU_i[iSpecies] - P_i/rho_i));
//      daR[iSpecies] = 1.0/(2.0*aF) * (1/rhoCvtr_j*(Ru/Ms[iSpecies] - Cvtrs*dPdU_j[nSpecies+nDim])*P_j/rho_j
//          + 1.0/rho_j*(1.0+dPdU_j[nSpecies+nDim])*(dPdU_j[iSpecies] - P_j/rho_j));
//    }
//    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
//      daL[nSpecies-1] = 1.0/(2.0*aF*rho_i) * (1+dPdU_i[nSpecies+nDim])*(dPdU_i[nSpecies-1] - P_i/rho_i);
//      daR[nSpecies-1] = 1.0/(2.0*aF*rho_j) * (1+dPdU_j[nSpecies+nDim])*(dPdU_j[nSpecies-1] - P_j/rho_j);
//    }
//
//    /*--- Sound speed derivatives: Momentum ---*/
//    for (iDim = 0; iDim < nDim; iDim++) {
//      daL[nSpecies+iDim] = -1.0/(2.0*rho_i*aF) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim])*u_i[iDim];
//      daR[nSpecies+iDim] = -1.0/(2.0*rho_j*aF) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim])*u_j[iDim];
//    }
//
//    /*--- Sound speed derivatives: Energy ---*/
//    daL[nSpecies+nDim]   = 1.0/(2.0*rho_i*aF) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim]);
//    daR[nSpecies+nDim]   = 1.0/(2.0*rho_j*aF) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim]);
//
//    /*--- Sound speed derivatives: Vib-el energy ---*/
//    daL[nSpecies+nDim+1] = 1.0/(2.0*rho_i*aF) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim+1]);
//    daR[nSpecies+nDim+1] = 1.0/(2.0*rho_j*aF) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim+1]);
//
//    /*--- Left state Jacobian ---*/
//    if (mF >= 0) {
//
//      /*--- Jacobian contribution: dFc terms ---*/
//      for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
//        for (jVar = 0; jVar < nVar; jVar++) {
//          val_Jacobian_i[iVar][jVar] += mF * FcL[iVar]/aF * daL[jVar];
//        }
//        val_Jacobian_i[iVar][iVar] += mF * aF;
//      }
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        val_Jacobian_i[nSpecies+nDim][iSpecies] += mF * (dPdU_i[iSpecies]*aF + rho_i*h_i*daL[iSpecies]);
//      }
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_Jacobian_i[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_i[nSpecies+nDim]*u_i[iDim]*aF + rho_i*h_i*daL[nSpecies+iDim]);
//      }
//      val_Jacobian_i[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_i[nSpecies+nDim])*aF + rho_i*h_i*daL[nSpecies+nDim]);
//      val_Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_i[nSpecies+nDim+1]*aF + rho_i*h_i*daL[nSpecies+nDim+1]);
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_i[nSpecies+nDim+1][jVar] +=  mF * FcL[nSpecies+nDim+1]/aF * daL[jVar];
//      }
//      val_Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += mF * aF;
//    }
//
//
//    /*--- Calculate derivatives of the split pressure flux ---*/
//    if ( (mF >= 0) || ((mF < 0)&&(fabs(mF) <= 1.0)) ) {
//      if (fabs(mL) <= 1.0) {
//
//        /*--- Mach number ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmLP[iSpecies] = 0.5*(mL+1.0) * (-ProjVel_i/(rho_i*aF) - ProjVel_i*daL[iSpecies]/(aF*aF));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmLP[nSpecies+iDim] = 0.5*(mL+1.0) * (-ProjVel_i/(aF*aF) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*aF));
//        dmLP[nSpecies+nDim]   = 0.5*(mL+1.0) * (-ProjVel_i/(aF*aF) * daL[nSpecies+nDim]);
//        dmLP[nSpecies+nDim+1] = 0.5*(mL+1.0) * (-ProjVel_i/(aF*aF) * daL[nSpecies+nDim+1]);
//
//        /*--- Pressure ---*/
//        for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dpLP[iSpecies] = 0.25*(mL+1.0) * (dPdU_i[iSpecies]*(mL+1.0)*(2.0-mL)
//                                            + P_i*(-ProjVel_i/(rho_i*aF)
//                                                   -ProjVel_i*daL[iSpecies]/(aF*aF))*(3.0-3.0*mL));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dpLP[nSpecies+iDim] = 0.25*(mL+1.0) * (-u_i[iDim]*dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
//              + P_i*( -ProjVel_i/(aF*aF) * daL[nSpecies+iDim]
//              + UnitNormal[iDim]/(rho_i*aF))*(3.0-3.0*mL));
//        dpLP[nSpecies+nDim]   = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
//            + P_i*(-ProjVel_i/(aF*aF) * daL[nSpecies+nDim])*(3.0-3.0*mL));
//        dpLP[nSpecies+nDim+1] = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim+1]*(mL+1.0)*(2.0-mL)
//            + P_i*(-ProjVel_i/(aF*aF) * daL[nSpecies+nDim+1])*(3.0-3.0*mL));
//      } else {
//
//        /*--- Mach number ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmLP[iSpecies]      = -ProjVel_i/(rho_i*aF) - ProjVel_i*daL[iSpecies]/(aF*aF);
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmLP[nSpecies+iDim] = -ProjVel_i/(aF*aF) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*aF);
//        dmLP[nSpecies+nDim]   = -ProjVel_i/(aF*aF) * daL[nSpecies+nDim];
//        dmLP[nSpecies+nDim+1] = -ProjVel_i/(aF*aF) * daL[nSpecies+nDim+1];
//
//        /*--- Pressure ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dpLP[iSpecies] = dPdU_i[iSpecies];
//        for (iDim = 0; iDim < nDim; iDim++)
//          dpLP[nSpecies+iDim] = (-u_i[iDim]*dPdU_i[nSpecies+nDim]);
//        dpLP[nSpecies+nDim]   = dPdU_i[nSpecies+nDim];
//        dpLP[nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1];
//      }
//
//      /*--- dM contribution ---*/
//      for (iVar = 0; iVar < nVar; iVar++) {
//        for (jVar = 0; jVar < nVar; jVar++) {
//          val_Jacobian_i[iVar][jVar] += dmLP[jVar]*FcLR[iVar];
//        }
//      }
//
//      /*--- Jacobian contribution: dP terms ---*/
//      for (iDim = 0; iDim < nDim; iDim++) {
//        for (iVar = 0; iVar < nVar; iVar++) {
//          val_Jacobian_i[nSpecies+iDim][iVar] += dpLP[iVar]*UnitNormal[iDim];
//        }
//      }
//    }
//
//    /*--- Right state Jacobian ---*/
//    if (mF < 0) {
//
//      /*--- Jacobian contribution: dFc terms ---*/
//      for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
//        for (jVar = 0; jVar < nVar; jVar++) {
//          val_Jacobian_j[iVar][jVar] += mF * FcR[iVar]/aF * daR[jVar];
//        }
//        val_Jacobian_j[iVar][iVar] += mF * aF;
//      }
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        val_Jacobian_j[nSpecies+nDim][iSpecies] += mF * (dPdU_j[iSpecies]*aF + rho_j*h_j*daR[iSpecies]);
//      }
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_Jacobian_j[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_j[nSpecies+nDim]*u_j[iDim]*aF + rho_j*h_j*daR[nSpecies+iDim]);
//      }
//      val_Jacobian_j[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_j[nSpecies+nDim])*aF + rho_j*h_j*daR[nSpecies+nDim]);
//      val_Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_j[nSpecies+nDim+1]*aF + rho_j*h_j*daR[nSpecies+nDim+1]);
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_j[nSpecies+nDim+1][jVar] +=  mF * FcR[nSpecies+nDim+1]/aF * daR[jVar];
//      }
//      val_Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += mF * aF;
//    }
//
//    /*--- Calculate derivatives of the split pressure flux ---*/
//    if ( (mF < 0) || ((mF >= 0)&&(fabs(mF) <= 1.0)) ) {
//      if (fabs(mR) <= 1.0) {
//
//        /*--- Mach ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmRM[iSpecies] = -0.5*(mR-1.0) * (-ProjVel_j/(rho_j*aF) - ProjVel_j*daR[iSpecies]/(aF*aF));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmRM[nSpecies+iDim] = -0.5*(mR-1.0) * (-ProjVel_j/(aF*aF) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*aF));
//        dmRM[nSpecies+nDim]   = -0.5*(mR-1.0) * (-ProjVel_j/(aF*aF) * daR[nSpecies+nDim]);
//        dmRM[nSpecies+nDim+1] = -0.5*(mR-1.0) * (-ProjVel_j/(aF*aF) * daR[nSpecies+nDim+1]);
//
//        /*--- Pressure ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dpRM[iSpecies] = 0.25*(mR-1.0) * (dPdU_j[iSpecies]*(mR-1.0)*(2.0+mR)
//                                            + P_j*(-ProjVel_j/(rho_j*aF)
//                                                   -ProjVel_j*daR[iSpecies]/(aF*aF))*(3.0+3.0*mR));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dpRM[nSpecies+iDim] = 0.25*(mR-1.0) * ((-u_j[iDim]*dPdU_j[nSpecies+nDim])*(mR-1.0)*(2.0+mR)
//              + P_j*( -ProjVel_j/(aF*aF) * daR[nSpecies+iDim]
//              + UnitNormal[iDim]/(rho_j*aF))*(3.0+3.0*mR));
//        dpRM[nSpecies+nDim]   = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim]*(mR-1.0)*(2.0+mR)
//            + P_j*(-ProjVel_j/(aF*aF)*daR[nSpecies+nDim])*(3.0+3.0*mR));
//        dpRM[nSpecies+nDim+1] = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim+1]*(mR-1.0)*(2.0+mR)
//            + P_j*(-ProjVel_j/(aF*aF) * daR[nSpecies+nDim+1])*(3.0+3.0*mR));
//
//      } else {
//
//        /*--- Mach ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmRM[iSpecies]      = -ProjVel_j/(rho_j*aF) - ProjVel_j*daR[iSpecies]/(aF*aF);
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmRM[nSpecies+iDim] = -ProjVel_j/(aF*aF) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*aF);
//        dmRM[nSpecies+nDim]   = -ProjVel_j/(aF*aF) * daR[nSpecies+nDim];
//        dmRM[nSpecies+nDim+1] = -ProjVel_j/(aF*aF) * daR[nSpecies+nDim+1];
//
//        /*--- Pressure ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dpRM[iSpecies] = dPdU_j[iSpecies];
//        for (iDim = 0; iDim < nDim; iDim++)
//          dpRM[nSpecies+iDim] = -u_j[iDim]*dPdU_j[nSpecies+nDim];
//        dpRM[nSpecies+nDim]   = dPdU_j[nSpecies+nDim];
//        dpRM[nSpecies+nDim+1] = dPdU_j[nSpecies+nDim+1];
//      }
//
//      /*--- Jacobian contribution: dM terms ---*/
//      for (iVar = 0; iVar < nVar; iVar++) {
//        for (jVar = 0; jVar < nVar; jVar++) {
//          val_Jacobian_j[iVar][jVar] += dmRM[jVar] * FcLR[iVar];
//        }
//      }
//
//      /*--- Jacobian contribution: dP terms ---*/
//      for (iDim = 0; iDim < nDim; iDim++) {
//        for (iVar = 0; iVar < nVar; iVar++) {
//          val_Jacobian_j[nSpecies+iDim][iVar] += dpRM[iVar]*UnitNormal[iDim];
//        }
//      }
//    }
//
//    /*--- Integrate over dual-face area ---*/
//    for (iVar = 0; iVar < nVar; iVar++) {
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_i[iVar][jVar] *= Area;
//        val_Jacobian_j[iVar][jVar] *= Area;
//      }
//    }
//  }

  return ResidualType<>(Flux, nullptr, nullptr);
}
