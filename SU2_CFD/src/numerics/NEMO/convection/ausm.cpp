/*!
 * \file ausm.cpp
 * \brief Implementations of the AUSM-family of schemes in NEMO.
 * \author F. Palacios, S.R. Copeland, W. Maier, C. Garbacz
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

#include "../../../../include/numerics/NEMO/convection/ausm.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSM_NEMO::CUpwAUSM_NEMO(unsigned short val_nDim, unsigned short val_nVar, 
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad, 
                             CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  FcL    = new su2double [nVar];
  FcR    = new su2double [nVar];
  //dmLP   = new su2double [nVar];
  //dmRM   = new su2double [nVar];
  //dpLP   = new su2double [nVar];
  //dpRM   = new su2double [nVar];
  rhos_i = new su2double [nSpecies];
  rhos_j = new su2double [nSpecies];
  u_i    = new su2double [nDim];
  u_j    = new su2double [nDim];

  Flux   = new su2double[nVar];

}

CUpwAUSM_NEMO::~CUpwAUSM_NEMO(void) {
  
  delete [] FcL;
  delete [] FcR;
  //delete [] dmLP;
  //delete [] dmRM;
  //delete [] dpLP;
  //delete [] dpRM;
  delete [] rhos_i;
  delete [] rhos_j;
  delete [] u_i;
  delete [] u_j;
  delete [] Flux;
}

CNumerics::ResidualType<> CUpwAUSM_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iDim, iVar, iSpecies;
  su2double rho_i, rho_j, 
  e_ve_i, e_ve_j, mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;

  /*--- Compute geometric quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    u_i[iDim] = V_i[VEL_INDEX+iDim];
    u_j[iDim] = V_j[VEL_INDEX+iDim];
  }

  P_i   = V_i[P_INDEX];   P_j   = V_j[P_INDEX];
  h_i   = V_i[H_INDEX];   h_j   = V_j[H_INDEX];
  a_i   = V_i[A_INDEX];   a_j   = V_j[A_INDEX];
  rho_i = V_i[RHO_INDEX]; rho_j = V_j[RHO_INDEX];
  
  e_ve_i  = 0; e_ve_j  = 0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    e_ve_i += (V_i[RHOS_INDEX+iSpecies]*eve_i[iSpecies])/rho_i;
    e_ve_j += (V_j[RHOS_INDEX+iSpecies]*eve_j[iSpecies])/rho_j;
  }

  /*--- Projected velocities ---*/
  ProjVel_i = 0.0; ProjVel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += u_i[iDim]*UnitNormal[iDim];
    ProjVel_j += u_j[iDim]*UnitNormal[iDim];
  }

  /*--- Calculate L/R Mach numbers ---*/
  mL = ProjVel_i/a_i;
  mR = ProjVel_j/a_j;

  /*--- Calculate split numerical fluxes ---*/
  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
  else                 mLP = 0.5*(mL+fabs(mL));

  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
  else                 mRM = 0.5*(mR-fabs(mR));

  mF = mLP + mRM;

  if (fabs(mL) <= 1.0) pLP = 0.25*P_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
  else                 pLP = 0.5*P_i*(mL+fabs(mL))/mL;

  if (fabs(mR) <= 1.0) pRM = 0.25*P_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
  else                 pRM = 0.5*P_j*(mR-fabs(mR))/mR;

  pF = pLP + pRM;
  Phi = fabs(mF);

  /*--- Assign left & right convective vectors ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    FcL[iSpecies] = rhos_i[iSpecies]*a_i;
    FcR[iSpecies] = rhos_j[iSpecies]*a_j;
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    FcL[nSpecies+iDim] = rho_i*a_i*u_i[iDim];
    FcR[nSpecies+iDim] = rho_j*a_j*u_j[iDim];
  }
  FcL[nSpecies+nDim]   = rho_i*a_i*h_i;
  FcR[nSpecies+nDim]   = rho_j*a_j*h_j;
  FcL[nSpecies+nDim+1] = rho_i*a_i*e_ve_i;
  FcR[nSpecies+nDim+1] = rho_j*a_j*e_ve_j;

  /*--- Compute numerical flux ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = 0.5*((mF+Phi)*FcL[iVar]+(mF-Phi)*FcR[iVar])*Area;

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[nSpecies+iDim] += pF*UnitNormal[iDim]*Area;

//  if (implicit)

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
//      daL[iSpecies] = 1.0/(2.0*a_i) * (1/rhoCvtr_i*(Ru/Ms[iSpecies] - Cvtrs*dPdU_i[nSpecies+nDim])*P_i/rho_i
//          + 1.0/rho_i*(1.0+dPdU_i[nSpecies+nDim])*(dPdU_i[iSpecies] - P_i/rho_i));
//      daR[iSpecies] = 1.0/(2.0*a_j) * (1/rhoCvtr_j*(Ru/Ms[iSpecies] - Cvtrs*dPdU_j[nSpecies+nDim])*P_j/rho_j
//          + 1.0/rho_j*(1.0+dPdU_j[nSpecies+nDim])*(dPdU_j[iSpecies] - P_j/rho_j));
//    }
//    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
//      daL[nSpecies-1] = 1.0/(2.0*a_i*rho_i) * (1+dPdU_i[nSpecies+nDim])*(dPdU_i[nSpecies-1] - P_i/rho_i);
//      daR[nSpecies-1] = 1.0/(2.0*a_j*rho_j) * (1+dPdU_j[nSpecies+nDim])*(dPdU_j[nSpecies-1] - P_j/rho_j);
//    }
//
//    /*--- Sound speed derivatives: Momentum ---*/
//    for (iDim = 0; iDim < nDim; iDim++) {
//      daL[nSpecies+iDim] = -1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim])*u_i[iDim];
//      daR[nSpecies+iDim] = -1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim])*u_j[iDim];
//    }
//
//    /*--- Sound speed derivatives: Energy ---*/
//    daL[nSpecies+nDim]   = 1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim]);
//    daR[nSpecies+nDim]   = 1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim]);
//
//    /*--- Sound speed derivatives: Vib-el energy ---*/
//    daL[nSpecies+nDim+1] = 1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim+1]);
//    daR[nSpecies+nDim+1] = 1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim+1]);
//
//    /*--- Left state Jacobian ---*/
//    if (mF >= 0) {
//
//      /*--- Jacobian contribution: dFc terms ---*/
//      for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
//        for (jVar = 0; jVar < nVar; jVar++) {
//          val_Jacobian_i[iVar][jVar] += mF * FcL[iVar]/a_i * daL[jVar];
//        }
//        val_Jacobian_i[iVar][iVar] += mF * a_i;
//      }
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        val_Jacobian_i[nSpecies+nDim][iSpecies] += mF * (dPdU_i[iSpecies]*a_i + rho_i*h_i*daL[iSpecies]);
//      }
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_Jacobian_i[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_i[nSpecies+nDim]*u_i[iDim]*a_i + rho_i*h_i*daL[nSpecies+iDim]);
//      }
//      val_Jacobian_i[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_i[nSpecies+nDim])*a_i + rho_i*h_i*daL[nSpecies+nDim]);
//      val_Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_i[nSpecies+nDim+1]*a_i + rho_i*h_i*daL[nSpecies+nDim+1]);
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_i[nSpecies+nDim+1][jVar] +=  mF * FcL[nSpecies+nDim+1]/a_i * daL[jVar];
//      }
//      val_Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += mF * a_i;
//    }
//
//    /*--- Calculate derivatives of the split pressure flux ---*/
//    if ( (mF >= 0) || ((mF < 0)&&(fabs(mF) <= 1.0)) ) {
//      if (fabs(mL) <= 1.0) {
//
//        /*--- Mach number ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmLP[iSpecies] = 0.5*(mL+1.0) * (-ProjVel_i/(rho_i*a_i) - ProjVel_i*daL[iSpecies]/(a_i*a_i));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmLP[nSpecies+iDim] = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*a_i));
//        dmLP[nSpecies+nDim]   = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim]);
//        dmLP[nSpecies+nDim+1] = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1]);
//
//        /*--- Pressure ---*/
//        for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dpLP[iSpecies] = 0.25*(mL+1.0) * (dPdU_i[iSpecies]*(mL+1.0)*(2.0-mL)
//                                            + P_i*(-ProjVel_i/(rho_i*a_i)
//                                                   -ProjVel_i*daL[iSpecies]/(a_i*a_i))*(3.0-3.0*mL));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dpLP[nSpecies+iDim] = 0.25*(mL+1.0) * (-u_i[iDim]*dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
//              + P_i*( -ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim]
//              + UnitNormal[iDim]/(rho_i*a_i))*(3.0-3.0*mL));
//        dpLP[nSpecies+nDim]   = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
//            + P_i*(-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim])*(3.0-3.0*mL));
//        dpLP[nSpecies+nDim+1] = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim+1]*(mL+1.0)*(2.0-mL)
//            + P_i*(-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1])*(3.0-3.0*mL));
//      } else {
//
//        /*--- Mach number ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmLP[iSpecies]      = -ProjVel_i/(rho_i*a_i) - ProjVel_i*daL[iSpecies]/(a_i*a_i);
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmLP[nSpecies+iDim] = -ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*a_i);
//        dmLP[nSpecies+nDim]   = -ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim];
//        dmLP[nSpecies+nDim+1] = -ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1];
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
//          val_Jacobian_j[iVar][jVar] += mF * FcR[iVar]/a_j * daR[jVar];
//        }
//        val_Jacobian_j[iVar][iVar] += mF * a_j;
//      }
//      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//        val_Jacobian_j[nSpecies+nDim][iSpecies] += mF * (dPdU_j[iSpecies]*a_j + rho_j*h_j*daR[iSpecies]);
//      }
//      for (iDim = 0; iDim < nDim; iDim++) {
//        val_Jacobian_j[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_j[nSpecies+nDim]*u_j[iDim]*a_j + rho_j*h_j*daR[nSpecies+iDim]);
//      }
//      val_Jacobian_j[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_j[nSpecies+nDim])*a_j + rho_j*h_j*daR[nSpecies+nDim]);
//      val_Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_j[nSpecies+nDim+1]*a_j + rho_j*h_j*daR[nSpecies+nDim+1]);
//      for (jVar = 0; jVar < nVar; jVar++) {
//        val_Jacobian_j[nSpecies+nDim+1][jVar] +=  mF * FcR[nSpecies+nDim+1]/a_j * daR[jVar];
//      }
//      val_Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += mF * a_j;
//    }
//
//    /*--- Calculate derivatives of the split pressure flux ---*/
//    if ( (mF < 0) || ((mF >= 0)&&(fabs(mF) <= 1.0)) ) {
//      if (fabs(mR) <= 1.0) {
//
//        /*--- Mach ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmRM[iSpecies] = -0.5*(mR-1.0) * (-ProjVel_j/(rho_j*a_j) - ProjVel_j*daR[iSpecies]/(a_j*a_j));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmRM[nSpecies+iDim] = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*a_j));
//        dmRM[nSpecies+nDim]   = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim]);
//        dmRM[nSpecies+nDim+1] = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1]);
//
//        /*--- Pressure ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dpRM[iSpecies] = 0.25*(mR-1.0) * (dPdU_j[iSpecies]*(mR-1.0)*(2.0+mR)
//                                            + P_j*(-ProjVel_j/(rho_j*a_j)
//                                                   -ProjVel_j*daR[iSpecies]/(a_j*a_j))*(3.0+3.0*mR));
//        for (iDim = 0; iDim < nDim; iDim++)
//          dpRM[nSpecies+iDim] = 0.25*(mR-1.0) * ((-u_j[iDim]*dPdU_j[nSpecies+nDim])*(mR-1.0)*(2.0+mR)
//              + P_j*( -ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim]
//              + UnitNormal[iDim]/(rho_j*a_j))*(3.0+3.0*mR));
//        dpRM[nSpecies+nDim]   = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim]*(mR-1.0)*(2.0+mR)
//            + P_j*(-ProjVel_j/(a_j*a_j)*daR[nSpecies+nDim])*(3.0+3.0*mR));
//        dpRM[nSpecies+nDim+1] = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim+1]*(mR-1.0)*(2.0+mR)
//            + P_j*(-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1])*(3.0+3.0*mR));
//
//      } else {
//
//        /*--- Mach ---*/
//        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//          dmRM[iSpecies]      = -ProjVel_j/(rho_j*a_j) - ProjVel_j*daR[iSpecies]/(a_j*a_j);
//        for (iDim = 0; iDim < nDim; iDim++)
//          dmRM[nSpecies+iDim] = -ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*a_j);
//        dmRM[nSpecies+nDim]   = -ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim];
//        dmRM[nSpecies+nDim+1] = -ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1];
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
