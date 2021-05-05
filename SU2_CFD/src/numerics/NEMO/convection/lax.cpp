/*!
 * \file lax.cpp
 * \brief Implementations of Lax centered scheme.
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

#include "../../../../include/numerics/NEMO/convection/lax.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CCentLax_NEMO::CCentLax_NEMO(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad,
                             CConfig *config) : CNEMONumerics(val_nDim, val_nVar,
                                                              val_nPrimVar,
                                                              val_nPrimVarGrad,
                                                              config) {

  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();

  /*--- Allocate some structures ---*/
  Diff_U   = new su2double[nVar];
  MeanU    = new su2double[nVar];
  MeanV    = new su2double[nPrimVar];
  ProjFlux = new su2double[nVar];
  Flux     = new su2double[nVar];

}

CCentLax_NEMO::~CCentLax_NEMO(void) {

  delete [] Diff_U;
  delete [] MeanU;
  delete [] MeanV;
  delete [] ProjFlux;
  delete [] Flux;
}

CNumerics::ResidualType<> CCentLax_NEMO::ComputeResidual(const CConfig *config) {

  su2double rho_i, rho_j, h_i, h_j, a_i, a_j;
  su2double ProjVel_i, ProjVel_j;

  /*--- Calculate geometrical quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Rename for convenience ---*/
  rho_i = V_i[RHO_INDEX]; rho_j = V_j[RHO_INDEX];
  h_i   = V_i[H_INDEX];   h_j   = V_j[H_INDEX];
  a_i   = V_i[A_INDEX];   a_j   = V_j[A_INDEX];

  /*--- Compute mean quantities for the variables ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    MeanU[iVar] = 0.5*(U_i[iVar]+U_j[iVar]);
  for (iVar = 0; iVar < nPrimVar; iVar++)
    MeanV[iVar] = 0.5*(V_i[iVar]+V_j[iVar]);

  /*--- Compute NonEq specific variables ---*/
  //vector<su2double> mean_eves = fluidmodel->ComputeSpeciesEve(MeanV[TVE_INDEX]);
  //fluidmodel->ComputedPdU(MeanV, mean_eves, MeandPdU);

  /*--- Get projected flux tensor ---*/
  GetInviscidProjFlux(MeanU, MeanV, Normal, ProjFlux);

  /*--- Compute the local spectral radius and the stretching factor ---*/
  ProjVel_i = 0; ProjVel_j = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVel_i += V_i[VEL_INDEX+iDim]*Normal[iDim];
    ProjVel_j += V_j[VEL_INDEX+iDim]*Normal[iDim];
  }
  Area = GeometryToolbox::Norm(nDim, Normal);

  /*--- Dissipation --*/
  Local_Lambda_i = (fabs(ProjVel_i)+a_i*Area);
  Local_Lambda_j = (fabs(ProjVel_j)+a_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);

  /*--- Computes differences btw. conservative variables ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = U_i[iVar] - U_j[iVar];
  Diff_U[nSpecies+nDim] = rho_i*h_i - rho_j*h_j;

  /*--- Compute dissipation coefficient ---*/
  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;

  /*--- Compute viscous part of the residual ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Flux[iVar] = ProjFlux[iVar]+Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  }

  return ResidualType<>(Flux, nullptr, nullptr);

}
