/*!
 * \file turb_sources.cpp
 * \brief Implementation of numerics classes for integration of
 *        turbulence source-terms.
 * \author F. Palacios, T. Economon
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

#include "../../../include/numerics/turbulent/turb_sources.hpp"

CSourceBase_TurbSA::CSourceBase_TurbSA(unsigned short val_nDim,
                                       unsigned short val_nVar,
                                       const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE),
  rotating_frame(config->GetRotating_Frame())
{
  /*--- Spalart-Allmaras closure constants ---*/

  cv1_3 = pow(7.1, 3.0);
  k2    = pow(0.41, 2.0);
  cb1   = 0.1355;
  cw2   = 0.3;
  ct3   = 1.2;
  ct4   = 0.5;
  cw3_6 = pow(2.0, 6.0);
  sigma = 2./3.;
  cb2   = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1.0+cb2)/sigma;
  cr1 = 0.5;

  /*--- Setup the Jacobian pointer, we need to return su2double** but
   *    we know the Jacobian is 1x1 so we use this "trick" to avoid
   *    having to dynamically allocate. ---*/

  Jacobian_i = &Jacobian_Buffer;

}

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 const CConfig* config) :
                         CSourceBase_TurbSA(val_nDim, val_nVar, config) {

  transition = (config->GetKind_Trans_Model() == BC);
}

CNumerics::ResidualType<> CSourcePieceWise_TurbSA::ComputeResidual(const CConfig* config) {

//  AD::StartPreacc();
//  AD::SetPreaccIn(V_i, nDim+6);
//  AD::SetPreaccIn(Vorticity_i, nDim);
//  AD::SetPreaccIn(StrainMag_i);
//  AD::SetPreaccIn(TurbVar_i[0]);
//  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
//  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  // Set the boolean here depending on whether the point is closest to a rough wall or not.
  roughwall = (roughness_i > 0.0);

  if (incompressible) {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }

  Residual        = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  Jacobian_i[0]   = 0.0;

  /*--- Evaluate Omega ---*/

  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);

  /*--- Rotational correction term ---*/

  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

  if (dist_i > 1e-10) {

    /*--- Production term ---*/

    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;

    /*--- Modified values for roughness ---*/
    /*--- Ref: Aupoix, B. and Spalart, P. R., "Extensions of the Spalart-Allmaras Turbulence Model to Account for Wall Roughness,"
     * International Journal of Heat and Fluid Flow, Vol. 24, 2003, pp. 454-462. ---*/
    /* --- See https://turbmodels.larc.nasa.gov/spalart.html#sarough for detailed explanation. ---*/

    Ji = TurbVar_i[0]/nu  + cr1*(roughness_i/(dist_i+EPS)); //roughness_i = 0 for smooth walls and Ji remains the same, changes only if roughness is specified.
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);

    /*--- Using a modified relation so as to not change the Shat that depends on fv2. ---*/
    fv2 = 1.0 - TurbVar_i[0]/(nu+TurbVar_i[0]*fv1);   // From NASA turb modeling resource and 2003 paper

    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);

    Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;

//    Original SA model
//    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;

    if (transition) {

      /*--- BC model constants (2020 revision). ---*/
      const su2double chi_1 = 0.002;
      const su2double chi_2 = 50.0;

      su2double tu = config->GetTurbulenceIntensity_FreeStream();

      su2double nu_t = (TurbVar_i[0]*fv1); //S-A variable

      su2double re_v = ((Density_i*pow(dist_i,2.))/(Laminar_Viscosity_i))*Omega;
      su2double re_theta = re_v/2.193;
      su2double re_theta_t = (803.73 * pow((tu + 0.6067),-1.027)); //MENTER correlation
      //re_theta_t = 163.0 + exp(6.91-tu); //ABU-GHANNAM & SHAW correlation

      su2double term1 = sqrt(max(re_theta-re_theta_t,0.)/(chi_1*re_theta_t));
      su2double term2 = sqrt(max((nu_t*chi_2)/nu,0.));
      su2double term_exponential = (term1 + term2);

      Gamma_BC = 1.0 - exp(-term_exponential);

      Production = Gamma_BC*cb1*Shat*TurbVar_i[0]*Volume;
    }
    else {
      Production = cb1*Shat*TurbVar_i[0]*Volume;
    }

    /*--- Destruction term ---*/

    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    g = r + cw2*(pow(r,6.0)-r);
    g_6 =  pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;

//    Original SA model
//    Destruction = (cw1*fw-cb1*ft2/k2)*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/

    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];

    CrossProduction = cb2_sigma*norm2_Grad*Volume;

    Residual = Production - Destruction + CrossProduction;

    /*--- Implicit part, production term ---*/

    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;

    if (transition) {
      Jacobian_i[0] += Gamma_BC*cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    }
    else {
      Jacobian_i[0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    }

    /*--- Implicit part, destruction term ---*/

    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    Jacobian_i[0] -= cw1*(dfw*TurbVar_i[0] +  2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;

  }

//  AD::SetPreaccOut(Residual);
//  AD::EndPreacc();

  return ResidualType<>(&Residual, &Jacobian_i, nullptr);

}

CSourcePieceWise_TurbSA_COMP::CSourcePieceWise_TurbSA_COMP(unsigned short val_nDim,
                                                           unsigned short val_nVar,
                                                           const CConfig* config) :
                              CSourceBase_TurbSA(val_nDim, val_nVar, config), c5(3.5) { }

CNumerics::ResidualType<> CSourcePieceWise_TurbSA_COMP::ComputeResidual(const CConfig* config) {

  //  AD::StartPreacc();
  //  AD::SetPreaccIn(V_i, nDim+6);
  //  AD::SetPreaccIn(Vorticity_i, nDim);
  //  AD::SetPreaccIn(StrainMag_i);
  //  AD::SetPreaccIn(TurbVar_i[0]);
  //  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
  //  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }

  Residual        = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  Jacobian_i[0]   = 0.0;

  /*--- Evaluate Omega ---*/

  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);

  /*--- Rotational correction term ---*/

  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

  if (dist_i > 1e-10) {

    /*--- Production term ---*/

    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);

    Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;

    /*--- Production term ---*/;

    Production = cb1*Shat*TurbVar_i[0]*Volume;

    /*--- Destruction term ---*/

    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    g = r + cw2*(pow(r,6.0)-r);
    g_6 = pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;

    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/

    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];

    CrossProduction = cb2_sigma*norm2_Grad*Volume;

    Residual = Production - Destruction + CrossProduction;

    /*--- Compressibility Correction term ---*/
    Pressure_i = V_i[nDim+1];
    SoundSpeed_i = sqrt(Pressure_i*Gamma/Density_i);
    aux_cc=0;
    for(iDim=0;iDim<nDim;++iDim){
      for(jDim=0;jDim<nDim;++jDim){
        aux_cc+=PrimVar_Grad_i[1+iDim][jDim]*PrimVar_Grad_i[1+iDim][jDim];}}
    CompCorrection=c5*(TurbVar_i[0]*TurbVar_i[0]/(SoundSpeed_i*SoundSpeed_i))*aux_cc*Volume;

    Residual -= CompCorrection;

    /*--- Implicit part, production term ---*/

    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    Jacobian_i[0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

    /*--- Implicit part, destruction term ---*/

    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    Jacobian_i[0] -= cw1*(dfw*TurbVar_i[0] + 2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;

    /* Compressibility Correction */
    Jacobian_i[0] -= 2.0*c5*(TurbVar_i[0]/(SoundSpeed_i*SoundSpeed_i))*aux_cc*Volume;

  }

  //  AD::SetPreaccOut(Residual);
  //  AD::EndPreacc();

  return ResidualType<>(&Residual, &Jacobian_i, nullptr);

}

CSourcePieceWise_TurbSA_E::CSourcePieceWise_TurbSA_E(unsigned short val_nDim,
                                                     unsigned short val_nVar,
                                                     const CConfig* config) :
                           CSourceBase_TurbSA(val_nDim, val_nVar, config) { }

CNumerics::ResidualType<> CSourcePieceWise_TurbSA_E::ComputeResidual(const CConfig* config) {

  unsigned short iDim, jDim;

  //  AD::StartPreacc();
  //  AD::SetPreaccIn(V_i, nDim+6);
  //  AD::SetPreaccIn(Vorticity_i, nDim);
  //  AD::SetPreaccIn(StrainMag_i);
  //  AD::SetPreaccIn(TurbVar_i[0]);
  //  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
  //  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }

  Residual        = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  Jacobian_i[0]   = 0.0;

  /*
   From NASA Turbulence model site. http://turbmodels.larc.nasa.gov/spalart.html
   This form was developed primarily to improve the near-wall numerical behavior of the model (i.e., the goal was to improve the convergence behavior). The reference is:
   Edwards, J. R. and Chandra, S. "Comparison of Eddy Viscosity-Transport Turbulence Models for Three-Dimensional, Shock-Separated Flowfields," AIAA Journal, Vol. 34, No. 4, 1996, pp. 756-763.
   In this modificaton Omega is replaced by Strain Rate
   */

  /*--- Evaluate Omega, here Omega is the Strain Rate ---*/

  Sbar = 0.0;
  for(iDim=0;iDim<nDim;++iDim){
    for(jDim=0;jDim<nDim;++jDim){
      Sbar+= (PrimVar_Grad_i[1+iDim][jDim]+PrimVar_Grad_i[1+jDim][iDim])*(PrimVar_Grad_i[1+iDim][jDim]);}}
  for(iDim=0;iDim<nDim;++iDim){
    Sbar-= ((2.0/3.0)*pow(PrimVar_Grad_i[1+iDim][iDim],2.0));}

  Omega= sqrt(max(Sbar,0.0));

  /*--- Rotational correction term ---*/

  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

  if (dist_i > 1e-10) {

    /*--- Production term ---*/

    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);

    //Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    Shat = max(S*((1.0/max(Ji,1.0e-16))+fv1),1.0e-16);

    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;

    /*--- Production term ---*/;

    Production = cb1*Shat*TurbVar_i[0]*Volume;

    /*--- Destruction term ---*/

    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    r=tanh(r)/tanh(1.0);

    g = r + cw2*(pow(r,6.0)-r);
    g_6 = pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;

    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/

    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];

    CrossProduction = cb2_sigma*norm2_Grad*Volume;

    Residual = Production - Destruction + CrossProduction;

    /*--- Implicit part, production term ---*/

    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);

    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = -S*pow(Ji,-2.0)/nu + S*dfv1;
    Jacobian_i[0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

    /*--- Implicit part, destruction term ---*/

    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    dr=(1-pow(tanh(r),2.0))*(dr)/tanh(1.0);
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    Jacobian_i[0] -= cw1*(dfw*TurbVar_i[0] + 2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;

  }

  //  AD::SetPreaccOut(Residual);
  //  AD::EndPreacc();

  return ResidualType<>(&Residual, &Jacobian_i, nullptr);

}

CSourcePieceWise_TurbSA_E_COMP::CSourcePieceWise_TurbSA_E_COMP(unsigned short val_nDim,
                                                               unsigned short val_nVar,
                                                               const CConfig* config) :
                                CSourceBase_TurbSA(val_nDim, val_nVar, config) { }

CNumerics::ResidualType<> CSourcePieceWise_TurbSA_E_COMP::ComputeResidual(const CConfig* config) {

  unsigned short iDim;

  //  AD::StartPreacc();
  //  AD::SetPreaccIn(V_i, nDim+6);
  //  AD::SetPreaccIn(Vorticity_i, nDim);
  //  AD::SetPreaccIn(StrainMag_i);
  //  AD::SetPreaccIn(TurbVar_i[0]);
  //  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
  //  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }

  Residual        = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  Jacobian_i[0]   = 0.0;

  /*
   From NASA Turbulence model site. http://turbmodels.larc.nasa.gov/spalart.html
   This form was developed primarily to improve the near-wall numerical behavior of the model (i.e., the goal was to improve the convergence behavior). The reference is:
   Edwards, J. R. and Chandra, S. "Comparison of Eddy Viscosity-Transport Turbulence Models for Three-Dimensional, Shock-Separated Flowfields," AIAA Journal, Vol. 34, No. 4, 1996, pp. 756-763.
   In this modificaton Omega is replaced by Strain Rate
   */

  /*--- Evaluate Omega, here Omega is the Strain Rate ---*/

  Sbar = 0.0;
  for(iDim=0;iDim<nDim;++iDim){
    for(jDim=0;jDim<nDim;++jDim){
      Sbar+= (PrimVar_Grad_i[1+iDim][jDim]+PrimVar_Grad_i[1+jDim][iDim])*(PrimVar_Grad_i[1+iDim][jDim]);}}
  for(iDim=0;iDim<nDim;++iDim){
    Sbar-= ((2.0/3.0)*pow(PrimVar_Grad_i[1+iDim][iDim],2.0));}

  Omega= sqrt(max(Sbar,0.0));

  /*--- Rotational correction term ---*/

  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

  if (dist_i > 1e-10) {

    /*--- Production term ---*/

    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);

    Shat = max(S*((1.0/max(Ji,1.0e-16))+fv1),1.0e-16);

    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;

    /*--- Production term ---*/;

    Production = cb1*Shat*TurbVar_i[0]*Volume;

    /*--- Destruction term ---*/

    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    r=tanh(r)/tanh(1.0);

    g = r + cw2*(pow(r,6.0)-r);
    g_6 = pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;

    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/

    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];

    CrossProduction = cb2_sigma*norm2_Grad*Volume;

    Residual = Production - Destruction + CrossProduction;

    /*--- Compressibility Correction term ---*/
    Pressure_i = V_i[nDim+1];
    SoundSpeed_i = sqrt(Pressure_i*Gamma/Density_i);
    aux_cc=0;
    for(iDim=0;iDim<nDim;++iDim){
        for(jDim=0;jDim<nDim;++jDim){
            aux_cc+=PrimVar_Grad_i[1+iDim][jDim]*PrimVar_Grad_i[1+iDim][jDim];}}
    CompCorrection=c5*(TurbVar_i[0]*TurbVar_i[0]/(SoundSpeed_i*SoundSpeed_i))*aux_cc*Volume;

    Residual -= CompCorrection;

    /*--- Implicit part, production term ---*/

    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);

    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = -S*pow(Ji,-2.0)/nu + S*dfv1;
    Jacobian_i[0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

    /*--- Implicit part, destruction term ---*/

    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    dr=(1-pow(tanh(r),2.0))*(dr)/tanh(1.0);
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    Jacobian_i[0] -= cw1*(dfw*TurbVar_i[0] + 2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;

    /* Compressibility Correction */
    Jacobian_i[0] -= 2.0*c5*(TurbVar_i[0]/(SoundSpeed_i*SoundSpeed_i))*aux_cc*Volume;

  }

  //  AD::SetPreaccOut(Residual);
  //  AD::EndPreacc();

  return ResidualType<>(&Residual, &Jacobian_i, nullptr);

}

CSourcePieceWise_TurbSA_Neg::CSourcePieceWise_TurbSA_Neg(unsigned short val_nDim,
                                                         unsigned short val_nVar,
                                                         const CConfig* config) :
                             CSourceBase_TurbSA(val_nDim, val_nVar, config) { }

CNumerics::ResidualType<> CSourcePieceWise_TurbSA_Neg::ComputeResidual(const CConfig* config) {

  unsigned short iDim;

//  AD::StartPreacc();
//  AD::SetPreaccIn(V_i, nDim+6);
//  AD::SetPreaccIn(Vorticity_i, nDim);
//  AD::SetPreaccIn(StrainMag_i);
//  AD::SetPreaccIn(TurbVar_i[0]);
//  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
//  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }

  Residual        = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  Jacobian_i[0]   = 0.0;

  /*--- Evaluate Omega ---*/

  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);

  /*--- Rotational correction term ---*/

  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }

  if (dist_i > 1e-10) {

    if (TurbVar_i[0] > 0.0) {

      /*--- Production term ---*/

      dist_i_2 = dist_i*dist_i;
      nu = Laminar_Viscosity_i/Density_i;
      Ji = TurbVar_i[0]/nu;
      Ji_2 = Ji*Ji;
      Ji_3 = Ji_2*Ji;
      fv1 = Ji_3/(Ji_3+cv1_3);
      fv2 = 1.0 - Ji/(1.0+Ji*fv1);
      ft2 = ct3*exp(-ct4*Ji_2);
      S = Omega;
      inv_k2_d2 = 1.0/(k2*dist_i_2);

      Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
      Shat = max(Shat, 1.0e-10);
      inv_Shat = 1.0/Shat;

      /*--- Production term ---*/;

      //    Original SA model
      //    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;

      Production = cb1*Shat*TurbVar_i[0]*Volume;

      /*--- Destruction term ---*/

      r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
      g = r + cw2*(pow(r,6.0)-r);
      g_6 =  pow(g,6.0);
      glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
      fw = g*glim;

      Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

      /*--- Diffusion term ---*/

      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];

      CrossProduction = cb2_sigma*norm2_Grad*Volume;

      Residual = Production - Destruction + CrossProduction;

      /*--- Implicit part, production term ---*/

      dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
      dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
      if ( Shat <= 1.0e-10 ) dShat = 0.0;
      else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
      Jacobian_i[0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

      /*--- Implicit part, destruction term ---*/

      dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
      if (r == 10.0) dr = 0.0;
      dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
      dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
      Jacobian_i[0] -= cw1*(dfw*TurbVar_i[0] +  2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;

    }

    else {

      /*--- Production term ---*/

      dist_i_2 = dist_i*dist_i;

      /*--- Production term ---*/;

      Production = cb1*(1.0-ct3)*Omega*TurbVar_i[0]*Volume;

      /*--- Destruction term ---*/

      Destruction = cw1*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

      /*--- Diffusion term ---*/

      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];

      CrossProduction = cb2_sigma*norm2_Grad*Volume;

      Residual = Production + Destruction + CrossProduction;

      /*--- Implicit part, production term ---*/

      Jacobian_i[0] += cb1*(1.0-ct3)*Omega*Volume;

      /*--- Implicit part, destruction term ---*/

      Jacobian_i[0] += 2.0*cw1*TurbVar_i[0]/dist_i_2*Volume;

    }

  }

//  AD::SetPreaccOut(Residual);
//  AD::EndPreacc();

  return ResidualType<>(&Residual, &Jacobian_i, nullptr);

}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim,
                                                   unsigned short val_nVar,
                                                   const su2double *constants,
                                                   su2double val_kine_Inf,
                                                   su2double val_omega_Inf,
                                                   const CConfig* config) :
                          CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  sustaining_terms = (config->GetKind_Turb_Model() == SST_SUST);
  axisymmetric = config->GetAxisymmetric();

  /*--- Closure constants ---*/
  sigma_k_1     = constants[0];
  sigma_k_2     = constants[1];
  sigma_w_1     = constants[2];
  sigma_w_2     = constants[3];
  beta_1        = constants[4];
  beta_2        = constants[5];
  beta_star     = constants[6];
  a1            = constants[7];
  alfa_1        = constants[8];
  alfa_2        = constants[9];

  /*--- Set the ambient values of k and omega to the free stream values. ---*/
  kAmb     = val_kine_Inf;
  omegaAmb = val_omega_Inf;

  /*--- "Allocate" the Jacobian using the static buffer. ---*/
  Jacobian_i[0] = Jacobian_Buffer;
  Jacobian_i[1] = Jacobian_Buffer+2;

}

CNumerics::ResidualType<> CSourcePieceWise_TurbSST::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F2_i); AD::SetPreaccIn(CDkw_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  AD::SetPreaccIn(Vorticity_i, 3);

  unsigned short iDim;
  su2double alfa_blended, beta_blended;
  su2double diverg, pk, pw, zeta;
  su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                Vorticity_i[1]*Vorticity_i[1] +
                                Vorticity_i[2]*Vorticity_i[2]);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];
    Eddy_Viscosity_i = V_i[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }

  Residual[0] = 0.0;       Residual[1] = 0.0;
  Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;

  /*--- Computation of blended constants for the source terms---*/

  alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
  beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;

  if (dist_i > 1e-10) {

   /*--- Production ---*/

   diverg = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
     diverg += PrimVar_Grad_i[iDim+1][iDim];

   /* if using UQ methodolgy, calculate production using perturbed Reynolds stress matrix */

   if (using_uq){
     ComputePerturbedRSM(nDim, Eig_Val_Comp, uq_permute, uq_delta_b, uq_urlx,
                         PrimVar_Grad_i+1, Density_i, Eddy_Viscosity_i,
                         TurbVar_i[0], MeanPerturbedRSM);
     SetPerturbedStrainMag(TurbVar_i[0]);
     pk = Eddy_Viscosity_i*PerturbedStrainMag*PerturbedStrainMag
          - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
   }
   else {
     pk = Eddy_Viscosity_i*StrainMag_i*StrainMag_i - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
   }

   pk = min(pk,20.0*beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]);
   pk = max(pk,0.0);

   zeta = max(TurbVar_i[1], VorticityMag*F2_i/a1);

   /* if using UQ methodolgy, calculate production using perturbed Reynolds stress matrix */

   if (using_uq){
     pw = PerturbedStrainMag * PerturbedStrainMag - 2.0/3.0*zeta*diverg;
   }
   else {
     pw = StrainMag_i*StrainMag_i - 2.0/3.0*zeta*diverg;
   }
   pw = alfa_blended*Density_i*max(pw,0.0);

   /*--- Sustaining terms, if desired. Note that if the production terms are
         larger equal than the sustaining terms, the original formulation is
         obtained again. This is in contrast to the version in literature
         where the sustaining terms are simply added. This latter approach could
         lead to problems for very big values of the free-stream turbulence
         intensity. ---*/

   if ( sustaining_terms ) {
     const su2double sust_k = beta_star*Density_i*kAmb*omegaAmb;
     const su2double sust_w = beta_blended*Density_i*omegaAmb*omegaAmb;

     pk = max(pk, sust_k);
     pw = max(pw, sust_w);
   }

   /*--- Add the production terms to the residuals. ---*/

   Residual[0] += pk*Volume;
   Residual[1] += pw*Volume;

   /*--- Dissipation ---*/

   Residual[0] -= beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]*Volume;
   Residual[1] -= beta_blended*Density_i*TurbVar_i[1]*TurbVar_i[1]*Volume;

   /*--- Cross diffusion ---*/

   Residual[1] += (1.0 - F1_i)*CDkw_i*Volume;
   
   /*--- Contribution due to 2D axisymmetric formulation ---*/
  
   if (axisymmetric) ResidualAxisymmetric(alfa_blended,zeta);

   /*--- Implicit part ---*/

   Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;
   Jacobian_i[0][1] = -beta_star*TurbVar_i[0]*Volume;
   Jacobian_i[1][0] = 0.0;
   Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
  }

  AD::SetPreaccOut(Residual, nVar);
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);

}

void CSourcePieceWise_TurbSST::SetPerturbedStrainMag(su2double turb_ke){

  /*--- Compute norm of perturbed strain rate tensor. ---*/

  PerturbedStrainMag = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    for (unsigned short jDim = 0; jDim < nDim; jDim++){
      su2double StrainRate_ij = MeanPerturbedRSM[iDim][jDim] - TWO3 * turb_ke * delta[iDim][jDim];
      StrainRate_ij = - StrainRate_ij * Density_i / (2 * Eddy_Viscosity_i);

      PerturbedStrainMag += pow(StrainRate_ij, 2.0);
    }
  }
  PerturbedStrainMag = sqrt(2.0*PerturbedStrainMag);

}
