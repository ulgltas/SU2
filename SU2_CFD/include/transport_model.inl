/*!
 * \file transport_model.inl
 * \brief In-Line subroutines of the <i>transport_model.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 7.0.5 "Blackbird"
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

#pragma once

inline su2double CViscosityModel::GetViscosity() const { return Mu; }
inline su2double CViscosityModel::Getdmudrho_T () const { return dmudrho_T; }
inline su2double CViscosityModel::GetdmudT_rho() const { return dmudT_rho; }
inline void CViscosityModel::SetViscosity(su2double T, su2double rho) {}
inline void CViscosityModel::SetDerViscosity(su2double T, su2double rho) {}

inline su2double CConductivityModel::GetConductivity() const { return Kt; }
inline su2double CConductivityModel::Getdktdrho_T () const { return dktdrho_T; }
inline su2double CConductivityModel::GetdktdT_rho () const { return dktdT_rho; }
inline void CConductivityModel::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {}
inline void CConductivityModel::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {}
