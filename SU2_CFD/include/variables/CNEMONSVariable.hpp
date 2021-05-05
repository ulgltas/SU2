﻿/*!
 * \file CNEMONSVariable.hpp
 * \brief Class for defining the variables of the compressible NEMO Navier-Stokes solver.
 * \author C. Garbacz, W. Maier, S.R. Copeland.
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

#pragma once

#include "CNEMOEulerVariable.hpp"

/*!
 * \class CNEMONSVariable
 * \brief Main class for defining the variables of the NEMO Navier-Stokes' solver.
 * \ingroup Navier_Stokes_Equations
 * \author C. Garbacz, W. Maier, S.R. Copeland.
 * \version 7.0.8
 */
class CNEMONSVariable final : public CNEMOEulerVariable {
private:
  VectorType Prandtl_Lam;       /*!< \brief Laminar Prandtl number. */
  VectorType Temperature_Ref;   /*!< \brief Reference temperature of the fluid. */
  VectorType Viscosity_Ref;     /*!< \brief Reference viscosity of the fluid. */
  VectorType Viscosity_Inf;     /*!< \brief Viscosity of the fluid at the infinity. */
  MatrixType DiffusionCoeff;    /*!< \brief Diffusion coefficient of the mixture. */
  CVectorOfMatrix Dij;            /*!< \brief Binary diffusion coefficients. */
  VectorType LaminarViscosity;  /*!< \brief Viscosity of the fluid. */
  VectorType ThermalCond;       /*!< \brief T-R thermal conductivity of the gas mixture. */
  VectorType ThermalCond_ve;    /*!< \brief V-E thermal conductivity of the gas mixture. */
  vector<su2double> thermalconductivities;
  vector<su2double> Ds;

  su2double inv_TimeScale;      /*!< \brief Inverse of the reference time scale. */

  MatrixType Vorticity;         /*!< \brief Vorticity of the fluid. */
  VectorType StrainMag;         /*!< \brief Magnitude of rate of strain tensor. */
  VectorType Tau_Wall;          /*!< \brief Magnitude of the wall shear stress from a wall function. */
  VectorType DES_LengthScale;   /*!< \brief DES Length Scale. */
  VectorType Roe_Dissipation;   /*!< \brief Roe low dissipation coefficient. */
  VectorType Vortex_Tilting;    /*!< \brief Value of the vortex tilting variable for DES length scale computation. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_density - Value of the flow density (initialization value).
   * \param[in] val_massfrac - Value of the flow mass fraction (initialization value).
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_temperature - Value of the flow temperature (initialization value).
   * \param[in] val_temperature_ve - Value of the flow temperature_ve (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nPrimVar - Number of primitive variables.
   * \param[in] val_nPrimVargrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMONSVariable(su2double val_density, const su2double *val_massfrac, su2double *val_velocity,
                  su2double val_temperature, su2double val_temperature_ve, unsigned long npoint,
                  unsigned long val_nDim, unsigned long val_nVar, unsigned long val_nPrimVar,
                  unsigned long val_nPrimVarGrad, CConfig *config, CNEMOGas *fluidmodel);

  /*!
   * \brief Constructor of the class.
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nPrimVar - Number of primitive variables.
   * \param[in] val_nPrimgVarGrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMONSVariable(su2double *val_solution, unsigned long val_nDim, unsigned long val_nVar,
                  unsigned long val_nPrimVar, unsigned long val_nPrimVarGrad, unsigned long npoint,
                  CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CNEMONSVariable() = default;

    /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive_Aux(void) const { return Primitive_Aux; }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) final;

  /*!
   * \brief Set the vorticity value.
   */
  bool SetVorticity(void);

  /*!
   * \overload
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline void SetEddyViscosity(unsigned long iPoint, su2double eddy_visc) override { Primitive(iPoint,EDDY_VISC_INDEX) = eddy_visc; }

  /*!
   * \brief Get the species diffusion coefficient.
   * \return Value of the species diffusion coefficient.
   */
  inline su2double* GetDiffusionCoeff(unsigned long iPoint) override { return DiffusionCoeff[iPoint]; }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(unsigned long iPoint) const override { return LaminarViscosity(iPoint); }

  /*!
   * \brief Get the eddy viscosity of the flow.
   * \return The eddy viscosity of the flow.
   */
  inline su2double GetEddyViscosity(unsigned long iPoint) const override { return Primitive(iPoint,EDDY_VISC_INDEX); }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(unsigned long iPoint) const override {return ThermalCond(iPoint); }

  /*!
   * \brief Get the vib-el. thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity_ve(unsigned long iPoint) const override { return ThermalCond_ve(iPoint); }

  /*!
   * \brief Set the temperature at the wall
   */
  inline void SetWallTemperature(unsigned long iPoint, su2double temperature_wall) override {
    Primitive(iPoint,T_INDEX) = temperature_wall;
  }

  /*!
   * \brief Get the value of the vorticity.
   * \return Value of the vorticity.
   */
  inline su2double *GetVorticity(unsigned long iPoint) override { return Vorticity[iPoint]; }


};
