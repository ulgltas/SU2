/*!
 * \file CSU2TCLib.cpp
 * \brief Source of user defined 2T nonequilibrium gas model.
 * \author C. Garbacz, W. Maier, S. R. Copeland
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

#include "../../include/fluid/CSU2TCLib.hpp"
#include "../../../Common/include/option_structure.hpp"

CSU2TCLib::CSU2TCLib(const CConfig* config, unsigned short val_nDim, bool viscous): CNEMOGas(config, val_nDim){

  unsigned short maxEl = 0;
  su2double mf = 0.0;

  const auto MassFrac_Freestream = config->GetGas_Composition();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mf += MassFrac_Freestream[iSpecies];

  /*--- Allocate vectors for gas properties ---*/
  nElStates.resize(nSpecies,0);
  Wall_Catalycity.resize(nSpecies,0.0);
  CharVibTemp.resize(nSpecies,0.0);
  RotationModes.resize(nSpecies,0.0);
  Diss.resize(nSpecies,0.0);
  A.resize(5,0.0);
  Omega00.resize(nSpecies,nSpecies,4,0.0);
  Omega11.resize(nSpecies,nSpecies,4,0.0);
  RxnConstantTable.resize(6,5) = su2double(0.0);
  Blottner.resize(nSpecies,3)  = su2double(0.0);

  if(viscous){
    MolarFracWBE.resize(nSpecies,0.0);
    phis.resize(nSpecies,0.0);
    mus.resize(nSpecies,0.0);
  }

  if (gas_model =="ARGON"){
    if (nSpecies != 1) {
      cout << "CONFIG ERROR: nSpecies mismatch between gas model & gas composition" << endl;
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      cout << "CONFIG ERROR: Intial gas mass fractions do not sum to 1!" << " mf is equal to "<< mf <<endl;
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.667;
    nReactions  = 0;

    // Molar mass [kg/kmol]
    MolarMass[0] = 39.948;
    // Rotational modes of energy storage
    RotationModes[0] = 0.0;
    // Characteristic vibrational temperatures
    CharVibTemp[0] = 0.0;

    Enthalpy_Formation[0] = 0.0;
    Ref_Temperature[0] = 0.0;
    nElStates[0] = 7;

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    /*--- AR: Blottner coefficients. ---*/
    Blottner(0,0) = 3.83444322E-03;   Blottner(0,1) = 6.74718764E-01;   Blottner(0,2) = -1.24290388E+01;

    /*--- AR: 7 states ---*/
    CharElTemp(0,0) = 0.000000000000000E+00;
    CharElTemp(0,1) = 1.611135736988230E+05;
    CharElTemp(0,2) = 1.625833076870950E+05;
    CharElTemp(0,3) = 1.636126382960720E+05;
    CharElTemp(0,4) = 1.642329518358000E+05;
    CharElTemp(0,5) = 1.649426852542080E+05;
    CharElTemp(0,6) = 1.653517702884570E+05;
    ElDegeneracy(0,0) = 1;
    ElDegeneracy(0,1) = 9;
    ElDegeneracy(0,2) = 21;
    ElDegeneracy(0,3) = 7;
    ElDegeneracy(0,4) = 3;
    ElDegeneracy(0,5) = 5;
    ElDegeneracy(0,6) = 15;

  } else if (gas_model == "N2"){
    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 2) {
      cout << "CONFIG ERROR: nSpecies mismatch between gas model & gas composition" << endl;
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      cout << "CONFIG ERROR: Intial gas mass fractions do not sum to 1!" << endl;
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 2;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    /*--- Assign gas properties ---*/
    // Wall mass fractions for catalytic boundaries
    Wall_Catalycity[0] = 0.999;
    Wall_Catalycity[1] = 0.001;
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 0.0;
    // Molar mass [kg/kmol]
    MolarMass[0] = 2.0*14.0067;
    MolarMass[1] = 14.0067;
    // Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 0.0;
    // Formation enthalpy: (JANAF values [KJ/Kmol])
    // J/kg - from Scalabrin
    Enthalpy_Formation[0] = 0.0;      //N2
    Enthalpy_Formation[1] = 3.36E7;   //N
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    // Blottner viscosity coefficients
    // A                       // B                       // C
    Blottner(0,0) = 2.68E-2;   Blottner(0,1) = 3.18E-1;   Blottner(0,2) = -1.13E1;  // N2
    Blottner(1,0) = 1.16E-2;   Blottner(1,1) = 6.03E-1;   Blottner(1,2) = -1.24E1;  // N
    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 3;                     // N
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    /*--- Assign values to data structures ---*/
    // N2: 15 states
    CharElTemp(0,0)  = 0.000000000000000E+00;
    CharElTemp(0,1)  = 7.223156514095200E+04;
    CharElTemp(0,2)  = 8.577862640384000E+04;
    CharElTemp(0,3)  = 8.605026716160000E+04;
    CharElTemp(0,4)  = 9.535118627874400E+04;
    CharElTemp(0,5)  = 9.805635702203200E+04;
    CharElTemp(0,6)  = 9.968267656935200E+04;
    CharElTemp(0,7)  = 1.048976467715200E+05;
    CharElTemp(0,8)  = 1.116489555200000E+05;
    CharElTemp(0,9)  = 1.225836470400000E+05;
    CharElTemp(0,10) = 1.248856873600000E+05;
    CharElTemp(0,11) = 1.282476158188320E+05;
    CharElTemp(0,12) = 1.338060936000000E+05;
    CharElTemp(0,13) = 1.404296391107200E+05;
    CharElTemp(0,14) = 1.504958859200000E+05;
    ElDegeneracy(0,0)  = 1;
    ElDegeneracy(0,1)  = 3;
    ElDegeneracy(0,2)  = 6;
    ElDegeneracy(0,3)  = 6;
    ElDegeneracy(0,4)  = 3;
    ElDegeneracy(0,5)  = 1;
    ElDegeneracy(0,6)  = 2;
    ElDegeneracy(0,7)  = 2;
    ElDegeneracy(0,8)  = 5;
    ElDegeneracy(0,9)  = 1;
    ElDegeneracy(0,10) = 6;
    ElDegeneracy(0,11) = 6;
    ElDegeneracy(0,12) = 10;
    ElDegeneracy(0,13) = 6;
    ElDegeneracy(0,14) = 6;
    // N: 3 states
    CharElTemp(1,0) = 0.000000000000000E+00;
    CharElTemp(1,1) = 2.766469645581980E+04;
    CharElTemp(1,2) = 4.149309313560210E+04;
    ElDegeneracy(1,0) = 4;
    ElDegeneracy(1,1) = 10;
    ElDegeneracy(1,2) = 6;
    /*--- Set Arrhenius coefficients for chemical reactions ---*/
    // Note: Data lists coefficients in (cm^3/mol-s) units, need to convert
    //       to (m^3/kmol-s) to be consistent with the rest of the code
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 3.0E22;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    // Characteristic temperature
    ArrheniusTheta[0] = 113200.0;
    ArrheniusTheta[1] = 113200.0;
    /*--- Set reaction maps ---*/
    // N2 + N2 -> 2N + N2
    Reactions(0,0,0)=0;   Reactions(0,0,1)=0;   Reactions(0,0,2)=nSpecies;
    Reactions(0,1,0)=1;   Reactions(0,1,1)=1;   Reactions(0,1,2) =0;
    // N2 + N -> 2N + N
    Reactions(1,0,0)=0;   Reactions(1,0,1)=1;   Reactions(1,0,2)=nSpecies;
    Reactions(1,1,0)=1;   Reactions(1,1,1)=1;   Reactions(1,1,2)=1;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0] = 0.5; Tcf_b[0] = 0.5; Tcb_a[0] = 1;  Tcb_b[0] = 0;
    Tcf_a[1] = 0.5; Tcf_b[1] = 0.5; Tcb_a[1] = 1;  Tcb_b[1] = 0;
    /*--- Dissociation potential [KJ/kg] ---*/
    Diss[0] = 3.36E4;
    Diss[1] = 0.0;
    /*--- Collision integral data ---*/
    Omega00(0,0,0) = -6.0614558E-03;  Omega00(0,0,1) = 1.2689102E-01;   Omega00(0,0,2) = -1.0616948E+00;  Omega00(0,0,3) = 8.0955466E+02;
    Omega00(0,1,0) = -1.0796249E-02;  Omega00(0,1,1) = 2.2656509E-01;   Omega00(0,1,2) = -1.7910602E+00;  Omega00(0,1,3) = 4.0455218E+03;
    Omega00(1,0,0) = -1.0796249E-02;  Omega00(1,0,1) = 2.2656509E-01;   Omega00(1,0,2) = -1.7910602E+00;  Omega00(1,0,3) = 4.0455218E+03;
    Omega00(1,1,0) = -9.6083779E-03;  Omega00(1,1,1) = 2.0938971E-01;   Omega00(1,1,2) = -1.7386904E+00;  Omega00(1,1,3) = 3.3587983E+03;
    Omega11(0,0,0) = -7.6303990E-03;  Omega11(0,0,1) = 1.6878089E-01;   Omega11(0,0,2) = -1.4004234E+00;  Omega11(0,0,3) = 2.1427708E+03;
    Omega11(0,1,0) = -8.3493693E-03;  Omega11(0,1,1) = 1.7808911E-01;   Omega11(0,1,2) = -1.4466155E+00;  Omega11(0,1,3) = 1.9324210E+03;
    Omega11(1,0,0) = -8.3493693E-03;  Omega11(1,0,1) = 1.7808911E-01;   Omega11(1,0,2) = -1.4466155E+00;  Omega11(1,0,3) = 1.9324210E+03;
    Omega11(1,1,0) = -7.7439615E-03;  Omega11(1,1,1) = 1.7129007E-01;   Omega11(1,1,2) = -1.4809088E+00;  Omega11(1,1,3) = 2.1284951E+03;

  } else if (gas_model == "AIR-5"){

    /*--- Check for errors in the initialization ---*/
    if (nSpecies != 5) {
      cout << "CONFIG ERROR: nSpecies mismatch between gas model & gas composition" << endl;
    }
    mf = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      mf += MassFrac_Freestream[iSpecies];
    if (mf != 1.0) {
      cout << "CONFIG ERROR: Initial gas mass fractions do not sum to 1!" << endl;
    }

    /*--- Define parameters of the gas model ---*/
    gamma       = 1.4;
    nReactions  = 17;

    Reactions.resize(nReactions,2,6,0.0);
    ArrheniusCoefficient.resize(nReactions,0.0);
    ArrheniusEta.resize(nReactions,0.0);
    ArrheniusTheta.resize(nReactions,0.0);
    Tcf_a.resize(nReactions,0.0);
    Tcf_b.resize(nReactions,0.0);
    Tcb_a.resize(nReactions,0.0);
    Tcb_b.resize(nReactions,0.0);

    // Wall mass fractions for catalytic boundaries
    Wall_Catalycity[0] = 0.4;
    Wall_Catalycity[1] = 0.4;
    Wall_Catalycity[2] = 0.1;
    Wall_Catalycity[3] = 0.05;
    Wall_Catalycity[4] = 0.05;
    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 2.0;
    RotationModes[2] = 2.0;
    RotationModes[3] = 0.0;
    RotationModes[4] = 0.0;
    // Molar mass [kg/kmol]
    MolarMass[0] = 2.0*14.0067;
    MolarMass[1] = 2.0*15.9994;
    MolarMass[2] = 14.0067+15.9994;
    MolarMass[3] = 14.0067;
    MolarMass[4] = 15.9994;
    //Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 2239.0;
    CharVibTemp[2] = 2817.0;
    CharVibTemp[3] = 0.0;
    CharVibTemp[4] = 0.0;
    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;      //N2
    Enthalpy_Formation[1] = 0.0;      //O2
    Enthalpy_Formation[2] = 3.0E6;    //NO
    Enthalpy_Formation[3] = 3.36E7;   //N
    Enthalpy_Formation[4] = 1.54E7;   //O
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    //        Ref_Temperature[2] = 298.15;
    //        Ref_Temperature[3] = 298.15;
    //        Ref_Temperature[4] = 298.15
    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner(0,0) = 2.68E-2;   Blottner(0,1) =  3.18E-1;  Blottner(0,2) = -1.13E1;  // N2
    Blottner(1,0) = 4.49E-2;   Blottner(1,1) = -8.26E-2;  Blottner(1,2) = -9.20E0;  // O2
    Blottner(2,0) = 4.36E-2;   Blottner(2,1) = -3.36E-2;  Blottner(2,2) = -9.58E0;  // NO
    Blottner(3,0) = 1.16E-2;   Blottner(3,1) =  6.03E-1;  Blottner(3,2) = -1.24E1;  // N
    Blottner(4,0) = 2.03E-2;   Blottner(4,1) =  4.29E-1;  Blottner(4,2) = -1.16E1;  // O
    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 7;                     // O2
    nElStates[2] = 16;                    // NO
    nElStates[3] = 3;                     // N
    nElStates[4] = 5;                     // O
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);
    /*--- Allocate and initialize electron data arrays ---*/
    CharElTemp.resize(nSpecies,maxEl) = su2double(0.0);
    ElDegeneracy.resize(nSpecies,maxEl) = su2double(0.0);

    //N2: 15 states
    CharElTemp(0,0)  = 0.000000000000000E+00;
    CharElTemp(0,1)  = 7.223156514095200E+04;
    CharElTemp(0,2)  = 8.577862640384000E+04;
    CharElTemp(0,3)  = 8.605026716160000E+04;
    CharElTemp(0,4)  = 9.535118627874400E+04;
    CharElTemp(0,5)  = 9.805635702203200E+04;
    CharElTemp(0,6)  = 9.968267656935200E+04;
    CharElTemp(0,7)  = 1.048976467715200E+05;
    CharElTemp(0,8)  = 1.116489555200000E+05;
    CharElTemp(0,9)  = 1.225836470400000E+05;
    CharElTemp(0,10) = 1.248856873600000E+05;
    CharElTemp(0,11) = 1.282476158188320E+05;
    CharElTemp(0,12) = 1.338060936000000E+05;
    CharElTemp(0,13) = 1.404296391107200E+05;
    CharElTemp(0,14) = 1.504958859200000E+05;
    ElDegeneracy(0,0)  = 1;
    ElDegeneracy(0,1)  = 3;
    ElDegeneracy(0,2)  = 6;
    ElDegeneracy(0,3)  = 6;
    ElDegeneracy(0,4)  = 3;
    ElDegeneracy(0,5)  = 1;
    ElDegeneracy(0,6)  = 2;
    ElDegeneracy(0,7)  = 2;
    ElDegeneracy(0,8)  = 5;
    ElDegeneracy(0,9)  = 1;
    ElDegeneracy(0,10) = 6;
    ElDegeneracy(0,11) = 6;
    ElDegeneracy(0,12) = 10;
    ElDegeneracy(0,13) = 6;
    ElDegeneracy(0,14) = 6;
    // O2: 7 states
    CharElTemp(1,0) = 0.000000000000000E+00;
    CharElTemp(1,1) = 1.139156019700800E+04;
    CharElTemp(1,2) = 1.898473947826400E+04;
    CharElTemp(1,3) = 4.755973576639200E+04;
    CharElTemp(1,4) = 4.991242097343200E+04;
    CharElTemp(1,5) = 5.092268575561600E+04;
    CharElTemp(1,6) = 7.189863255967200E+04;
    ElDegeneracy(1,0) = 3;
    ElDegeneracy(1,1) = 2;
    ElDegeneracy(1,2) = 1;
    ElDegeneracy(1,3) = 1;
    ElDegeneracy(1,4) = 6;
    ElDegeneracy(1,5) = 3;
    ElDegeneracy(1,6) = 3;
    // NO: 16 states
    CharElTemp(2,0)  = 0.000000000000000E+00;
    CharElTemp(2,1)  = 5.467345760000000E+04;
    CharElTemp(2,2)  = 6.317139627802400E+04;
    CharElTemp(2,3)  = 6.599450342445600E+04;
    CharElTemp(2,4)  = 6.906120960000000E+04;
    CharElTemp(2,5)  = 7.049998480000000E+04;
    CharElTemp(2,6)  = 7.491055017560000E+04;
    CharElTemp(2,7)  = 7.628875293968000E+04;
    CharElTemp(2,8)  = 8.676188537552000E+04;
    CharElTemp(2,9)  = 8.714431182368000E+04;
    CharElTemp(2,10) = 8.886077063728000E+04;
    CharElTemp(2,11) = 8.981755614528000E+04;
    CharElTemp(2,12) = 8.988445919208000E+04;
    CharElTemp(2,13) = 9.042702132000000E+04;
    CharElTemp(2,14) = 9.064283760000000E+04;
    CharElTemp(2,15) = 9.111763341600000E+04;
    ElDegeneracy(2,0)  = 4;
    ElDegeneracy(2,1)  = 8;
    ElDegeneracy(2,2)  = 2;
    ElDegeneracy(2,3)  = 4;
    ElDegeneracy(2,4)  = 4;
    ElDegeneracy(2,5)  = 4;
    ElDegeneracy(2,6)  = 4;
    ElDegeneracy(2,7)  = 2;
    ElDegeneracy(2,8)  = 4;
    ElDegeneracy(2,9)  = 2;
    ElDegeneracy(2,10) = 4;
    ElDegeneracy(2,11) = 4;
    ElDegeneracy(2,12) = 2;
    ElDegeneracy(2,13) = 2;
    ElDegeneracy(2,14) = 2;
    ElDegeneracy(2,15) = 4;
    // N: 3 states
    CharElTemp(3,0) = 0.000000000000000E+00;
    CharElTemp(3,1) = 2.766469645581980E+04;
    CharElTemp(3,2) = 4.149309313560210E+04;
    ElDegeneracy(3,0)= 4;
    ElDegeneracy(3,1)= 10;
    ElDegeneracy(3,2)= 6;
    // O: 5 states
    CharElTemp(4,0) = 0.000000000000000E+00;
    CharElTemp(4,1) = 2.277077570280000E+02;
    CharElTemp(4,2) = 3.265688785704000E+02;
    CharElTemp(4,3) = 2.283028632262240E+04;
    CharElTemp(4,4) = 4.861993036434160E+04;
    ElDegeneracy(4,0) = 5;
    ElDegeneracy(4,1) = 3;
    ElDegeneracy(4,2) = 1;
    ElDegeneracy(4,3) = 5;
    ElDegeneracy(4,4) = 1;
    /*--- Set reaction maps ---*/
    // N2 dissociation
    Reactions(0,0,0)=0;    Reactions(0,0,1)=0;   Reactions(0,0,2)=nSpecies;    Reactions(0,1,0)=3;   Reactions(0,1,1)=3;   Reactions(0,1,2) =0;
    Reactions(1,0,0)=0;    Reactions(1,0,1)=1;   Reactions(1,0,2)=nSpecies;    Reactions(1,1,0)=3;   Reactions(1,1,1)=3;   Reactions(1,1,2) =1;
    Reactions(2,0,0)=0;    Reactions(2,0,1)=2;   Reactions(2,0,2)=nSpecies;    Reactions(2,1,0)=3;   Reactions(2,1,1)=3;   Reactions(2,1,2) =2;
    Reactions(3,0,0)=0;    Reactions(3,0,1)=3;   Reactions(3,0,2)=nSpecies;    Reactions(3,1,0)=3;   Reactions(3,1,1)=3;   Reactions(3,1,2) =3;
    Reactions(4,0,0)=0;    Reactions(4,0,1)=4;   Reactions(4,0,2)=nSpecies;    Reactions(4,1,0)=3;   Reactions(4,1,1)=3;   Reactions(4,1,2) =4;
    // O2 dissociation
    Reactions(5,0,0)=1;    Reactions(5,0,1)=0;   Reactions(5,0,2)=nSpecies;    Reactions(5,1,0)=4;   Reactions(5,1,1)=4;   Reactions(5,1,2) =0;
    Reactions(6,0,0)=1;    Reactions(6,0,1)=1;   Reactions(6,0,2)=nSpecies;    Reactions(6,1,0)=4;   Reactions(6,1,1)=4;   Reactions(6,1,2) =1;
    Reactions(7,0,0)=1;    Reactions(7,0,1)=2;   Reactions(7,0,2)=nSpecies;    Reactions(7,1,0)=4;   Reactions(7,1,1)=4;   Reactions(7,1,2) =2;
    Reactions(8,0,0)=1;    Reactions(8,0,1)=3;   Reactions(8,0,2)=nSpecies;    Reactions(8,1,0)=4;   Reactions(8,1,1)=4;   Reactions(8,1,2) =3;
    Reactions(9,0,0)=1;    Reactions(9,0,1)=4;   Reactions(9,0,2)=nSpecies;    Reactions(9,1,0)=4;   Reactions(9,1,1)=4;   Reactions(9,1,2) =4;
    // NO dissociation
    Reactions(10,0,0)=2;   Reactions(10,0,1)=0;  Reactions(10,0,2)=nSpecies;   Reactions(10,1,0)=3;  Reactions(10,1,1)=4;    Reactions(10,1,2) =0;
    Reactions(11,0,0)=2;   Reactions(11,0,1)=1;  Reactions(11,0,2)=nSpecies;   Reactions(11,1,0)=3;  Reactions(11,1,1)=4;    Reactions(11,1,2) =1;
    Reactions(12,0,0)=2;   Reactions(12,0,1)=2;  Reactions(12,0,2)=nSpecies;   Reactions(12,1,0)=3;  Reactions(12,1,1)=4;    Reactions(12,1,2) =2;
    Reactions(13,0,0)=2;   Reactions(13,0,1)=3;  Reactions(13,0,2)=nSpecies;   Reactions(13,1,0)=3;  Reactions(13,1,1)=4;    Reactions(13,1,2) =3;
    Reactions(14,0,0)=2;   Reactions(14,0,1)=4;  Reactions(14,0,2)=nSpecies;   Reactions(14,1,0)=3;  Reactions(14,1,1)=4;    Reactions(14,1,2) =4;
    // N2 + O -> NO + N
    Reactions(15,0,0)=0;   Reactions(15,0,1)=4;  Reactions(15,0,2)=nSpecies;   Reactions(15,1,0)=2;  Reactions(15,1,1)=3;    Reactions(15,1,2)= nSpecies;
    // NO + O -> O2 + N
    Reactions(16,0,0)=2;   Reactions(16,0,1)=4;  Reactions(16,0,2)=nSpecies;   Reactions(16,1,0)=1;  Reactions(16,1,1)=3;    Reactions(16,1,2)= nSpecies;
    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 7.0E21;
    ArrheniusCoefficient[2]  = 7.0E21;
    ArrheniusCoefficient[3]  = 3.0E22;
    ArrheniusCoefficient[4]  = 3.0E22;
    ArrheniusCoefficient[5]  = 2.0E21;
    ArrheniusCoefficient[6]  = 2.0E21;
    ArrheniusCoefficient[7]  = 2.0E21;
    ArrheniusCoefficient[8]  = 1.0E22;
    ArrheniusCoefficient[9]  = 1.0E22;
    ArrheniusCoefficient[10] = 5.0E15;
    ArrheniusCoefficient[11] = 5.0E15;
    ArrheniusCoefficient[12] = 5.0E15;
    ArrheniusCoefficient[13] = 1.1E17;
    ArrheniusCoefficient[14] = 1.1E17;
    ArrheniusCoefficient[15] = 6.4E17;
    ArrheniusCoefficient[16] = 8.4E12;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    ArrheniusEta[2]  = -1.60;
    ArrheniusEta[3]  = -1.60;
    ArrheniusEta[4]  = -1.60;
    ArrheniusEta[5]  = -1.50;
    ArrheniusEta[6]  = -1.50;
    ArrheniusEta[7]  = -1.50;
    ArrheniusEta[8]  = -1.50;
    ArrheniusEta[9]  = -1.50;
    ArrheniusEta[10] = 0.0;
    ArrheniusEta[11] = 0.0;
    ArrheniusEta[12] = 0.0;
    ArrheniusEta[13] = 0.0;
    ArrheniusEta[14] = 0.0;
    ArrheniusEta[15] = -1.0;
    ArrheniusEta[16] = 0.0;
    // Characteristic temperature
    ArrheniusTheta[0]  = 113200.0;
    ArrheniusTheta[1]  = 113200.0;
    ArrheniusTheta[2]  = 113200.0;
    ArrheniusTheta[3]  = 113200.0;
    ArrheniusTheta[4]  = 113200.0;
    ArrheniusTheta[5]  = 59500.0;
    ArrheniusTheta[6]  = 59500.0;
    ArrheniusTheta[7]  = 59500.0;
    ArrheniusTheta[8]  = 59500.0;
    ArrheniusTheta[9]  = 59500.0;
    ArrheniusTheta[10] = 75500.0;
    ArrheniusTheta[11] = 75500.0;
    ArrheniusTheta[12] = 75500.0;
    ArrheniusTheta[13] = 75500.0;
    ArrheniusTheta[14] = 75500.0;
    ArrheniusTheta[15] = 38400.0;
    ArrheniusTheta[16] = 19450.0;
    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;  Tcb_b[0] = 0;
    Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;  Tcb_b[1] = 0;
    Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;  Tcb_b[2] = 0;
    Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;  Tcb_b[3] = 0;
    Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;  Tcb_b[4] = 0;
    Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;  Tcb_b[5] = 0;
    Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;  Tcb_b[6] = 0;
    Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;  Tcb_b[7] = 0;
    Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;  Tcb_b[8] = 0;
    Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;  Tcb_b[9] = 0;
    Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;  Tcb_b[10] = 0;
    Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;  Tcb_b[11] = 0;
    Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;  Tcb_b[12] = 0;
    Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;  Tcb_b[13] = 0;
    Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;  Tcb_b[14] = 0;
    Tcf_a[15] = 1.0; Tcf_b[15] = 0.0; Tcb_a[15] = 1;  Tcb_b[15] = 0;
    Tcf_a[16] = 1.0; Tcf_b[16] = 0.0; Tcb_a[16] = 1;  Tcb_b[16] = 0;
    /*--- Collision integral data ---*/
    // Omega(0,0) ----------------------
    //N2
    Omega00(0,0,0) = -6.0614558E-03;  Omega00(0,0,1) = 1.2689102E-01;   Omega00(0,0,2) = -1.0616948E+00;  Omega00(0,0,3) = 8.0955466E+02;
    Omega00(0,1,0) = -3.7959091E-03;  Omega00(0,1,1) = 9.5708295E-02;   Omega00(0,1,2) = -1.0070611E+00;  Omega00(0,1,3) = 8.9392313E+02;
    Omega00(0,2,0) = -1.9295666E-03;  Omega00(0,2,1) = 2.7995735E-02;   Omega00(0,2,2) = -3.1588514E-01;  Omega00(0,2,3) = 1.2880734E+02;
    Omega00(0,3,0) = -1.0796249E-02;  Omega00(0,3,1) = 2.2656509E-01;   Omega00(0,3,2) = -1.7910602E+00;  Omega00(0,3,3) = 4.0455218E+03;
    Omega00(0,4,0) = -2.7244269E-03;  Omega00(0,4,1) = 6.9587171E-02;   Omega00(0,4,2) = -7.9538667E-01;  Omega00(0,4,3) = 4.0673730E+02;
    //O2
    Omega00(1,0,0) = -3.7959091E-03;  Omega00(1,0,1) = 9.5708295E-02;   Omega00(1,0,2) = -1.0070611E+00;  Omega00(1,0,3) = 8.9392313E+02;
    Omega00(1,1,0) = -8.0682650E-04;  Omega00(1,1,1) = 1.6602480E-02;   Omega00(1,1,2) = -3.1472774E-01;  Omega00(1,1,3) = 1.4116458E+02;
    Omega00(1,2,0) = -6.4433840E-04;  Omega00(1,2,1) = 8.5378580E-03;   Omega00(1,2,2) = -2.3225102E-01;  Omega00(1,2,3) = 1.1371608E+02;
    Omega00(1,3,0) = -1.1453028E-03;  Omega00(1,3,1) = 1.2654140E-02;   Omega00(1,3,2) = -2.2435218E-01;  Omega00(1,3,3) = 7.7201588E+01;
    Omega00(1,4,0) = -4.8405803E-03;  Omega00(1,4,1) = 1.0297688E-01;   Omega00(1,4,2) = -9.6876576E-01;  Omega00(1,4,3) = 6.1629812E+02;
    //NO
    Omega00(2,0,0) = -1.9295666E-03;  Omega00(2,0,1) = 2.7995735E-02;   Omega00(2,0,2) = -3.1588514E-01;  Omega00(2,0,3) = 1.2880734E+02;
    Omega00(2,1,0) = -6.4433840E-04;  Omega00(2,1,1) = 8.5378580E-03;   Omega00(2,1,2) = -2.3225102E-01;  Omega00(2,1,3) = 1.1371608E+02;
    Omega00(2,2,0) = -0.0000000E+00;  Omega00(2,2,1) = -1.1056066E-02;  Omega00(2,2,2) = -5.9216250E-02;  Omega00(2,2,3) = 7.2542367E+01;
    Omega00(2,3,0) = -1.5770918E-03;  Omega00(2,3,1) = 1.9578381E-02;   Omega00(2,3,2) = -2.7873624E-01;  Omega00(2,3,3) = 9.9547944E+01;
    Omega00(2,4,0) = -1.0885815E-03;  Omega00(2,4,1) = 1.1883688E-02;   Omega00(2,4,2) = -2.1844909E-01;  Omega00(2,4,3) = 7.5512560E+01;
    //N
    Omega00(3,0,0) = -1.0796249E-02;  Omega00(3,0,1) = 2.2656509E-01;   Omega00(3,0,2) = -1.7910602E+00;  Omega00(3,0,3) = 4.0455218E+03;
    Omega00(3,1,0) = -1.1453028E-03;  Omega00(3,1,1) = 1.2654140E-02;   Omega00(3,1,2) = -2.2435218E-01;  Omega00(3,1,3) = 7.7201588E+01;
    Omega00(3,2,0) = -1.5770918E-03;  Omega00(3,2,1) = 1.9578381E-02;   Omega00(3,2,2) = -2.7873624E-01;  Omega00(3,2,3) = 9.9547944E+01;
    Omega00(3,3,0) = -9.6083779E-03;  Omega00(3,3,1) = 2.0938971E-01;   Omega00(3,3,2) = -1.7386904E+00;  Omega00(3,3,3) = 3.3587983E+03;
    Omega00(3,4,0) = -7.8147689E-03;  Omega00(3,4,1) = 1.6792705E-01;   Omega00(3,4,2) = -1.4308628E+00;  Omega00(3,4,3) = 1.6628859E+03;
    //O
    Omega00(4,0,0) = -2.7244269E-03;  Omega00(4,0,1) = 6.9587171E-02;   Omega00(4,0,2) = -7.9538667E-01;  Omega00(4,0,3) = 4.0673730E+02;
    Omega00(4,1,0) = -4.8405803E-03;  Omega00(4,1,1) = 1.0297688E-01;   Omega00(4,1,2) = -9.6876576E-01;  Omega00(4,1,3) = 6.1629812E+02;
    Omega00(4,2,0) = -1.0885815E-03;  Omega00(4,2,1) = 1.1883688E-02;   Omega00(4,2,2) = -2.1844909E-01;  Omega00(4,2,3) = 7.5512560E+01;
    Omega00(4,3,0) = -7.8147689E-03;  Omega00(4,3,1) = 1.6792705E-01;   Omega00(4,3,2) = -1.4308628E+00;  Omega00(4,3,3) = 1.6628859E+03;
    Omega00(4,4,0) = -6.4040535E-03;  Omega00(4,4,1) = 1.4629949E-01;   Omega00(4,4,2) = -1.3892121E+00;  Omega00(4,4,3) = 2.0903441E+03;

    // Omega(1,1) ----------------------
    //N2
    Omega11(0,0,0) = -7.6303990E-03;  Omega11(0,0,1) = 1.6878089E-01;   Omega11(0,0,2) = -1.4004234E+00;  Omega11(0,0,3) = 2.1427708E+03;
    Omega11(0,1,0) = -8.0457321E-03;  Omega11(0,1,1) = 1.9228905E-01;   Omega11(0,1,2) = -1.7102854E+00;  Omega11(0,1,3) = 5.2213857E+03;
    Omega11(0,2,0) = -6.8237776E-03;  Omega11(0,2,1) = 1.4360616E-01;   Omega11(0,2,2) = -1.1922240E+00;  Omega11(0,2,3) = 1.2433086E+03;
    Omega11(0,3,0) = -8.3493693E-03;  Omega11(0,3,1) = 1.7808911E-01;   Omega11(0,3,2) = -1.4466155E+00;  Omega11(0,3,3) = 1.9324210E+03;
    Omega11(0,4,0) = -8.3110691E-03;  Omega11(0,4,1) = 1.9617877E-01;   Omega11(0,4,2) = -1.7205427E+00;  Omega11(0,4,3) = 4.0812829E+03;
    //O2
    Omega11(1,0,0) = -8.0457321E-03;  Omega11(1,0,1) = 1.9228905E-01;   Omega11(1,0,2) = -1.7102854E+00;  Omega11(1,0,3) = 5.2213857E+03;
    Omega11(1,1,0) = -6.2931612E-03;  Omega11(1,1,1) = 1.4624645E-01;   Omega11(1,1,2) = -1.3006927E+00;  Omega11(1,1,3) = 1.8066892E+03;
    Omega11(1,2,0) = -6.8508672E-03;  Omega11(1,2,1) = 1.5524564E-01;   Omega11(1,2,2) = -1.3479583E+00;  Omega11(1,2,3) = 2.0037890E+03;
    Omega11(1,3,0) = -1.0608832E-03;  Omega11(1,3,1) = 1.1782595E-02;   Omega11(1,3,2) = -2.1246301E-01;  Omega11(1,3,3) = 8.4561598E+01;
    Omega11(1,4,0) = -3.7969686E-03;  Omega11(1,4,1) = 7.6789981E-02;   Omega11(1,4,2) = -7.3056809E-01;  Omega11(1,4,3) = 3.3958171E+02;
    //NO
    Omega11(2,0,0) = -6.8237776E-03;  Omega11(2,0,1) = 1.4360616E-01;   Omega11(2,0,2) = -1.1922240E+00;  Omega11(2,0,3) = 1.2433086E+03;
    Omega11(2,1,0) = -6.8508672E-03;  Omega11(2,1,1) = 1.5524564E-01;   Omega11(2,1,2) = -1.3479583E+00;  Omega11(2,1,3) = 2.0037890E+03;
    Omega11(2,2,0) = -7.4942466E-03;  Omega11(2,2,1) = 1.6626193E-01;   Omega11(2,2,2) = -1.4107027E+00;  Omega11(2,2,3) = 2.3097604E+03;
    Omega11(2,3,0) = -1.4719259E-03;  Omega11(2,3,1) = 1.8446968E-02;   Omega11(2,3,2) = -2.6460411E-01;  Omega11(2,3,3) = 1.0911124E+02;
    Omega11(2,4,0) = -1.0066279E-03;  Omega11(2,4,1) = 1.1029264E-02;   Omega11(2,4,2) = -2.0671266E-01;  Omega11(2,4,3) = 8.2644384E+01;
    //N
    Omega11(3,0,0) = -8.3493693E-03;  Omega11(3,0,1) = 1.7808911E-01;   Omega11(3,0,2) = -1.4466155E+00;  Omega11(3,0,3) = 1.9324210E+03;
    Omega11(3,1,0) = -1.0608832E-03;  Omega11(3,1,1) = 1.1782595E-02;   Omega11(3,1,2) = -2.1246301E-01;  Omega11(3,1,3) = 8.4561598E+01;
    Omega11(3,2,0) = -1.4719259E-03;  Omega11(3,2,1) = 1.8446968E-02;   Omega11(3,2,2) = -2.6460411E-01;  Omega11(3,2,3) = 1.0911124E+02;
    Omega11(3,3,0) = -7.7439615E-03;  Omega11(3,3,1) = 1.7129007E-01;   Omega11(3,3,2) = -1.4809088E+00;  Omega11(3,3,3) = 2.1284951E+03;
    Omega11(3,4,0) = -5.0478143E-03;  Omega11(3,4,1) = 1.0236186E-01;   Omega11(3,4,2) = -9.0058935E-01;  Omega11(3,4,3) = 4.4472565E+02;
    //O
    Omega11(4,0,0) = -8.3110691E-03;  Omega11(4,0,1) = 1.9617877E-01;   Omega11(4,0,2) = -1.7205427E+00;  Omega11(4,0,3) = 4.0812829E+03;
    Omega11(4,1,0) = -3.7969686E-03;  Omega11(4,1,1) = 7.6789981E-02;   Omega11(4,1,2) = -7.3056809E-01;  Omega11(4,1,3) = 3.3958171E+02;
    Omega11(4,2,0) = -1.0066279E-03;  Omega11(4,2,1) = 1.1029264E-02;   Omega11(4,2,2) = -2.0671266E-01;  Omega11(4,2,3) = 8.2644384E+01;
    Omega11(4,3,0) = -5.0478143E-03;  Omega11(4,3,1) = 1.0236186E-01;   Omega11(4,3,2) = -9.0058935E-01;  Omega11(4,3,3) = 4.4472565E+02;
    Omega11(4,4,0) = -4.2451096E-03;  Omega11(4,4,1) = 9.6820337E-02;   Omega11(4,4,2) = -9.9770795E-01;  Omega11(4,4,3) = 8.3320644E+02;
  }

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }
}

CSU2TCLib::~CSU2TCLib(){}

void CSU2TCLib::SetTDStateRhosTTv(vector<su2double>& val_rhos, su2double val_temperature, su2double val_temperature_ve){

  rhos = val_rhos;
  T    = val_temperature;
  Tve  = val_temperature_ve;

  Density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Density += rhos[iSpecies];

  Pressure = ComputePressure();

}

vector<su2double>& CSU2TCLib::GetSpeciesCvTraRot(){

  // NOTE: This will be non-const with future models.
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    Cvtrs[iSpecies] = (3.0/2.0 + RotationModes[iSpecies]/2.0) * Ru/MolarMass[iSpecies];

  return Cvtrs;
}

vector<su2double>& CSU2TCLib::ComputeSpeciesCvVibEle(){

  su2double thoTve, exptv, num, num2, num3, denom, Cvvs, Cves;
  unsigned short iElectron = nSpecies-1;

  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    /*--- If requesting electron specific heat ---*/
    if (ionization && iSpecies == iElectron) {
      Cvvs = 0.0;
      Cves = 3.0/2.0 * Ru/MolarMass[nSpecies-1];
    }

    /*--- Heavy particle specific heat ---*/
    else {

      /*--- Vibrational energy ---*/
      if (CharVibTemp[iSpecies] != 0.0) {
        thoTve = CharVibTemp[iSpecies]/Tve;
        exptv = exp(CharVibTemp[iSpecies]/Tve);
        Cvvs  = Ru/MolarMass[iSpecies] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
      } else {
        Cvvs = 0.0;
      }

      /*--- Electronic energy ---*/
      if (nElStates[iSpecies] != 0) {
        num = 0.0; num2 = 0.0;
        denom = ElDegeneracy[iSpecies][0] * exp(-CharElTemp[iSpecies][0]/Tve);
        num3  = ElDegeneracy[iSpecies][0] * (CharElTemp[iSpecies][0]/(Tve*Tve))*exp(-CharElTemp[iSpecies][0]/Tve);
        for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
          thoTve = CharElTemp[iSpecies][iEl]/Tve;
          exptv = exp(-CharElTemp[iSpecies][iEl]/Tve);

          num   += ElDegeneracy[iSpecies][iEl] * CharElTemp[iSpecies][iEl] * exptv;
          denom += ElDegeneracy[iSpecies][iEl] * exptv;
          num2  += ElDegeneracy[iSpecies][iEl] * (thoTve*thoTve) * exptv;
          num3  += ElDegeneracy[iSpecies][iEl] * thoTve/Tve * exptv;
        }
        Cves = Ru/MolarMass[iSpecies] * (num2/denom - num*num3/(denom*denom));
      } else {
        Cves = 0.0;
      }
    }

    Cvves[iSpecies] = Cvvs + Cves;
  }

  return Cvves;

}

vector<su2double>& CSU2TCLib::ComputeMixtureEnergies(){

  su2double rhoEmix, rhoEve, Ef, Ev, Ee, num, denom;

  rhoEmix = 0.0;
  rhoEve  = 0.0;
  denom   = 0.0;

  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++){

    // Species formation energy
    Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies];

    // Species vibrational energy
    if (CharVibTemp[iSpecies] != 0.0)
      Ev = Ru/MolarMass[iSpecies] * CharVibTemp[iSpecies] / (exp(CharVibTemp[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = ElDegeneracy(iSpecies,0) * exp(CharElTemp(iSpecies,0)/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += ElDegeneracy(iSpecies,iEl) * CharElTemp(iSpecies,iEl) * exp(-CharElTemp(iSpecies,iEl)/Tve);
      denom += ElDegeneracy(iSpecies,iEl) * exp(-CharElTemp(iSpecies,iEl)/Tve);
    }
    Ee = Ru/MolarMass[iSpecies] * (num/denom);

    // Mixture total energy
    rhoEmix += rhos[iSpecies] * ((3.0/2.0+RotationModes[iSpecies]/2.0) * Ru/MolarMass[iSpecies] * (T-Ref_Temperature[iSpecies]) + Ev + Ee + Ef);

    // Mixture vibrational-electronic energy
    rhoEve += rhos[iSpecies] * (Ev + Ee);

  }

  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {

    // Species formation energy
    Ef = Enthalpy_Formation[nSpecies-1] - Ru/MolarMass[nSpecies-1] * Ref_Temperature[nSpecies-1];

    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += (3.0/2.0) * Ru/MolarMass[nSpecies-1] * (Tve - Ref_Temperature[nSpecies-1]); //bug? not multiplying by rhos[iSpecies]
  }

  energies[0] = rhoEmix/Density;
  energies[1] = rhoEve/Density;

  return energies;

}

vector<su2double>& CSU2TCLib::ComputeSpeciesEve(su2double val_T){

  su2double Ev, Eel, Ef, num, denom;
  unsigned short iElectron = nSpecies-1;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){

    /*--- Electron species energy ---*/
    if ( ionization && (iSpecies == iElectron)) {
      /*--- Calculate formation energy ---*/
      Ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies] * Ref_Temperature[iSpecies];

      /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
      Eel = (3.0/2.0) * Ru/MolarMass[iSpecies] * (val_T - Ref_Temperature[iSpecies]) + Ef;
      Ev  = 0.0;
    }
    /*--- Heavy particle energy ---*/
    else {
      /*--- Calculate vibrational energy (harmonic-oscillator model) ---*/
      if (CharVibTemp[iSpecies] != 0.0)
        Ev = Ru/MolarMass[iSpecies] * CharVibTemp[iSpecies] / (exp(CharVibTemp[iSpecies]/val_T)-1.0);
      else
        Ev = 0.0;
      /*--- Calculate electronic energy ---*/
      num = 0.0;
      denom = ElDegeneracy[iSpecies][0] * exp(-CharElTemp[iSpecies][0]/val_T);
      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        num   += ElDegeneracy[iSpecies][iEl] * CharElTemp[iSpecies][iEl] * exp(-CharElTemp[iSpecies][iEl]/val_T);
        denom += ElDegeneracy[iSpecies][iEl] * exp(-CharElTemp[iSpecies][iEl]/val_T);
      }
      Eel = Ru/MolarMass[iSpecies] * (num/denom);
    }

    eves[iSpecies] = Ev + Eel;
  }

  return eves;
}

vector<su2double>& CSU2TCLib::ComputeNetProductionRates(){

  /*--- Nonequilibrium chemistry ---*/
  unsigned short ii, iReaction;
  su2double T_min, epsilon, Thf, Thb, Trxnf, Trxnb, Keq, kf, kb, kfb, fwdRxn, bkwRxn, af, bf, ab, bb;

  /*--- Define artificial chemistry parameters ---*/
  // Note: These parameters artificially increase the rate-controlling reaction
  //       temperature.  This relaxes some of the stiffness in the chemistry
  //       source term.
  T_min   = 800.0;
  epsilon = 80;
  /*--- Define preferential dissociation coefficient ---*/
  //alpha = 0.3;

  for( iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    ws[iSpecies] = 0.0;

  for (iReaction = 0; iReaction < nReactions; iReaction++) {

    /*--- Determine the rate-controlling temperature ---*/
    af = Tcf_a[iReaction];
    bf = Tcf_b[iReaction];
    ab = Tcb_a[iReaction];
    bb = Tcb_b[iReaction];
    Trxnf = pow(T, af)*pow(Tve, bf);
    Trxnb = pow(T, ab)*pow(Tve, bb);

    /*--- Calculate the modified temperature ---*/
    Thf = 0.5 * (Trxnf+T_min + sqrt((Trxnf-T_min)*(Trxnf-T_min)+epsilon*epsilon));
    Thb = 0.5 * (Trxnb+T_min + sqrt((Trxnb-T_min)*(Trxnb-T_min)+epsilon*epsilon));

    /*--- Get the Keq & Arrhenius coefficients ---*/
    ComputeKeqConstants(iReaction);

    /*--- Calculate Keq ---*/
    Keq = exp(  A[0]*(Thb/1E4) + A[1] + A[2]*log(1E4/Thb)
        + A[3]*(1E4/Thb) + A[4]*(1E4/Thb)*(1E4/Thb) );

    /*--- Calculate rate coefficients ---*/
    kf  = ArrheniusCoefficient[iReaction] * exp(ArrheniusEta[iReaction]*log(Thf)) * exp(-ArrheniusTheta[iReaction]/Thf);
    kfb = ArrheniusCoefficient[iReaction] * exp(ArrheniusEta[iReaction]*log(Thb)) * exp(-ArrheniusTheta[iReaction]/Thb);
    kb  = kfb / Keq;

    /*--- Determine production & destruction of each species ---*/
    fwdRxn = 1.0;
    bkwRxn = 1.0;
    for (ii = 0; ii < 3; ii++) {

      /*--- Reactants ---*/
      iSpecies = Reactions(iReaction,0,ii);
      if ( iSpecies != nSpecies)
        fwdRxn *= 0.001*rhos[iSpecies]/MolarMass[iSpecies];

      /*--- Products ---*/
      jSpecies = Reactions(iReaction,1,ii);
      if (jSpecies != nSpecies) {
        bkwRxn *= 0.001*rhos[jSpecies]/MolarMass[jSpecies];
      }
    }

    fwdRxn = 1000.0 * kf * fwdRxn;
    bkwRxn = 1000.0 * kb * bkwRxn;

    for (ii = 0; ii < 3; ii++) {

      /*--- Products ---*/
      iSpecies = Reactions(iReaction,1,ii);
      if (iSpecies != nSpecies)
        ws[iSpecies] += MolarMass[iSpecies] * (fwdRxn-bkwRxn);

      /*--- Reactants ---*/
      iSpecies = Reactions(iReaction,0,ii);
      if (iSpecies != nSpecies)
        ws[iSpecies] -= MolarMass[iSpecies] * (fwdRxn-bkwRxn);
    }
  }
  return ws;
}

void CSU2TCLib::ComputeKeqConstants(unsigned short val_Reaction) {

  unsigned short ii, iIndex, tbl_offset, pwr;
  su2double N, tmp1, tmp2;

  /*--- Acquire database constants from CConfig ---*/
  GetChemistryEquilConstants(val_Reaction);

  /*--- Calculate mixture number density ---*/
  N = 0.0;
  for (iSpecies =0 ; iSpecies < nSpecies; iSpecies++) {
    N += rhos[iSpecies]/MolarMass[iSpecies]*AVOGAD_CONSTANT;
  }

  /*--- Convert number density from 1/m^3 to 1/cm^3 for table look-up ---*/
  N = N*(1E-6);

  /*--- Determine table index based on mixture N ---*/
  tbl_offset = 14;
  pwr        = floor(log10(N));

  /*--- Bound the interpolation to table limit values ---*/
  iIndex = int(pwr) - tbl_offset;
  if (iIndex <= 0) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable(0,ii);
    return;
  } else if (iIndex >= 5) {
    for (ii = 0; ii < 5; ii++)
      A[ii] = RxnConstantTable(5,ii);
    return;
  }

  /*--- Calculate interpolation denominator terms avoiding pow() ---*/
  tmp1 = 1.0;
  tmp2 = 1.0;
  for (ii = 0; ii < pwr; ii++) {
    tmp1 *= 10.0;
    tmp2 *= 10.0;
  }
  tmp2 *= 10.0;

  /*--- Interpolate ---*/
  for (ii = 0; ii < 5; ii++) {
    A[ii] =  (RxnConstantTable(iIndex+1,ii) - RxnConstantTable(iIndex,ii))
        / (tmp2 - tmp1) * (N - tmp1)
        + RxnConstantTable(iIndex,ii);
  }
}

su2double CSU2TCLib::ComputeEveSourceTerm(){

  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section
  su2double conc, N, mu, A_sr, B_sr, num, denom, Cs, sig_s, tau_sr, tauP, tauMW, taus, omegaVT, omegaCV;
  vector<su2double> MolarFrac, eve_eq, eve;

  MolarFrac.resize(nSpecies,0.0);
  eve_eq.resize(nSpecies,0.0);
  eve.resize(nSpecies,0.0);

  omegaVT = 0.0;
  omegaCV = 0.0;


  /*--- Calculate mole fractions ---*/
  N    = 0.0;
  conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    conc += rhos[iSpecies] / MolarMass[iSpecies];
    N    += rhos[iSpecies] / MolarMass[iSpecies] * AVOGAD_CONSTANT;
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFrac[iSpecies] = (rhos[iSpecies] / MolarMass[iSpecies]) / conc;

  eve_eq = ComputeSpeciesEve(T);
  eve    = ComputeSpeciesEve(Tve);

  /*--- Loop over species to calculate source term --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Millikan & White relaxation time ---*/
    num   = 0.0;
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      mu     = MolarMass[iSpecies]*MolarMass[jSpecies] / (MolarMass[iSpecies] + MolarMass[jSpecies]);
      A_sr   = 1.16 * 1E-3 * sqrt(mu) * pow(CharVibTemp[iSpecies], 4.0/3.0);
      B_sr   = 0.015 * pow(mu, 0.25);
      tau_sr = 101325.0/Pressure * exp(A_sr*(pow(T,-1.0/3.0) - B_sr) - 18.42);
      num   += MolarFrac[jSpecies];
      denom += MolarFrac[jSpecies] / tau_sr;
    }

    tauMW = num / denom;

    /*--- Park limiting cross section ---*/
    Cs    = sqrt((8.0*Ru*T)/(PI_NUMBER*MolarMass[iSpecies]));
    sig_s = 1E-20*(5E4*5E4)/(T*T);

    tauP = 1/(sig_s*Cs*N);

    /*--- Species relaxation time ---*/
    taus = tauMW + tauP;

    /*--- Add species contribution to residual ---*/
    omegaVT += rhos[iSpecies] * (eve_eq[iSpecies] -
                                 eve[iSpecies]) / taus;
  }

  /*--- Vibrational energy change due to chemical reactions ---*/
  if(!frozen){
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      omegaCV += ws[iSpecies]*eve[iSpecies];
  }

  omega = omegaVT + omegaCV;

  return omega;

}

vector<su2double>& CSU2TCLib::ComputeSpeciesEnthalpy(su2double val_T, su2double val_Tve, su2double *val_eves){

  vector<su2double> cvtrs;

  cvtrs = GetSpeciesCvTraRot();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    eves[iSpecies] = val_eves[iSpecies];
    hs[iSpecies] = Ru/MolarMass[iSpecies]*val_T + cvtrs[iSpecies]*val_T + Enthalpy_Formation[iSpecies] + eves[iSpecies];
  }

  return hs;

}

vector<su2double>& CSU2TCLib::GetDiffusionCoeff(){

  if(Kind_TransCoeffModel == WILKE)
   DiffusionCoeffWBE();
  if(Kind_TransCoeffModel == GUPTAYOS)
   DiffusionCoeffGY();

  return DiffusionCoeff;

}

su2double CSU2TCLib::GetViscosity(){

  if(Kind_TransCoeffModel == WILKE)
    ViscosityWBE();
  if(Kind_TransCoeffModel == GUPTAYOS)
    ViscosityGY();

  return Mu;

}

vector<su2double>& CSU2TCLib::GetThermalConductivities(){

  if(Kind_TransCoeffModel == WILKE)
    ThermalConductivitiesWBE();
  if(Kind_TransCoeffModel == GUPTAYOS)
    ThermalConductivitiesGY();

  return ThermalConductivities;

}

void CSU2TCLib::DiffusionCoeffWBE(){

  su2double conc, Mi, Mj, M, Omega_ij, denom;
  su2activematrix Dij;

  Dij.resize(nSpecies, nSpecies) = su2double(0.0);

  /*--- Calculate species mole fraction ---*/
  conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    MolarFracWBE[iSpecies] = rhos[iSpecies]/MolarMass[iSpecies];
    conc               += MolarFracWBE[iSpecies];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFracWBE[iSpecies] = MolarFracWBE[iSpecies]/conc;
  /*--- Calculate mixture molar mass (kg/mol) ---*/
  // Note: Species molar masses stored as kg/kmol, need 1E-3 conversion
  M = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    M += MolarMass[iSpecies]*MolarFracWBE[iSpecies];
  M = M*1E-3;
  /*---+++                  +++---*/
  /*--- Diffusion coefficients ---*/
  /*---+++                  +++---*/
  /*--- Solve for binary diffusion coefficients ---*/
  // Note: Dij = Dji, so only loop through req'd indices
  // Note: Correlation requires kg/mol, hence 1E-3 conversion from kg/kmol
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Mi = MolarMass[iSpecies]*1E-3;
    for (jSpecies = iSpecies; jSpecies < nSpecies; jSpecies++) {
      Mj = MolarMass[jSpecies]*1E-3;
      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20/PI_NUMBER * Omega00(iSpecies,jSpecies,3)
          * pow(T, Omega00(iSpecies,jSpecies,0)*log(T)*log(T)
          +  Omega00(iSpecies,jSpecies,1)*log(T)
          +  Omega00(iSpecies,jSpecies,2));
      Dij(iSpecies,jSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(Density*Omega_ij);
      Dij(jSpecies,iSpecies) = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(Density*Omega_ij);
    }
  }
  /*--- Calculate species-mixture diffusion coefficient --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    DiffusionCoeff[iSpecies] = 0.0;
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      if (jSpecies != iSpecies) {
        denom += MolarFracWBE[jSpecies]/Dij(iSpecies,jSpecies);
      }
    }

    if (nSpecies==1) DiffusionCoeff[0] = 0;
    else DiffusionCoeff[iSpecies] = (1-MolarFracWBE[iSpecies])/denom;
  }
}

void CSU2TCLib::ViscosityWBE(){

  su2double tmp1, tmp2, conc;

  /*--- Calculate species mole fraction ---*/
  conc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    MolarFracWBE[iSpecies] = rhos[iSpecies]/MolarMass[iSpecies];
    conc               += MolarFracWBE[iSpecies];
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MolarFracWBE[iSpecies] = MolarFracWBE[iSpecies]/conc;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
                             Blottner[iSpecies][1])*log(T) +
                             Blottner[iSpecies][2]);

  /*--- Determine species 'phi' value for Blottner model ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    phis[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(MolarMass[jSpecies]/MolarMass[iSpecies], 0.25);
      tmp2 = sqrt(8.0*(1.0+MolarMass[iSpecies]/MolarMass[jSpecies]));
      phis[iSpecies] += MolarFracWBE[jSpecies]*tmp1*tmp1/tmp2;
    }
  }

  /*--- Calculate mixture laminar viscosity ---*/
  Mu = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    Mu += MolarFracWBE[iSpecies]*mus[iSpecies]/phis[iSpecies];
  }
}

void CSU2TCLib::ThermalConductivitiesWBE(){

  vector<su2double> ks, kves;

  ks.resize(nSpecies,0.0);
  kves.resize(nSpecies,0.0);


  Cvves = ComputeSpeciesCvVibEle();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ks[iSpecies] = mus[iSpecies]*(15.0/4.0 + RotationModes[iSpecies]/2.0)*Ru/MolarMass[iSpecies];
    kves[iSpecies] = mus[iSpecies]*Cvves[iSpecies];
  }
  /*--- Calculate mixture tr & ve conductivities ---*/
  ThermalCond_tr = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ThermalCond_tr += MolarFracWBE[iSpecies]*ks[iSpecies]/phis[iSpecies];
    ThermalCond_ve += MolarFracWBE[iSpecies]*kves[iSpecies]/phis[iSpecies];
  }

  ThermalConductivities[0] = ThermalCond_tr;
  ThermalConductivities[1] = ThermalCond_ve;
}

void CSU2TCLib::DiffusionCoeffGY(){

  su2double Mi, Mj, pi, kb, gam_i, gam_j, gam_t, denom, d1_ij, D_ij, Omega_ij;

  pi   = PI_NUMBER;
  kb   = BOLTZMANN_CONSTANT;

  /*--- Calculate mixture gas constant ---*/
  gam_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    gam_t += rhos[iSpecies] / (Density*MolarMass[iSpecies]);
  }
  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;
    /*--- Calculate molar concentration ---*/
    Mi      = MolarMass[iSpecies];
    gam_i   = rhos[iSpecies] / (Density*Mi);
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      if (jSpecies != iSpecies) {
        Mj    = MolarMass[jSpecies];
        gam_j = rhos[iSpecies] / (Density*Mj);
        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00(iSpecies,jSpecies,3)
            * pow(T, Omega00(iSpecies,jSpecies,0)*log(T)*log(T)
            + Omega00(iSpecies,jSpecies,1)*log(T)
            + Omega00(iSpecies,jSpecies,2));
        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
        /*--- Calculate heavy-particle binary diffusion coefficient ---*/
        D_ij = kb*T/(Pressure*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    if (ionization) {
      jSpecies = nSpecies-1;
      Mj       = MolarMass[jSpecies];
      gam_j    = rhos[iSpecies] / (Density*Mj);

      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega00(iSpecies,jSpecies,3)
          * pow(Tve, Omega00(iSpecies,jSpecies,0)*log(Tve)*log(Tve)
          + Omega00(iSpecies,jSpecies,1)*log(Tve)
          + Omega00(iSpecies,jSpecies,2));

      /*--- Calculate "delta1_ij" ---*/
      d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
    }

    /*--- Assign species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = gam_t*gam_t*Mi*(1-Mi*gam_i) / denom;
  }
  if (ionization) {
    iSpecies = nSpecies-1;

    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;

    /*--- Calculate molar concentration ---*/
    Mi      = MolarMass[iSpecies];
    gam_i   = rhos[iSpecies] / (Density*Mi);
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      if (iSpecies != jSpecies) {
        Mj    = MolarMass[jSpecies];
        gam_j = rhos[iSpecies] / (Density*Mj);

        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00(iSpecies,jSpecies,3)
            * pow(Tve, Omega00(iSpecies,jSpecies,0)*log(Tve)*log(Tve)
            + Omega00(iSpecies,jSpecies,1)*log(Tve)
            + Omega00(iSpecies,jSpecies,2));

        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

        /*--- Calculate heavy-particle binary diffusion coefficient ---*/
        D_ij = kb*Tve/(Pressure*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    DiffusionCoeff[iSpecies] = gam_t*gam_t*MolarMass[iSpecies]*(1-MolarMass[iSpecies]*gam_i) / denom;
  }
}

void CSU2TCLib::ViscosityGY(){

  su2double Mi, Mj, pi, Na, gam_i, gam_j, denom, Omega_ij, d2_ij;

  pi   = PI_NUMBER;
  Na   = AVOGAD_CONSTANT;
  Mu = 0.0;
  /*--- Mixture viscosity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    denom = 0.0;
    /*--- Calculate molar concentration ---*/
    Mi    = MolarMass[iSpecies];
    gam_i = rhos[iSpecies] / (Density*Mi);
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      Mj    = MolarMass[jSpecies];
      gam_j = rhos[jSpecies] / (Density*Mj);
      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11(iSpecies,jSpecies,3)
          * pow(T, Omega11(iSpecies,jSpecies,0)*log(T)*log(T)
          + Omega11(iSpecies,jSpecies,1)*log(T)
          + Omega11(iSpecies,jSpecies,2));
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
      /*--- Add to denominator of viscosity ---*/
      denom += gam_j*d2_ij;
    }
    if (ionization) {
      jSpecies = nSpecies-1;
      Mj    = MolarMass[jSpecies];
      gam_j = rhos[jSpecies] / (Density*Mj);
      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11(iSpecies,jSpecies,3)
          * pow(Tve, Omega11(iSpecies,jSpecies,0)*log(Tve)*log(Tve)
          + Omega11(iSpecies,jSpecies,1)*log(Tve)
          + Omega11(iSpecies,jSpecies,2));
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
      denom += gam_j*d2_ij;
    }
    /*--- Calculate species laminar viscosity ---*/
    Mu += (Mi/Na * gam_i) / denom;
  }
  if (ionization) {
    iSpecies = nSpecies-1;
    denom = 0.0;
    /*--- Calculate molar concentration ---*/
    Mi    = MolarMass[iSpecies];
    gam_i = rhos[iSpecies] / (Density*Mi);
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      Mj    = MolarMass[jSpecies];
      gam_j = rhos[jSpecies] / (Density*Mj);
      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11(iSpecies,jSpecies,3)
          * pow(Tve, Omega11(iSpecies,jSpecies,0)*log(Tve)*log(Tve)
          + Omega11(iSpecies,jSpecies,1)*log(Tve)
          + Omega11(iSpecies,jSpecies,2));
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
      /*--- Add to denominator of viscosity ---*/
      denom += gam_j*d2_ij;
    }
    Mu += (Mi/Na * gam_i) / denom;
  }
}

void CSU2TCLib::ThermalConductivitiesGY(){

  su2double Cvve, Mi, Mj, mi, mj, pi, R, Na, kb, gam_i, gam_j, denom_t, denom_r, d1_ij, d2_ij, a_ij, Omega_ij, rhoCvve;

  pi   = PI_NUMBER;
  Na   = AVOGAD_CONSTANT;
  kb   = BOLTZMANN_CONSTANT;

  if (ionization) {
    SU2_MPI::Error("NEEDS REVISION w/ IONIZATION",CURRENT_FUNCTION);
  }

  /*--- Mixture vibrational-electronic specific heat ---*/
  Cvves = ComputeSpeciesCvVibEle();
  rhoCvve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoCvve += rhos[iSpecies]*Cvves[iSpecies];
  Cvve = rhoCvve/Density;

  /*--- Calculate mixture gas constant ---*/
  R = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    R += Ru * rhos[iSpecies]/Density;
  }
  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  ThermalCond_tr    = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    /*--- Calculate molar concentration ---*/
    Mi      = MolarMass[iSpecies];
    mi      = Mi/Na;
    gam_i   = rhos[iSpecies] / (Density*Mi);
    denom_t = 0.0;
    denom_r = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      Mj    = MolarMass[jSpecies];
      mj    = Mj/Na;
      gam_j = rhos[iSpecies] / (Density*Mj);
      a_ij = 1.0 + (1.0 - mi/mj)*(0.45 - 2.54*mi/mj) / ((1.0 + mi/mj)*(1.0 + mi/mj));
      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega00(iSpecies,jSpecies,3)
          * pow(T, Omega00(iSpecies,jSpecies,0)*log(T)*log(T)
          + Omega00(iSpecies,jSpecies,1)*log(T)
          + Omega00(iSpecies,jSpecies,2));
      /*--- Calculate "delta1_ij" ---*/
      d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
      /*--- Calculate the Omega^(1,1)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega11(iSpecies,jSpecies,3)
          * pow(T, Omega11(iSpecies,jSpecies,0)*log(T)*log(T)
          + Omega11(iSpecies,jSpecies,1)*log(T)
          + Omega11(iSpecies,jSpecies,2));
      /*--- Calculate "delta2_ij" ---*/
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;
      denom_t += a_ij*gam_j*d2_ij;
      denom_r += gam_j*d1_ij;
    }
    /*--- Translational contribution to thermal conductivity ---*/
    ThermalCond_tr    += (15.0/4.0)*kb*gam_i/denom_t;
    /*--- Translational contribution to thermal conductivity ---*/
    if (RotationModes[iSpecies] != 0.0)
      ThermalCond_tr  += kb*gam_i/denom_r;
    /*--- Vibrational-electronic contribution to thermal conductivity ---*/
    ThermalCond_ve += kb*Cvve/R*gam_i / denom_r;
  }

  ThermalConductivities[0] = ThermalCond_tr;
  ThermalConductivities[1] = ThermalCond_ve;

}

vector<su2double>& CSU2TCLib::ComputeTemperatures(vector<su2double>& val_rhos, su2double rhoE, su2double rhoEve, su2double rhoEvel){

  vector<su2double> val_eves;
  su2double rhoCvtr, rhoE_f, rhoE_ref, rhoEve_t, Tve2, Tve_o, Btol, Tmin, Tmax;
  bool Bconvg;
  unsigned short iIter, maxBIter;

  rhos = val_rhos;

  /*----------Translational temperature----------*/

  rhoE_f   = 0.0;
  rhoE_ref = 0.0;
  rhoCvtr  = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr  += rhos[iSpecies] * Cvtrs[iSpecies];
    rhoE_ref += rhos[iSpecies] * Cvtrs[iSpecies] * Ref_Temperature[iSpecies];
    rhoE_f   += rhos[iSpecies] * (Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies]);
  }

  T = (rhoE - rhoEve - rhoE_f + rhoE_ref - rhoEvel) / rhoCvtr;

  /*--- Set temperature clipping values ---*/
  Tmin  = 50.0; Tmax = 8E4;
  Tve_o = 50.0; Tve2 = 8E4;

  /* Determine if the temperature lies within the acceptable range */
  if      (T < Tmin) T = Tmin;
  else if (T > Tmax) T = Tmax;

  /*--- Set vibrational temperature algorithm parameters ---*/
  Btol     = 1.0E-6;    // Tolerance for the Bisection method
  maxBIter = 50;        // Maximum Bisection method iterations

  //Initialize solution
  Tve   = T;

  // Execute the root-finding method
  Bconvg = false;

  for (iIter = 0; iIter < maxBIter; iIter++) {
    Tve      = (Tve_o+Tve2)/2.0;
    val_eves = ComputeSpeciesEve(Tve);
    rhoEve_t = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) rhoEve_t += rhos[iSpecies] * val_eves[iSpecies];
    if (fabs(rhoEve_t - rhoEve) < Btol) {
      Bconvg = true;
      break;
    } else {
      if (rhoEve_t > rhoEve) Tve2 = Tve;
      else                  Tve_o = Tve;
    }
  }
  // If absolutely no convergence, then assign to the TR temperature
  if (!Bconvg) Tve = T;

  temperatures[0] = T;
  temperatures[1] = Tve;

  return temperatures;

}

void CSU2TCLib::GetChemistryEquilConstants(unsigned short iReaction){

  if (gas_model == "O2"){

    //O2 + M -> 2O + M
    RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
    RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
    RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
    RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
    RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
    RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

  } else if (gas_model == "N2"){

    //N2 + M -> 2N + M
    RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
    RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
    RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
    RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
    RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
    RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

  } else if (gas_model == "ARGON_SID"){

    //N2 + M -> 2N + M
    RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
    RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
    RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
    RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
    RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
    RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

  } else if (gas_model == "AIR-5"){

    if (iReaction <= 4) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 4 && iReaction <= 9) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 9 && iReaction <= 14) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 15) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 16) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;
    }

  } else if (gas_model == "AIR-7"){

    if (iReaction <= 6) {

      //N2 + M -> 2N + M
      RxnConstantTable(0,0) = 3.4907;  RxnConstantTable(0,1) = 0.83133; RxnConstantTable(0,2) = 4.0978;  RxnConstantTable(0,3) = -12.728; RxnConstantTable(0,4) = 0.07487;   //n = 1E14
      RxnConstantTable(1,0) = 2.0723;  RxnConstantTable(1,1) = 1.38970; RxnConstantTable(1,2) = 2.0617;  RxnConstantTable(1,3) = -11.828; RxnConstantTable(1,4) = 0.015105;  //n = 1E15
      RxnConstantTable(2,0) = 1.6060;  RxnConstantTable(2,1) = 1.57320; RxnConstantTable(2,2) = 1.3923;  RxnConstantTable(2,3) = -11.533; RxnConstantTable(2,4) = -0.004543; //n = 1E16
      RxnConstantTable(3,0) = 1.5351;  RxnConstantTable(3,1) = 1.60610; RxnConstantTable(3,2) = 1.2993;  RxnConstantTable(3,3) = -11.494; RxnConstantTable(3,4) = -0.00698;  //n = 1E17
      RxnConstantTable(4,0) = 1.4766;  RxnConstantTable(4,1) = 1.62910; RxnConstantTable(4,2) = 1.2153;  RxnConstantTable(4,3) = -11.457; RxnConstantTable(4,4) = -0.00944;  //n = 1E18
      RxnConstantTable(5,0) = 1.4766;  RxnConstantTable(5,1) = 1.62910; RxnConstantTable(5,2) = 1.2153;  RxnConstantTable(5,3) = -11.457; RxnConstantTable(5,4) = -0.00944;  //n = 1E19

    } else if (iReaction > 6 && iReaction <= 13) {

      //O2 + M -> 2O + M
      RxnConstantTable(0,0) = 1.8103;  RxnConstantTable(0,1) = 1.9607;  RxnConstantTable(0,2) = 3.5716;  RxnConstantTable(0,3) = -7.3623;   RxnConstantTable(0,4) = 0.083861;
      RxnConstantTable(1,0) = 0.91354; RxnConstantTable(1,1) = 2.3160;  RxnConstantTable(1,2) = 2.2885;  RxnConstantTable(1,3) = -6.7969;   RxnConstantTable(1,4) = 0.046338;
      RxnConstantTable(2,0) = 0.64183; RxnConstantTable(2,1) = 2.4253;  RxnConstantTable(2,2) = 1.9026;  RxnConstantTable(2,3) = -6.6277;   RxnConstantTable(2,4) = 0.035151;
      RxnConstantTable(3,0) = 0.55388; RxnConstantTable(3,1) = 2.4600;  RxnConstantTable(3,2) = 1.7763;  RxnConstantTable(3,3) = -6.5720;   RxnConstantTable(3,4) = 0.031445;
      RxnConstantTable(4,0) = 0.52455; RxnConstantTable(4,1) = 2.4715;  RxnConstantTable(4,2) = 1.7342;  RxnConstantTable(4,3) = -6.55534;  RxnConstantTable(4,4) = 0.030209;
      RxnConstantTable(5,0) = 0.50989; RxnConstantTable(5,1) = 2.4773;  RxnConstantTable(5,2) = 1.7132;  RxnConstantTable(5,3) = -6.5441;   RxnConstantTable(5,4) = 0.029591;

    } else if (iReaction > 13 && iReaction <= 20) {

      //NO + M -> N + O + M
      RxnConstantTable(0,0) = 2.1649;  RxnConstantTable(0,1) = 0.078577;  RxnConstantTable(0,2) = 2.8508;  RxnConstantTable(0,3) = -8.5422; RxnConstantTable(0,4) = 0.053043;
      RxnConstantTable(1,0) = 1.0072;  RxnConstantTable(1,1) = 0.53545;   RxnConstantTable(1,2) = 1.1911;  RxnConstantTable(1,3) = -7.8098; RxnConstantTable(1,4) = 0.004394;
      RxnConstantTable(2,0) = 0.63817; RxnConstantTable(2,1) = 0.68189;   RxnConstantTable(2,2) = 0.66336; RxnConstantTable(2,3) = -7.5773; RxnConstantTable(2,4) = -0.011025;
      RxnConstantTable(3,0) = 0.55889; RxnConstantTable(3,1) = 0.71558;   RxnConstantTable(3,2) = 0.55396; RxnConstantTable(3,3) = -7.5304; RxnConstantTable(3,4) = -0.014089;
      RxnConstantTable(4,0) = 0.5150;  RxnConstantTable(4,1) = 0.73286;   RxnConstantTable(4,2) = 0.49096; RxnConstantTable(4,3) = -7.5025; RxnConstantTable(4,4) = -0.015938;
      RxnConstantTable(5,0) = 0.50765; RxnConstantTable(5,1) = 0.73575;   RxnConstantTable(5,2) = 0.48042; RxnConstantTable(5,3) = -7.4979; RxnConstantTable(5,4) = -0.016247;

    } else if (iReaction == 21) {

      //N2 + O -> NO + N
      RxnConstantTable(0,0) = 1.3261;  RxnConstantTable(0,1) = 0.75268; RxnConstantTable(0,2) = 1.2474;  RxnConstantTable(0,3) = -4.1857; RxnConstantTable(0,4) = 0.02184;
      RxnConstantTable(1,0) = 1.0653;  RxnConstantTable(1,1) = 0.85417; RxnConstantTable(1,2) = 0.87093; RxnConstantTable(1,3) = -4.0188; RxnConstantTable(1,4) = 0.010721;
      RxnConstantTable(2,0) = 0.96794; RxnConstantTable(2,1) = 0.89131; RxnConstantTable(2,2) = 0.7291;  RxnConstantTable(2,3) = -3.9555; RxnConstantTable(2,4) = 0.006488;
      RxnConstantTable(3,0) = 0.97646; RxnConstantTable(3,1) = 0.89043; RxnConstantTable(3,2) = 0.74572; RxnConstantTable(3,3) = -3.9642; RxnConstantTable(3,4) = 0.007123;
      RxnConstantTable(4,0) = 0.96188; RxnConstantTable(4,1) = 0.89617; RxnConstantTable(4,2) = 0.72479; RxnConstantTable(4,3) = -3.955;  RxnConstantTable(4,4) = 0.006509;
      RxnConstantTable(5,0) = 0.96921; RxnConstantTable(5,1) = 0.89329; RxnConstantTable(5,2) = 0.73531; RxnConstantTable(5,3) = -3.9596; RxnConstantTable(5,4) = 0.006818;

    } else if (iReaction == 22) {

      //NO + O -> O2 + N
      RxnConstantTable(0,0) = 0.35438;   RxnConstantTable(0,1) = -1.8821; RxnConstantTable(0,2) = -0.72111;  RxnConstantTable(0,3) = -1.1797;   RxnConstantTable(0,4) = -0.030831;
      RxnConstantTable(1,0) = 0.093613;  RxnConstantTable(1,1) = -1.7806; RxnConstantTable(1,2) = -1.0975;   RxnConstantTable(1,3) = -1.0128;   RxnConstantTable(1,4) = -0.041949;
      RxnConstantTable(2,0) = -0.003732; RxnConstantTable(2,1) = -1.7434; RxnConstantTable(2,2) = -1.2394;   RxnConstantTable(2,3) = -0.94952;  RxnConstantTable(2,4) = -0.046182;
      RxnConstantTable(3,0) = 0.004815;  RxnConstantTable(3,1) = -1.7443; RxnConstantTable(3,2) = -1.2227;   RxnConstantTable(3,3) = -0.95824;  RxnConstantTable(3,4) = -0.045545;
      RxnConstantTable(4,0) = -0.009758; RxnConstantTable(4,1) = -1.7386; RxnConstantTable(4,2) = -1.2436;   RxnConstantTable(4,3) = -0.949;    RxnConstantTable(4,4) = -0.046159;
      RxnConstantTable(5,0) = -0.002428; RxnConstantTable(5,1) = -1.7415; RxnConstantTable(5,2) = -1.2331;   RxnConstantTable(5,3) = -0.95365;  RxnConstantTable(5,4) = -0.04585;

    } else if (iReaction == 23) {

      //N + O -> NO+ + e-
      RxnConstantTable(0,0) = -2.1852;   RxnConstantTable(0,1) = -6.6709; RxnConstantTable(0,2) = -4.2968; RxnConstantTable(0,3) = -2.2175; RxnConstantTable(0,4) = -0.050748;
      RxnConstantTable(1,0) = -1.0276;   RxnConstantTable(1,1) = -7.1278; RxnConstantTable(1,2) = -2.637;  RxnConstantTable(1,3) = -2.95;   RxnConstantTable(1,4) = -0.0021;
      RxnConstantTable(2,0) = -0.65871;  RxnConstantTable(2,1) = -7.2742; RxnConstantTable(2,2) = -2.1096; RxnConstantTable(2,3) = -3.1823; RxnConstantTable(2,4) = 0.01331;
      RxnConstantTable(3,0) = -0.57924;  RxnConstantTable(3,1) = -7.3079; RxnConstantTable(3,2) = -1.9999; RxnConstantTable(3,3) = -3.2294; RxnConstantTable(3,4) = 0.016382;
      RxnConstantTable(4,0) = -0.53538;  RxnConstantTable(4,1) = -7.3252; RxnConstantTable(4,2) = -1.937;  RxnConstantTable(4,3) = -3.2572; RxnConstantTable(4,4) = 0.01823;
      RxnConstantTable(5,0) = -0.52801;  RxnConstantTable(5,1) = -7.3281; RxnConstantTable(5,2) = -1.9264; RxnConstantTable(5,3) = -3.2618; RxnConstantTable(5,4) = 0.01854;
    }
  }
}
