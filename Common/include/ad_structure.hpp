/*!
 * \file ad_structure.hpp
 * \brief Main routines for the algorithmic differentiation (AD) structure.
 * \author T. Albring
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

#pragma once

#include "../include/datatype_structure.hpp"

/*!
 * \namespace AD
 * \brief Contains routines for the reverse mode of AD.
 * In case there is no reverse type configured, they have no effect at all.
 */

namespace AD{
  /*!
   * \brief Start the recording of the operations and involved variables.
   * If called, the computational graph of all operations occuring after the call will be stored,
   * starting with the variables registered with RegisterInput.
   */
  void StartRecording();

  /*!
   * \brief Stops the recording of the operations and variables.
   */
  void StopRecording();

  /*!
   * \brief Check if the tape is active
   * \param[out] Boolean which determines whether the tape is active.
   */
  bool TapeActive();

  /*!
   * \brief Prints out tape statistics.
   */
  void PrintStatistics();

  /*!
   * \brief Registers the variable as an input and saves internal data (indices). I.e. as a leaf of the computational graph.
   * \param[in] data - The variable to be registered as input.
   * \param[in] push_index - boolean whether we also want to push the index.
   */
  void RegisterInput(su2double &data, bool push_index = true);

  /*!
   * \brief Registers the variable as an output. I.e. as the root of the computational graph.
   * \param[in] data - The variable to be registered as output.
   */
  void RegisterOutput(su2double &data);

  /*!
   * \brief Sets the adjoint value at index to val
   * \param[in] index - Position in the adjoint vector.
   * \param[in] val - adjoint value to be set.
   */
  void SetDerivative(int index, const double val);

  /*!
   * \brief Extracts the adjoint value at index
   * \param[in] index - position in the adjoint vector where the derivative will be extracted.
   * \return Derivative value.
   */
  double GetDerivative(int index);

  /*!
   * \brief Clears the currently stored adjoints but keeps the computational graph.
   */
  void ClearAdjoints();

  /*!
   * \brief Computes the adjoints, i.e. the derivatives of the output with respect to the input variables.
   */
  void ComputeAdjoint();  
  
  /*!
   * \brief Computes the adjoints, i.e. the derivatives of the output with respect to the input variables.
   * \param[in] enter - Position where we start evaluating the tape.
   * \param[in] leave - Position where we stop evaluating the tape.
   */
  void ComputeAdjoint(unsigned short enter, unsigned short leave);

  /*!
   * \brief Reset the tape structure to be ready for a new recording.
   */
  void Reset();

  /*!
   * \brief Reset the variable (set index to zero).
   * \param[in] data - the variable to be unregistered from the tape.
   */
  void ResetInput(su2double &data);

  /*!
   * \brief Sets the scalar input of a preaccumulation section.
   * \param[in] data - the scalar input variable.
   */
  void SetPreaccIn(const su2double &data);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 1D array.
   * \param[in] data - the input 1D array.
   * \param[in] size - size of the array.
   */
  void SetPreaccIn(const su2double* data, const int size);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the input 2D array.
   * \param[in] size_x - size of the array in x dimension.
   * \param[in] size_y - size of the array in y dimension.
   */
  void SetPreaccIn(const su2double* const *data, const int size_x, const int size_y);

  /*!
   * \brief Starts a new preaccumulation section and sets the input variables.
   *
   * The idea of preaccumulation is to store only the Jacobi matrix of a code section during
   * the taping process instead of all operations. This decreases the tape size and reduces runtime.
   *
   * Input/Output of the section are set with several calls to SetPreaccIn()/SetPreaccOut().
   *
   * Note: the call of this routine must be followed by a call of EndPreacc() and the end of the code section.
   */
  void StartPreacc();

  /*!
   * \brief Sets the scalar output of a preaccumulation section.
   * \param[in] data - the scalar output variable.
   */
  void SetPreaccOut(su2double &data);

  /*!
   * \brief Sets the output variables of a preaccumulation section using a 1D array.
   * \param[in] data - the output 1D array.
   */
  void SetPreaccOut(su2double* data, const int size);

  /*!
   * \brief Sets the input variables of a preaccumulation section using a 2D array.
   * \param[in] data - the output 1D array.
   */
  void SetPreaccOut(su2double** data, const int size_x, const int size_y);

  /*!
   * \brief Ends a preaccumulation section and computes the local Jacobi matrix
   * of a code section using the variables set with SetLocalInput(), SetLocalOutput() and pushes a statement
   * for each output variable to the AD tape.
   */
  void EndPreacc();
  
  /*!
   * \brief Initializes an externally differentiated function. Input and output variables are set with SetExtFuncIn/SetExtFuncOut
   * \param[in] storePrimalInput - Specifies whether the primal input values are stored for the reverse call of the external function.
   * \param[in] storePrimalOutput - Specifies whether the primal output values are stored for the reverse call of the external function.
   */
  void StartExtFunc(bool storePrimalInput, bool storePrimalOutput);
  
  /*!
   * \brief Sets the scalar input of a externally differentiated function.
   * \param[in] data - the scalar input variable.
   */
  void SetExtFuncIn(su2double &data);
  
  /*!
   * \brief Sets the input variables of a externally differentiated function using a 1D array.
   * \param[in] data - the input 1D array.
   * \param[in] size - number of rows.
   */
  void SetExtFuncIn(const su2double* data, const int size);
  
  /*!
  * \brief  Sets the input variables of a externally differentiated function using a 2D array.
  * \param[in] data - the input 2D array.
  * \param[in] size_x - number of rows.
  * \param[in] size_y - number of columns.
  */
  void SetExtFuncIn(const su2double* const *data, const int size_x, const int size_y);

  /*!
   * \brief Sets the scalar output of a externally differentiated function.
   * \param[in] data - the scalar output variable.
   */
  void SetExtFuncOut(su2double &data);

  /*!
   * \brief Sets the output variables of a externally differentiated function using a 1D array.
   * \param[in] data - the output 1D array.
   * \param[in] size - number of rows.
   */
  void SetExtFuncOut(su2double* data, const int size);

  /*!
  * \brief  Sets the output variables of a externally differentiated function using a 2D array.
  * \param[in] data - the output 2D array.
  * \param[in] size_x - number of rows.
  * \param[in] size_y - number of columns.
  */
  void SetExtFuncOut(su2double** data, const int size_x, const int size_y);
  
  /*!
   * \brief Ends an external function section by deleting the structures.
   */
  void EndExtFunc();

  /*!
   * \brief Evaluates and saves gradient data from a variable.
   * \param[in] data - variable whose gradient information will be extracted.
   * \param[in] index - where obtained gradient information will be stored.
   */
  void SetIndex(int &index, const su2double &data);

  /*!
   * \brief Pushes back the current tape position to the tape position's vector.
   */
  void Push_TapePosition();
}

/*--- Macro to begin and end sections with a passive tape ---*/

#ifdef CODI_REVERSE_TYPE
#define AD_BEGIN_PASSIVE         \
  if(AD::globalTape.isActive()) {\
     AD::globalTape.setPassive();\
     AD::Status = true;          \
  }
#define AD_END_PASSIVE           \
  if(AD::Status) {               \
     AD::globalTape.setActive(); \
     AD::Status = false;         \
  }
#else
#define AD_BEGIN_PASSIVE
#define AD_END_PASSIVE
#endif

/*--- If we compile under OSX we have to overload some of the operators for
 *   complex numbers to avoid the use of the standard operators
 *  (they use a lot of functions that are only defined for doubles) ---*/

#ifdef __APPLE__

namespace std{
  template<>
  su2double abs(const complex<su2double>& x);

  template<>
  complex<su2double> operator/(const complex<su2double>& x, const complex<su2double>& y);

  template<>
  complex<su2double> operator*(const complex<su2double>& x, const complex<su2double>& y);
}
#endif
#include "ad_structure.inl"
