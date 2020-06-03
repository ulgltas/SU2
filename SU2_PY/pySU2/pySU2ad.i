/*
################################################################################
#
# \file pySU2ad.i
# \brief Configuration file for the Swig compilation of the Python wrapper.
# \author D. Thomas, R. Sanchez
#  \version 7.0.0 "Falcon"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
*/

%feature("autodoc","1");

%module(docstring=
"'pysu2ad' module",
directors="1",
threads="1"
) pysu2ad
%{

#include "../../SU2_CFD/include/drivers/CDriver.hpp"
#include "../../SU2_CFD/include/drivers/CSinglezoneDriver.hpp"
#include "../../SU2_CFD/include/drivers/CMultizoneDriver.hpp"
#include "../../SU2_CFD/include/drivers/CDiscAdjSinglezoneDriver.hpp"

%}

// ----------- USED MODULES ------------
%import "../../Common/include/basic_types/datatype_structure.hpp"
%import "../../Common/include/mpi_structure.hpp"
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "typemaps.i"
//%include "numpy.i"
#ifdef HAVE_MPI                    //Need mpi4py only for a parallel build of the wrapper.
  %include "mpi4py/mpi4py.i"
  %mpi4py_typemap(Comm, MPI_Comm)
#endif

namespace std {
   %template() vector<int>;
   %template() vector<double>;
   %template() vector<string>;
   %template() map<string, int>;
   %template() map<string, string>;
}

// ----------- API CLASSES ----------------

//Constants definitions
/*!
 * \brief different software components of SU2
 */
enum SU2_COMPONENT {
  SU2_CFD = 1,	/*!< \brief Running the SU2_CFD software. */
  SU2_DEF = 2,	/*!< \brief Running the SU2_DEF software. */
  SU2_DOT = 3,	/*!< \brief Running the SU2_DOT software. */
  SU2_MSH = 4,	/*!< \brief Running the SU2_MSH software. */
  SU2_GEO = 5,	/*!< \brief Running the SU2_GEO software. */
  SU2_SOL = 6 	/*!< \brief Running the SU2_SOL software. */
};

const unsigned int MESH_0 = 0; /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1; /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0; /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1; /*!< \brief Definition of the first grid domain. */

// CDriver class
%include "../../SU2_CFD/include/drivers/CDriver.hpp"
%include "../../SU2_CFD/include/drivers/CSinglezoneDriver.hpp"
%include "../../SU2_CFD/include/drivers/CMultizoneDriver.hpp"
%include "../../SU2_CFD/include/drivers/CDiscAdjSinglezoneDriver.hpp"
