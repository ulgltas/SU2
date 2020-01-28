/*!
 * \file CFileWriter.cpp
 * \brief Filewriter base class.
 * \author T. Albring
 * \version 7.0.0 "Blackbird"
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

#include "../../../include/output/filewriter/CFileWriter.hpp"


CFileWriter::CFileWriter(vector<string> fields, string fileName,
                         CParallelDataSorter *dataSorter, string file_ext, unsigned short nDim):
  fieldnames(std::move(fields)),
  nDim(nDim),
  file_ext(file_ext),
  fileName(std::move(fileName)),
  dataSorter(dataSorter){

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  this->fileName += file_ext;

  file_size = 0.0;

}


CFileWriter::~CFileWriter(){

}
