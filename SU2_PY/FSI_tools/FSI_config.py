#!/usr/bin/env python

## \file FSI_config.py
#  \brief Python class for handling configuration file for FSI computation.
#  \authors Nicola Fonzi, Vittorio Cavalieri based on the work of David Thomas
#  \version 7.1.1 "Blackbird"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2020, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os
import sys
import shutil
import copy

# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class FSIConfig:
    """
    Class that contains all the parameters coming from the FSI configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self,FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.readConfig()
        self.applyDefaults()

    def __str__(self):
        tempString = str()
        for key, value in self._ConfigContent.items():
            tempString += "{} = {}\n".format(key,value)
        return tempString

    def __getitem__(self,key):
        return self._ConfigContent[key]

    def __setitem__(self, key, value):
        self._ConfigContent[key] = value

    def readConfig(self):
        input_file = open(self.ConfigFileName)
        while 1:
            line = input_file.readline()
            if not line:
                break
            # remove line returns
            line = line.strip('\r\n')
            # make sure it has useful data
            if (not "=" in line) or (line[0] == '%'):
                continue
            # split across equal sign
            line = line.split("=",1)
            this_param = line[0].strip()
            this_value = line[1].strip()

	        #integer values
            if (this_param == "NDIM")	or \
               (this_param == "RESTART_ITER")	or \
               (this_param == "TIME_TRESHOLD") or \
               (this_param == "NB_FSI_ITER") :
                self._ConfigContent[this_param] = int(this_value)

            #float values
            elif (this_param == "RBF_RADIUS") or \
                 (this_param == "AITKEN_PARAM") or \
                 (this_param == "UNST_TIMESTEP") or \
                 (this_param == "UNST_TIME") or \
                 (this_param == "FSI_TOLERANCE") :
                self._ConfigContent[this_param] = float(this_value)

            #string values
            elif (this_param == "CFD_CONFIG_FILE_NAME") or \
                 (this_param == "CSD_SOLVER") or \
                 (this_param == "CSD_CONFIG_FILE_NAME") or \
                 (this_param == "RESTART_SOL") or \
                 (this_param == "MATCHING_MESH") or \
                 (this_param == "MESH_INTERP_METHOD") or \
                 (this_param == "DISP_PRED") or \
                 (this_param == "AITKEN_RELAX") or \
                 (this_param == "TIME_MARCHING") :
                self._ConfigContent[this_param] = this_value

            else :
                print(this_param + " is an invalid option !")

    def applyDefaults(self):
        if self._ConfigContent["CSD_SOLVER"] == "IMPOSED":
            self._ConfigContent["AITKEN_RELAX"] = "STATIC"
            self._ConfigContent["AITKEN_PARAM"] = 1.0
