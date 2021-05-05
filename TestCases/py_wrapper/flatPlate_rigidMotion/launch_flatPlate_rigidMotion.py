#!/usr/bin/env python

## \file flatPlate_rigidMotion.py
#  \brief Python script to launch SU2_CFD with customized unsteady boundary conditions using the Python wrapper.
#  \author David Thomas
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

import sys
from optparse import OptionParser	# use a parser for configuration
import pysu2			            # imports the SU2 wrapped module
from math import *
import numpy as np

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

  # Command line options
  parser=OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
  parser.add_option("--parallel", action="store_true",
                    help="Specify if we need to initialize MPI", dest="with_MPI", default=False)

  (options, args) = parser.parse_args()
  options.nDim  = int(2)
  options.nZone = int(1)

  # Import mpi4py for parallel run
  if options.with_MPI == True:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
  else:
    comm = 0
    rank = 0

  # Initialize the corresponding driver of SU2, this includes solver preprocessing
  try:
      SU2Driver = pysu2.CSinglezoneDriver(options.filename, options.nZone, comm);
  except TypeError as exception:
    print('A TypeError occured in pysu2.CDriver : ',exception)
    if options.with_MPI == True:
      print('ERROR : You are trying to initialize MPI with a serial build of the wrapper. Please, remove the --parallel option that is incompatible with a serial build.')
    else:
      print('ERROR : You are trying to launch a computation without initializing MPI but the wrapper has been built in parallel. Please add the --parallel option in order to initialize MPI for the wrapper.')
    return


  MovingMarkerID = None
  MovingMarker = 'plate'       #specified by the user

  # Get all the tags with the moving option
  MovingMarkerList =  SU2Driver.GetAllDeformMeshMarkersTag()

  # Get all the markers defined on this rank and their associated indices.
  allMarkerIDs = SU2Driver.GetAllBoundaryMarkers()

  # Check if the specified marker has a moving option and if it exists on this rank.
  if MovingMarker in MovingMarkerList and MovingMarker in allMarkerIDs.keys():
    MovingMarkerID = allMarkerIDs[MovingMarker]

  # Number of vertices on the specified marker (per rank)
  nVertex_MovingMarker = 0         #total number of vertices (physical + halo)
  nVertex_MovingMarker_HALO = 0    #number of halo vertices
  nVertex_MovingMarker_PHYS = 0    #number of physical vertices

  if MovingMarkerID != None:
    nVertex_MovingMarker = SU2Driver.GetNumberVertices(MovingMarkerID)
    nVertex_MovingMarker_HALO = SU2Driver.GetNumberHaloVertices(MovingMarkerID)
    nVertex_MovingMarker_PHYS = nVertex_MovingMarker - nVertex_MovingMarker_HALO

  # Retrieve some control parameters from the driver
  deltaT = SU2Driver.GetUnsteady_TimeStep()
  TimeIter = SU2Driver.GetTime_Iter()
  nTimeIter = SU2Driver.GetnTimeIter()
  time = TimeIter*deltaT

  # Extract the initial position of each node on the moving marker
  CoordX = np.zeros(nVertex_MovingMarker)
  CoordY = np.zeros(nVertex_MovingMarker)
  CoordZ = np.zeros(nVertex_MovingMarker)
  for iVertex in range(nVertex_MovingMarker):
    CoordX[iVertex], CoordY[iVertex], CoordZ[iVertex] = SU2Driver.GetInitialMeshCoord(MovingMarkerID, iVertex)

  # Time loop is defined in Python so that we have acces to SU2 functionalities at each time step
  if rank == 0:
    print("\n------------------------------ Begin Solver -----------------------------\n")
  sys.stdout.flush()
  if options.with_MPI == True:
    comm.Barrier()

  while (TimeIter < nTimeIter):
    # Define the rigid body displacement and set the new coords of each node on the marker
    d_y = 0.0175*sin(2*pi*time)
    for iVertex in range(nVertex_MovingMarker):
      SU2Driver.SetMeshDisplacement(MovingMarkerID, int(iVertex), 0.0, d_y, 0.0)
    # Time iteration preprocessing
    SU2Driver.Preprocess(TimeIter)
    # Run one time iteration (e.g. dual-time)
    SU2Driver.Run()
    # Postprocess the solver
    SU2Driver.Postprocess()
    # Update the solver for the next time iteration
    SU2Driver.Update()
    # Monitor the solver and output solution to file if required
    stopCalc = SU2Driver.Monitor(TimeIter)
    SU2Driver.Output(TimeIter)
    if (stopCalc == True):
      break
    # Update control parameters
    TimeIter += 1
    time += deltaT

  # Postprocess the solver and exit cleanly
  SU2Driver.Postprocessing()

  if SU2Driver != None:
    del SU2Driver

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
