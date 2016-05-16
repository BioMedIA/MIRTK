#!/usr/bin/env python

# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

"""
Usage: convert-mris <input> <output>

Converts a FreeSurfer surface file to VTK using mris_convert.
When the output file name has extension ".vtp", the output file
of mris_convert is converted from legacy VTK file to a XML file.
The positions of the converted surface vertices are further corrected
for any transformation applied by FreeSurfer (i.e., RAS translation).

Requires:
- VTK IO Python wrappers
- FreeSurfer tools: mris_convert, mris_info

Arguments:
  input    Input FreeSurfer surface file (found in subject "surf" directory).
  output   Output VTK surface file. Must have extension .vtk or .vtp.

"""
from __future__ import absolute_import, print_function, unicode_literals

from subprocess import check_call, check_output, STDOUT
from sys import argv, exit
from os.path import basename, abspath
from os import remove
from vtk import vtkPolyDataReader, vtkPolyDataWriter, vtkXMLPolyDataWriter
import re


# ----------------------------------------------------------------------------------------------------------------------
def convert_mris(input_name, output_name):
  """
  Converts a FreeSurfer surface to VTK file format
  """
  # convert surface to VTK format
  if output_name.endswith('.vtp'): temp_name = output_name[0:-1] + 'k'
  else:                            temp_name = output_name
  if not temp_name.endswith('.vtk'):
    raise RuntimeError('Output file name extension must be either .vtk or .vtp')
  check_call(['mris_convert', input_name, temp_name])
  # get surface RAS translation
  out = check_output(["mris_info", input_name], stderr=STDOUT)
  m = re.search("c_\(ras\)\s:\s\((-?\d+\.\d+),\s(-?\d+\.\d+),\s(-?\d+\.\d+)\)", out)
  if m is None: raise RuntimeError('Could not find c_(ras) coordinates in mris_info output!')
  tx = float(m.group(1))
  ty = float(m.group(2))
  tz = float(m.group(3))
  # transform vertex positions to scanner RAS of orig.mgz
  reader = vtkPolyDataReader()
  reader.SetFileName(temp_name)
  reader.Update()
  surface = reader.GetOutput()
  points = surface.GetPoints()
  for i in range(points.GetNumberOfPoints()):
    x, y, z = points.GetPoint(i)
    points.SetPoint(i, x + tx, y + ty, z + tz)
  surface.SetPoints(points)
  if output_name.endswith('.vtp'): writer = vtkXMLPolyDataWriter()
  else:                            writer = vtkPolyDataWriter()
  writer.SetFileName(output_name)
  try:
    writer.SetInputData(surface)
  except:
    writer.SetInput(surface)
  writer.Write()
  if temp_name != output_name:
    remove(temp_name)

# ----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
  if len(argv) != 3:
    print("usage: " + basename(argv[0]) + " <surf> <output>")
    exit(1)
  convert_mris(abspath(argv[1]), abspath(argv[2]))
