/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/NeighborhoodOffsets.h"


namespace mirtk {


// -----------------------------------------------------------------------------
NeighborhoodOffsets::NeighborhoodOffsets()
:
  _Connectivity(CONNECTIVITY_26),
  _Size(0)
{
  for (int i = 0; i < 26; ++i) _Offsets[i] = 0;
}

// -----------------------------------------------------------------------------
NeighborhoodOffsets::NeighborhoodOffsets(const BaseImage* image, ConnectivityType connectivity)
{
  Initialize(image, connectivity);
}

// -----------------------------------------------------------------------------
void NeighborhoodOffsets::Initialize(const BaseImage* image, ConnectivityType connectivity)
{
  Initialize(image->X(), image->Y(), connectivity);
}

// -----------------------------------------------------------------------------
void NeighborhoodOffsets::Initialize(int xdim, int ydim, ConnectivityType connectivity)
{
  switch (connectivity){
    case CONNECTIVITY_4:
      _Size = 4;
      break;
    case CONNECTIVITY_6:
      _Size = 6;
      break;
    case CONNECTIVITY_18:
      _Size = 18;
      break;
    case CONNECTIVITY_26:
      _Size = 26;
      break;
    default:
      cerr << "NeighborhoodOffsets: Invalid connectivity type" << endl;
      exit(1);
  }

  _Connectivity = connectivity;

  // Bare minimum, 4 connectivity.
  // Get the face neighbours along X and Y.
  _Offsets[0] =  1;
  _Offsets[1] = -1;
  _Offsets[2] = xdim * -1;
  _Offsets[3] = xdim;

  if (_Connectivity == CONNECTIVITY_6 ||
      _Connectivity == CONNECTIVITY_18 ||
      _Connectivity == CONNECTIVITY_26){
    // Add the face neighbours along Z.
    _Offsets[4] = xdim * ydim * -1;
    _Offsets[5] = xdim * ydim;
  }

  if (_Connectivity == CONNECTIVITY_18 ||
      _Connectivity == CONNECTIVITY_26){
    // Add the edge neighbours.
    _Offsets[6]  =  1 + _Offsets[2];
    _Offsets[7]  =  1 + _Offsets[3];
    _Offsets[8]  =  1 + _Offsets[4];
    _Offsets[9]  =  1 + _Offsets[5];

    _Offsets[10] = -1 + _Offsets[2];
    _Offsets[11] = -1 + _Offsets[3];
    _Offsets[12] = -1 + _Offsets[4];
    _Offsets[13] = -1 + _Offsets[5];

    _Offsets[14] = _Offsets[2] + _Offsets[4];
    _Offsets[15] = _Offsets[2] + _Offsets[5];
    _Offsets[16] = _Offsets[3] + _Offsets[4];
    _Offsets[17] = _Offsets[3] + _Offsets[5];
  }

  if (_Connectivity == CONNECTIVITY_26){
    // Add the vertex neighbours for the 26 neighbourhood.
    _Offsets[18] =  1 + _Offsets[14];
    _Offsets[19] =  1 + _Offsets[15];
    _Offsets[20] =  1 + _Offsets[16];
    _Offsets[21] =  1 + _Offsets[17];

    _Offsets[22] = -1 + _Offsets[14];
    _Offsets[23] = -1 + _Offsets[15];
    _Offsets[24] = -1 + _Offsets[16];
    _Offsets[25] = -1 + _Offsets[17];
  }
}

// -----------------------------------------------------------------------------
NeighborhoodOffsets::~NeighborhoodOffsets()
{
}


} // namespace mirtk
