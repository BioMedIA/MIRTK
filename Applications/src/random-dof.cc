/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2019 Imperial College London
 * Copyright 2019 Andreas Schuh
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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"
#include "mirtk/Memory.h"
#include "mirtk/String.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/SimilarityTransformation.h"
#include "mirtk/AffineTransformation.h"
#include "mirtk/BSplineFreeFormTransformation3D.h"
#include "mirtk/BSplineFreeFormTransformation4D.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"
#include "mirtk/LinearInterpolateImageFunction.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <dofout> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  This tool creates a random transformation of the requested type.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -seed <uint>\n";
  cout << "    Seed for pseudo-random number generation engine.\n";
  cout << "    Use default random seed if argument is zero. (default: 0)\n";
  cout << "  -type Rigid|Similarity|Affine|FFD|MFFD|SVFFD|MSVFFD\n";
  cout << "    Type of transformation to create. (default: MFFD)\n";
  cout << "  -target <image>\n";
  cout << "    Target image. Required when output is a FFD. (default: none)\n";
  cout << "  -ndim 2|3|4\n";
  cout << "    Number of transformation domain dimensions. (default: target)\n";
  cout << "\n";
  cout << "Homogeneous transformation options:\n";
  cout << "  -translation <float>\n";
  cout << "    Maximum translation in mm. (default: 0)\n";
  cout << "  -rotation <float>\n";
  cout << "    Maximum rotation angle in degrees. (default: 0)\n";
  cout << "  -scaling <float>\n";
  cout << "    Maximum scaling factor in percentage. (default: 0)\n";
  cout << "  -shearing <float>\n";
  cout << "    Maximum shearing angle in degrees. (default: 0)\n";
  cout << "\n";
  cout << "Free-form deformation options:\n";
  cout << "  -displacement, -velocity <float>\n";
  cout << "    Maximum control point coefficient in mm. If a negative value is given,\n";
  cout << "    its absolute value is multiplied by the control point spacing. (default: 0)\n";
  cout << "  -spacing, -ds <float>\n";
  cout << "    Spatial control point spacing of FFD. (default: target voxel size)\n";
  cout << "  -dx <float>\n";
  cout << "    Control point spacing of FFD in x. (default: target voxel size in x)\n";
  cout << "  -dy <float>\n";
  cout << "    Control point spacing of FFD in y. (default: target voxel size in y)\n";
  cout << "  -dz <float>\n";
  cout << "    Control point spacing of FFD in z. (default: target voxel size in z)\n";
  cout << "  -dt <float>\n";
  cout << "    Control point spacing of FFD in t. (default: temporal target spacing)\n";
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Instantiate transformation of specified type
UniquePtr<Transformation>
NewTransformation(const char *type,
                  const ImageAttributes &attr, int ndim,
                  double sx, double sy, double sz, double st)
{
  const string ltype = ToLower(type);
  if (ltype == "rigid") {
    return NewUnique<RigidTransformation>();
  }
  if (ltype == "similarity") {
    return NewUnique<SimilarityTransformation>();
  }
  if (ltype == "affine") {
    return NewUnique<AffineTransformation>();
  }
  if (!attr) {
    return nullptr;
  }
  UniquePtr<Transformation> dof;
  if (ltype == "ffd" || ltype == "bsplineffd" || ltype == "mffd") {
    UniquePtr<FreeFormTransformation> affd;
    if (ndim < 4) {
      affd.reset(new BSplineFreeFormTransformation3D(attr, sx, sy, sz));
    } else {
      affd.reset(new BSplineFreeFormTransformation4D(attr, sx, sy, sz, st));
    }
    if (ltype == "mffd") {
      auto mffd = NewUnique<MultiLevelFreeFormTransformation>();
      mffd->PushLocalTransformation(affd.release());
      dof.reset(mffd.release());
    } else {
      dof.reset(affd.release());
    }
  }
  if (ltype == "svffd" || ltype == "bsplinesvffd" || ltype == "msvffd") {
    auto affd = NewUnique<BSplineFreeFormTransformationSV>(attr, sx, sy, sz);
    if (ltype == "msvffd") {
      auto mffd = NewUnique<MultiLevelFreeFormTransformation>();
      mffd->PushLocalTransformation(affd.release());
      dof.reset(mffd.release());
    } else {
      dof.reset(affd.release());
    }
  }
  return dof;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InitializeIOLibrary();
 
  EXPECTS_POSARGS(1);
  const char *dofout_name = POSARG(1);

  // Parse command-line arguments
  std::mt19937::result_type seed = 0;
  const char *dofout_type = "MFFD";
  const char *target_name = nullptr;
  const char *mask_name = nullptr;

  // Default settings for affine transformation
  double max_translation = .0;
  double max_rotation = .0;
  double max_scaling = .0;
  double max_shearing = .0;

  // Default settings for free-form deformation
  int ndim = 0;
  double sx = -1; // Control point spacing in x = target voxel size by default
  double sy = -1; // Control point spacing in y = target voxel size by default
  double sz = -1; // Control point spacing in z = target voxel size by default
  double st = -1; // Control point spacing in t = target voxel size by default
  double max_displacement = .0;

  for (ALL_OPTIONS) {
    // Output options
    if (OPTION("-type")) dofout_type = ARGUMENT;
    else if (OPTION("-seed")) PARSE_ARGUMENT(seed);
    else if (OPTION("-ndim")) PARSE_ARGUMENT(ndim);
    // FFD lattice attributes
    else if (OPTION("-dx")) PARSE_ARGUMENT(sx);
    else if (OPTION("-dy")) PARSE_ARGUMENT(sy);
    else if (OPTION("-dz")) PARSE_ARGUMENT(sz);
    else if (OPTION("-dt")) PARSE_ARGUMENT(st);
    else if (OPTION("-ds") || OPTION("-spacing")) {
      PARSE_ARGUMENT(sx);
      sy = sz = sx;
    }
    // Target/source image
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-mask")) mask_name = ARGUMENT;
    // Sampling range
    else if (OPTION("-translation")) PARSE_ARGUMENT(max_translation);
    else if (OPTION("-rotation")) PARSE_ARGUMENT(max_rotation);
    else if (OPTION("-scaling")) PARSE_ARGUMENT(max_scaling);
    else if (OPTION("-shearing")) PARSE_ARGUMENT(max_shearing);
    else if (OPTION("-displacement") || OPTION("-velocity")) {
      PARSE_ARGUMENT(max_displacement);
    }
    // Common or unknown option
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read attributes of target image
  RealImage target;
  if (target_name) target.Read(target_name);
  const ImageAttributes &attr = target.Attributes();
  if (ndim < 2) ndim = (attr._z == 1 ? 2 : 3);
  if (ndim > 4) FatalError("Only 2D, 3D, or 4D transformation domains supported");

  // Instantiate transformation
  auto dof = NewTransformation(dofout_type, attr, ndim, sx, sy, sz, st);
  if (!dof) FatalError("Invalid -type, -dim, or missing -target image");

  auto global = dynamic_cast<HomogeneousTransformation *>(dof.get());
  auto affd = dynamic_cast<FreeFormTransformation *>(dof.get());
  auto mffd = dynamic_cast<MultiLevelTransformation *>(dof.get());

  if (mffd) {
    global = mffd->GetGlobalTransformation();
    affd = mffd->GetLocalTransformation(0);
  }

  // Random number generator
  std::random_device rd;
  std::mt19937 mt(rd());
  if (seed != 0) mt.seed(seed);
  std::uniform_real_distribution<double> sampler(-1., 1.);

  if (global) {
    double params[12] = {.0};
    // Translation
    if (!IsZero(max_translation)) {
      params[TX] = sampler(mt) * max_translation;
      params[TY] = sampler(mt) * max_translation;
      if (ndim > 2) params[TZ] = sampler(mt) * max_translation;
    }
    // Rotation
    if (!IsZero(max_rotation)) {
      if (ndim > 2) params[RX] = sampler(mt) * max_rotation;
      if (ndim > 2) params[RY] = sampler(mt) * max_rotation;
      params[RZ] = sampler(mt) * max_rotation;
    }
    // Scaling
    if (IsZero(max_scaling)) {
      params[SX] = 100.;
      params[SY] = 100.;
      params[SZ] = 100.;
    } else {
      params[SX] = 100. + sampler(mt) * max_scaling;
      params[SY] = 100. + sampler(mt) * max_scaling;
      params[SZ] = 100. + (ndim > 2 ? sampler(mt) * max_scaling : .0);
    }
    // Shearing
    if (!IsZero(max_shearing)) {
      params[SXY] = sampler(mt) * max_shearing;
      if (ndim > 2) params[SYZ] = sampler(mt) * max_shearing;
      if (ndim > 2) params[SXZ] = sampler(mt) * max_shearing;
    }
    // Set parameters
    global->Put(params);
  }

  if (affd && !IsZero(max_displacement)) {

    const auto max_d = max_displacement;
    const auto max_displacement_x = (max_d > .0 ? max_d : -max_d * affd->GetXSpacing());
    const auto max_displacement_y = (max_d > .0 ? max_d : -max_d * affd->GetYSpacing());
    const auto max_displacement_z = (max_d > .0 ? max_d : -max_d * affd->GetZSpacing());

    int cp, xdof, ydof, zdof;
    const auto ndofs = affd->NumberOfDOFs();
    Array<double> params(ndofs, .0);
    if (mask_name) {
      double x, y, z, w;
      RealImage mask(mask_name);
      GenericLinearInterpolateImageFunction<RealImage> f_mask;
      f_mask.Input(&mask);
      f_mask.Initialize();
      for (int l = 0; l < affd->T(); ++l)
      for (int k = 0; k < affd->Z(); ++k)
      for (int j = 0; j < affd->Y(); ++j)
      for (int i = 0; i < affd->X(); ++i) {
        x = i, y = j, z = k;
        affd->LatticeToWorld(x, y, z);
        mask.WorldToImage(x, y, z);
        w = f_mask.Evaluate(x, y, z);
        if (!IsZero(w)) {
          cp = affd->LatticeToIndex(i, j, k, l);
          affd->IndexToDOFs(cp, xdof, ydof, zdof);
          params[xdof] = w * sampler(mt) * max_displacement_x;
          params[ydof] = w * sampler(mt) * max_displacement_y;
          if (ndim > 2) params[zdof] = w * sampler(mt) * max_displacement_z;
        }
      }
    } else {
      for (int l = 0; l < affd->T(); ++l)
      for (int k = 0; k < affd->Z(); ++k)
      for (int j = 0; j < affd->Y(); ++j)
      for (int i = 0; i < affd->X(); ++i) {
        cp = affd->LatticeToIndex(i, j, k, l);
        affd->IndexToDOFs(cp, xdof, ydof, zdof);
        params[xdof] = sampler(mt) * max_displacement_x;
        params[ydof] = sampler(mt) * max_displacement_y;
        if (ndim > 2) params[zdof] = sampler(mt) * max_displacement_z;
      }
    }
    affd->Put(params.data());
  }

  // Print transformation info
  if (verbose) {
    dof->Print();
  }

  // Write transformation to output file
  dof->Write(dofout_name);

  return 0;
}
