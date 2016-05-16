/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
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

#include "mirtk/GenericImage.h"
#include "mirtk/LieBracketImageFilter.h"

using namespace mirtk;


// ===========================================================================
// Help
// ===========================================================================

// ---------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <x> <y> <z>                                        [options]" << endl;
  cout << "       " << name << " <xx> <xy> [<xz>] <yx> <yy> [<yz>] <zx> <zy> [<zz>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes the Lie bracket of two vector fiels :math:`z = [x, y]`." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -jac   Use Jacobian definition of the Lie bracket. (default)" << endl;
  cout << "  -dif   Use difference of compositions: :math:`z = x(y) - y(x)`." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// ===========================================================================
// Image type
// ===========================================================================

typedef GenericImage<double> ImageType;

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Print help or version and exit if requested
  HANDLE_HELP_OR_VERSION();

  InitializeIOLibrary();

  // Parse positional arguments
  const char *x_fname = NULL;
  const char *y_fname = NULL;
  const char *z_fname = NULL;
  const char *xx_fname   = NULL;
  const char *xy_fname   = NULL;
  const char *xz_fname   = NULL;
  const char *yx_fname   = NULL;
  const char *yy_fname   = NULL;
  const char *yz_fname   = NULL;
  const char *zx_fname   = NULL;
  const char *zy_fname   = NULL;
  const char *zz_fname   = NULL;

  int posargc = 0;
  while (posargc+1 < argc && argv[posargc+1][0] != '-') posargc++;

  switch (posargc) {
    case 3:
      x_fname = argv[1];
      y_fname = argv[2];
      z_fname = argv[3];
      break;
    case 6:
      xx_fname = argv[1];
      xy_fname = argv[2];
      yx_fname = argv[3];
      yy_fname = argv[4];
      zx_fname = argv[5];
      zy_fname = argv[6];
      break;
    case 9:
      xx_fname = argv[1];
      xy_fname = argv[2];
      xz_fname = argv[3];
      yx_fname = argv[4];
      yy_fname = argv[5];
      yz_fname = argv[6];
      zx_fname = argv[7];
      zy_fname = argv[8];
      zz_fname = argv[9];
      break;
    default:
      PrintHelp(EXECNAME);
      exit(1);
  }

  // Parse optional arguments
  bool usejac = true;

  for (ARGUMENTS_AFTER(posargc)) {
    if      (OPTION("-jac")) usejac = true;
    else if (OPTION("-dif")) usejac = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input vector fields
  UniquePtr<ImageType> x, y;

  if (x_fname) {
    x.reset(new ImageType(x_fname));
  } else {
    ImageType xx(xx_fname);
    ImageType xy(xy_fname);
    UniquePtr<ImageType> xz(xz_fname ? new ImageType(xz_fname) : NULL);

    ImageAttributes attr = xx.Attributes();
    if (attr._t > 1) {
      cerr << "Separate input component images must be scalar fields." << endl;
      cerr << "Try calling this program with single vector field input and output images instead." << endl;
      exit(1);
    }
    if (xy.Attributes() != attr || (xz.get() && xz->Attributes() != attr)) {
      cerr << "Image attributes of input vector component images do not match!" << endl;
      exit(1);
    }
    attr._t = (xz.get() ? 3 : 2);

    x.reset(new ImageType(attr));
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      x->Put(i, j, k, 0, xx.Get(i, j, k));
      x->Put(i, j, k, 1, xy.Get(i, j, k));
      if (xz.get()) x->Put(i, j, k, 2, xz->Get(i, j, k));
    }
  }

  if (y_fname) {
    y.reset(new ImageType(y_fname));
  } else {
    ImageType yx(yx_fname);
    ImageType yy(yy_fname);
    UniquePtr<Image> yz(yz_fname ? new ImageType(yz_fname) : NULL);

    ImageAttributes attr = yx.Attributes();
    if (attr._t > 1) {
      cerr << "Separate input component images must be scalar fields." << endl;
      cerr << "Try calling this program with single vector field input and output images instead." << endl;
      exit(1);
    }
    if (yy.Attributes() != attr || (yz.get() && yz->Attributes() != attr)) {
      cerr << "Image attributes of input vector component images do not match!" << endl;
      exit(1);
    }
    attr._t = (yz.get() ? 3 : 2);

    y.reset(new ImageType(attr));
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      y->Put(i, j, k, 0, yx.Get(i, j, k));
      y->Put(i, j, k, 1, yy.Get(i, j, k));
      if (yz.get()) y->Put(i, j, k, 2, yz->Get(i, j, k));
    }
  }

  if (!y->HasSpatialAttributesOf(x.get())) {
    cerr << "Spatial attributes of input vector fields do not match!" << endl;
    exit(1);
  }
  if (y->T() != x->T()) {
    cerr << "Number of vector components of input vector fields differ!" << endl;
    exit(1);
  }

  // Compute Lie bracket
  ImageType z(x->Attributes());
  liebracket(&z, x.get(), y.get(), usejac);

  // Save output velocities
  if (z_fname) {
    z.Write(z_fname);
  } else {
    z.GetFrame(0).Write(zx_fname);
    z.GetFrame(1).Write(zy_fname);
    if (zz_fname) z.GetFrame(2).Write(zz_fname);
  }

  // Clean up
  return 0;
}
