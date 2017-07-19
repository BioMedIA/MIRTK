/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2017 Andreas Schuh
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
#include "mirtk/BaseImage.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input>... <output> [options]\n";
  cout << "       " << name << " <input>... -output <output> [options]\n";
  cout << "       " << name << " -input <volume> -image <slice> -output <volume> [options]\n";
  cout << "       " << name << " [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Concatenate two or more either 2D images to form a 3D volume,\n";
  cout << "  or 3D volumes to form a 3D+t temporal sequence. All input images\n";
  cout << "  must have the same image attributes, except in either the third (2D)\n";
  cout << "  or the third and fourth (3D) image dimension.\n";
  cout << "\n";
  cout << "  Moreover, given an :option:`-input` volume (sequence), additional\n";
  cout << "  slices (volumes) can be appended to it. Note that when :option:`-sort`\n";
  cout << "  is enabled, the additional slices (volumes) can be interleaved with the\n";
  cout << "  existing volume (sequence), allowing inserting a slice (volume) anywhere\n";
  cout << "  not only after the last slice (volume) of the input volume (sequence).\n";
  cout << "\n";
  cout << "  Note that the slice thickness of the output volume when concatenating slices\n";
  cout << "  is set equal the average slice thickness of the input images, unless it is\n";
  cout << "  overridden by the :option:`-spacing` value.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -input <file>       Read input -images from given volume/sequence.\n";
  cout << "  -output <path>      Output image file path. Last positional argument when option missing.\n";
  cout << "  -image <file>       Add one input image. Option can be given multiple times.\n";
  cout << "  -images <file>...   Add one or more input images. Option can be given multiple times.\n";
  cout << "  -sort [on|off]      Automatically determine correct order of input images\n";
  cout << "                      based on the distance along the third or fourth image axis.\n";
  cout << "  -nosort             Disable sorting of input images, same as :option:`-sort` off. (default)\n";
  cout << "  -start <float>      Position of first slice/frame of output volume/sequence.\n";
  cout << "  -spacing <float>    Slice/Temporal spacing of output volume/sequence.\n";
  PrintStandardOptions(cout);
  cout << "\n";
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
BaseImage *GetRegion(const UniquePtr<BaseImage> &input, int i)
{
  int i1 = 0, i2 = input->X();
  int j1 = 0, j2 = input->Y();
  int k1 = 0, k2 = input->Z();
  int l1 = 0, l2 = 1;
  if (input->T() > 1) {
    l1 = i;
    l2 = i + 1;
  } else {
    k1 = i;
    k2 = i + 1;
  }
  BaseImage *region = nullptr;
  input->GetRegion(region, i1, j1, k1, l1, i2, j2, k2, l2);
  return region;
}

// -----------------------------------------------------------------------------
bool CheckOrientation(const ImageAttributes &a, const ImageAttributes &b)
{
  if ((a._xaxis[0] * b._xaxis[0] + a._xaxis[1] * b._xaxis[1] + a._xaxis[2] * b._xaxis[2] - 1.0) > .001) return false;
  if ((a._yaxis[0] * b._yaxis[0] + a._yaxis[1] * b._yaxis[1] + a._yaxis[2] * b._yaxis[2] - 1.0) > .001) return false;
  if ((a._zaxis[0] * b._zaxis[0] + a._zaxis[1] * b._zaxis[1] + a._zaxis[2] * b._zaxis[2] - 1.0) > .001) return false;
  return true;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Parse arguments
  REQUIRES_POSARGS(0);

  Array<const char *> input_names;
  for (int i = 1; i <= NUM_POSARGS; ++i) {
    input_names.push_back(POSARG(i));
  }
  const char *input_name  = nullptr;
  const char *output_name = nullptr;

  bool   sort_input = false;
  double spacing    = NaN;
  double start      = NaN;

  for (ALL_OPTIONS) {
    if (OPTION("-image")) {
      input_names.push_back(ARGUMENT);
    }
    else if (OPTION("-images")) {
      do {
        input_names.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-input")) {
      input_name = ARGUMENT;
    }
    else if (OPTION("-output")) {
      output_name = ARGUMENT;
    }
    else if (OPTION("-sort")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(sort_input);
      else sort_input = true;
    }
    else if (OPTION("-nosort")) sort_input = false;
    else if (OPTION("-spacing")) PARSE_ARGUMENT(spacing);
    else if (OPTION("-start")) PARSE_ARGUMENT(start);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (!output_name) {
    if (input_names.empty()) {
      FatalError("No input images and -output file specified!");
    }
    output_name = input_names.back();
    input_names.pop_back();
  }

  int twod = -1;   // concatenate 0: volumes 1: slices
  int m = 0;      // sub-images coming from -input image
  int n = static_cast<int>(input_names.size());
  if (n < 1 || (!input_name && n < 2)) {
    FatalError("Not enough input images given!");
  }
  InitializeIOLibrary();
  UniquePtr<BaseImage> input;
  if (input_name) {
    input.reset(BaseImage::New(input_name));
    if (input->T() > 1) {
      m = input->T();
      twod = 0;
    } else {
      m = input->Z();
      twod = 1;
    }
    n += m;
  }

  // Get attributes of first input image
  UniquePtr<BaseImage> ref;
  if (input) {
    ref.reset(GetRegion(input, 0));
  } else {
    ref.reset(BaseImage::New(input_names[0]));
    if (ref->T() > 1) {
      FatalError("Input images must be either 2D (nz == 1) or 3D (nz == 1 and nt == 1)!");
    }
    twod = (ref->Z() == 1 ? 1 : 0);
  }

  // Sort input images
  //
  // A negative one-based index is associated with the sub-images of the -input image.
  // A positive zero-based index is associated with the (positional) -images.
  //
  // Note: Images are not first read all at once into memory to
  //       not be limited by the amount of memory we have available.
  //       We better read in the images twice from disk if necessary.
  Array<int> pos(n);
  for (int i = 0; i < m; ++i) pos[i] = - (i + 1);
  for (int i = m; i < n; ++i) pos[i] = i - m;

  if (sort_input) {
    if (verbose) {
      cout << "Sorting input images...";
      cout.flush();
    }
    // Compute distance to axis origin
    Point p;
    Array<double> origin(n);
    if (twod) {
      p = ref->GetOrigin();
      ref->WorldToImage(p);
      origin[0] = p._z;
    } else {
      origin[0] = ref->GetTOrigin();
    }
    UniquePtr<BaseImage> image;
    for (int i = 1; i < n; ++i) {
      if (i < m) {
        image.reset(GetRegion(input, i));
      } else {
        image.reset(BaseImage::New(input_names[i]));
      }
      p = image->GetOrigin();
      if (twod) {
        Point p = image->GetOrigin();
        ref->WorldToImage(p);
        origin[i] = p._z;
      } else {
        origin[i] = image->GetTOrigin();
      }
    }
    image.reset();
    // Determine order of images sorted by ascending distance value
    const int refpos = pos[0];
    pos = IncreasingOrder(origin);
    if (pos[0] != refpos) {
      if (pos[0] < 0) {
        ref.reset(GetRegion(input, -pos[0] - 1));
      } else {
        ref.reset(BaseImage::New(input_names[pos[0]]));
      }
    }
    if (verbose) {
      cout << " done\n";
    }
  }

  // Create output image
  if (verbose) {
    cout << "Making ";
    if (twod) cout << "volume";
    else      cout << "sequence";
    cout << " from " << n << " images";
    if (m > 0) {
      cout << ", including " << m << " ";
      if (twod) cout << "slice";
      else      cout << "volume";
      if (m > 1) cout << "s";
      cout << " from input ";
      if (twod) cout << "volume";
      else      cout << "sequence";
    }
    cout << "...";
    cout.flush();
  }

  ImageAttributes out_attr = ref->Attributes();
  if (twod) out_attr._z = n; // spacing and origin set later below
  else      out_attr._t = n;
  UniquePtr<BaseImage> output(BaseImage::New(ref->GetDataType()));
  output->Initialize(out_attr);

  // Concatenate images
  double sum_spacing = 0.;
  double min_zorigin = ref->GetOrigin()._z;

  if (twod) {

    sum_spacing = ref->ZSize();
    for (int j = 0; j < out_attr._y; ++j)
    for (int i = 0; i < out_attr._x; ++i) {
      output->PutAsDouble(i, j, 0, 0, ref->GetAsDouble(i, j));
    }
    UniquePtr<BaseImage> image;
    for (int k = 1; k < n; ++k) {
      if (pos[k] < 0) {
        image.reset(GetRegion(input, -pos[k] - 1));
      } else {
        image.reset(BaseImage::New(input_names[pos[k]]));
      }
      if (image->Z() != 1 || image->T() != 1) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << k+1 << " is not 2D!");
      }
      if (image->X() != ref->X() || image->Y() != ref->Y()) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << k+1 << " has differing number of voxels!");
      }
      if (!AreEqual(image->XSize(), ref->XSize()) ||
          !AreEqual(image->YSize(), ref->YSize())) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << k+1 << " has differing voxel size!");
      }
      if (!AreEqual(image->TSize(), ref->TSize())) {
        Warning("Input image " << k+1 << " has differing temporal spacing.");
      }
      if (!CheckOrientation(image->Attributes(), ref->Attributes())) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << k+1 << " has different orientation!");
      }
      for (int j = 0; j < out_attr._y; ++j)
      for (int i = 0; i < out_attr._x; ++i) {
        output->PutAsDouble(i, j, k, 0, image->GetAsDouble(i, j));
      }
      min_zorigin = min(min_zorigin, image->GetOrigin()._z);
      sum_spacing += image->ZSize();
    }

  } else {

    sum_spacing = ref->TSize();
    for (int k = 0; k < out_attr._z; ++k)
    for (int j = 0; j < out_attr._y; ++j)
    for (int i = 0; i < out_attr._x; ++i) {
      output->PutAsDouble(i, j, k, 0, ref->GetAsDouble(i, j, k));
    }
    UniquePtr<BaseImage> image;
    for (int l = 1; l < n; ++l) {
      if (pos[l] < 0) {
        image.reset(GetRegion(input, -pos[l] - 1));
      } else {
        image.reset(BaseImage::New(input_names[pos[l]]));
      }
      if (image->T() != 1) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << l+1 << " is not 3D!");
      }
      if (image->X() != ref->X() || image->Y() != ref->Y() || image->Z() != ref->Z()) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << l+1 << " has differing number of voxels!");
      }
      if (!AreEqual(image->XSize(), ref->XSize()) ||
          !AreEqual(image->YSize(), ref->YSize()) ||
          !AreEqual(image->ZSize(), ref->ZSize())) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << l+1 << " has differing voxel size!");
      }
      if (!CheckOrientation(image->Attributes(), ref->Attributes())) {
        if (verbose) cout << " failed" << endl;
        FatalError("Input image " << l+1 << " has different orientation!");
      }
      for (int k = 0; k < out_attr._z; ++k)
      for (int j = 0; j < out_attr._y; ++j)
      for (int i = 0; i < out_attr._x; ++i) {
        output->PutAsDouble(i, j, k, l, image->GetAsDouble(i, j, k));
      }
      sum_spacing += image->TSize();
    }

  }

  // Set output spacing to average of input images by default
  if (IsNaN(spacing)) {
    spacing = sum_spacing / n;
  }

  if (twod) {
    if (spacing < 0.) { // must be positive
      out_attr._zaxis[0] *= -1.;
      out_attr._zaxis[1] *= -1.;
      out_attr._zaxis[2] *= -1.;
      spacing            *= -1.;
    }
    double dx, dy, dz;
    output->GetPixelSize(dx, dy, dz);
    output->PutPixelSize(dx, dy, spacing);
  } else {
    output->PutTSize(spacing); // can be negative
  }

  // Set output origin
  if (twod) {
    if (IsNaN(start)) {
      start = min_zorigin;
    }
    Point origin = output->GetOrigin();
    origin._z = start;
    // Move to center of image stack, i.e. translate along z axis
    double s = .5 * (n - 1) * output->ZSize();
    origin._x += s * output->Attributes()._zaxis[0];
    origin._y += s * output->Attributes()._zaxis[1];
    origin._z += s * output->Attributes()._zaxis[2];
    output->PutOrigin(origin);
  } else {
    if (!IsNaN(start)) out_attr._torigin = start;
    // else torigin is time of first (sorted) input volume
  }

  if (verbose) cout << " done" << endl;

  // Write output image
  if (verbose) {
    cout << "Writing ";
    cout << (sort_input ? "combined" : "concatenated");
    cout << " images to " << output_name << "...";
    cout.flush();
  }
  output->Write(output_name);
  if (verbose) {
    cout << " done" << endl;
  }

  return 0;
}

