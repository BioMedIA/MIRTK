/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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
#include "mirtk/Transformation.h"
#include "mirtk/RegisteredImage.h"
#include "mirtk/EnergyMeasure.h"
#include "mirtk/SimilarityMeasure.h"
#include "mirtk/ImageSimilarity.h"
#include "mirtk/HistogramImageSimilarity.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"
#include "mirtk/HomogeneousTransformation.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <target> <source>... [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Evaluates the (dis-)similarity of two intensity images.\n";
  cout << "\n";
  cout << "  If more than one source image is given, the (dis-)similarity between\n";
  cout << "  each of these and the target is evaluated. By default, common image\n";
  cout << "  (dis-)similarity metrics are reported. One or more metrics for the\n";
  cout << "  evaluation can be chosen using :option:`-metric`.\n";
  cout << "\n";
  cout << "  The input source image can either be pre-aligned with the target image\n";
  cout << "  or the output of a previous registration can be specified using\n";
  cout << "  :option:`-dofin`.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  target   Target image.\n";
  cout << "  source   Source image.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -mask <file>       Target image region of interest mask. (default: :option:`-Tp`)\n";
  cout << "  -dofin <file>      Source image transformation. (default: Id)\n";
  cout << "  -interp <mode>     Interpolation mode. (default: linear with padding)\n";
  cout << "  -metric <sim>...   Image (dis-)similarity measure(s):\n";
  cout << "                     ";
  EnergyMeasure m = static_cast<EnergyMeasure>(SIM_Begin + 1);
  while (m < SIM_End) {
    if (m > SIM_Begin + 1) cout << ", ";
    cout << ToString(m);
    m = static_cast<EnergyMeasure>(m + 1);
  }
  cout << "\n";
  cout << "  -p, -padding <value>   Common image background value/threshold. (default: NaN)\n";
  cout << "  -Tp <value>        Target image background value/threshold. (default: NaN)\n";
  cout << "  -Sp <value>        Source image background value/threshold. (default: NaN)\n";
  cout << "  -bins <int>        Number of histogram bins.        (default: min(dynamic range, 255))\n";
  cout << "  -Tbins <int>       Number of target histogram bins. (default: min(dynamic range, 255))\n";
  cout << "  -Sbins <int>       Number of source histogram bins. (default: min(dynamic range, 255))\n";
  cout << "  -parzen            Smooth histogram samples using cubic B-spline Parzen window. (default: off)\n";
  cout << "  -Rx1 <int>         Leftmost  target voxel index along x axis.\n";
  cout << "  -Rx2 <int>         Rightmost target voxel index along x axis.\n";
  cout << "  -Ry1 <int>         Leftmost  target voxel index along y axis.\n";
  cout << "  -Ry2 <int>         Rightmost target voxel index along y axis.\n";
  cout << "  -Rz1 <int>         Leftmost  target voxel index along z axis.\n";
  cout << "  -Rz2 <int>         Rightmost target voxel index along z axis.\n";
  cout << "\n";
  cout << "Output options:\n";
  cout << "  -precision <int>   Number of significant digits. (default: 5)\n";
  cout << "  -delim <char>      Delimiter for output of multiple metric values. (default: ,)\n";
  cout << "  -table             Output in tabular format.\n";
  cout << "  -csv               Output as comma separated values table.\n";
  cout << "  -tsv               Output as tab   separated values table.\n";
  cout << "  -noid              Exclude image IDs from output. (default: off)\n";
  cout << "  -fullid            Use complete input image file path as ID.\n";
  cout << "                     (default: file name without image extension)\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

typedef RegisteredImage::DisplacementImageType DisplacementField;

// -----------------------------------------------------------------------------
int DefaultNumberOfBins(const BaseImage *, double min_intensity, double max_intensity)
{
  int nbins = iround(max_intensity - min_intensity) + 1;
  if (nbins > 256) nbins = 256;
  return nbins;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  typedef HistogramImageSimilarity::JointHistogramType JointHistogram;
  typedef RegisteredImage::VoxelType                       VoxelType;

  // Check command line
  REQUIRES_POSARGS(2);

  // Parse positional arguments
  const char *target_name = POSARG(1);
  Array<const char *> source_name;
  source_name.reserve(NUM_POSARGS-1);
  source_name.push_back(POSARG(2));
  for (OPTIONAL_POSARGS) {
    source_name.push_back(POSARG(ARGIDX));
  }

  // Parse optional arguments
  Array<SimilarityMeasure> metric;
  InterpolationMode interp   = Interpolation_LinearWithPadding;
  const char *mask_name      = nullptr;
  const char *dofin_name     = nullptr;
  const char *table_name     = nullptr;
  bool        invert         = false;
  double      target_padding = NaN;
  double      source_padding = NaN;
  int         target_bins    = 0;
  int         source_bins    = 0;
  bool        parzen_window  = false;
  bool        full_path_id   = false;
  bool        image_id_col   = true;
  bool        table_header   = true;
  int         digits         = 5;
  string      delim;

  int i1 =  0, j1 =  0, k1 =  0;
  int i2 = -1, j2 = -1, k2 = -1;

  for (ALL_OPTIONS) {
    if      (OPTION("-dofin")) dofin_name = ARGUMENT;
    else if (OPTION("-mask"))  mask_name  = ARGUMENT;
    else if (OPTION("-interp")) PARSE_ARGUMENT(interp);
    else if (OPTION("-Tp")) PARSE_ARGUMENT(target_padding);
    else if (OPTION("-Sp")) PARSE_ARGUMENT(source_padding);
    else if (OPTION("-padding") || OPTION("-p")) {
      PARSE_ARGUMENT(target_padding);
      source_padding = target_padding;
    }
    else if (OPTION("-tbins") || OPTION("-Tbins")) PARSE_ARGUMENT(target_bins);
    else if (OPTION("-sbins") || OPTION("-Sbins")) PARSE_ARGUMENT(target_bins);
    else if (OPTION("-bins")) {
      PARSE_ARGUMENT(target_bins);
      source_bins = target_bins;
    }
    else if (OPTION("-parzen") || OPTION("-parzen-window") || OPTION("-smooth-histogram")) {
      parzen_window = true;
    }
    else if (OPTION("-metric") || OPTION("-measure")) {
      do {
        SimilarityMeasure m;
        PARSE_ARGUMENT(m);
        metric.push_back(m);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-precision")) PARSE_ARGUMENT(digits);
    else if (OPTION("-delim")) {
      delim = ARGUMENT;
      size_t pos;
      while ((pos = delim.find("\\n")) != string::npos) {
        delim.replace(pos, 2, "\n");
      }
      while ((pos = delim.find("\\t")) != string::npos) {
        delim.replace(pos, 2, "\t");
      }
    }
    else if (OPTION("-table")) {
      verbose = -1;
      if (HAS_ARGUMENT) {
        table_name = ARGUMENT;
      }
    }
    else if (OPTION("-csv")) verbose = -1, delim = ',';
    else if (OPTION("-tsv")) verbose = -1, delim = '\t';
    else if (OPTION("-fullid")) full_path_id = true;
    else if (OPTION("-Rx1")) PARSE_ARGUMENT(i1);
    else if (OPTION("-Rx2")) PARSE_ARGUMENT(i2);
    else if (OPTION("-Ry1")) PARSE_ARGUMENT(j1);
    else if (OPTION("-Ry2")) PARSE_ARGUMENT(j2);
    else if (OPTION("-Rz1")) PARSE_ARGUMENT(k1);
    else if (OPTION("-Rz2")) PARSE_ARGUMENT(k2);
    else HANDLE_BOOL_OPTION(invert);
    else HANDLE_BOOLEAN_OPTION("id", image_id_col);
    else HANDLE_BOOLEAN_OPTION("header", table_header);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (metric.empty()) {
    metric.reserve(9);
    metric.push_back(SIM_SSD);
    metric.push_back(SIM_PSNR);
    metric.push_back(SIM_CoVar);
    metric.push_back(SIM_CC);
    metric.push_back(SIM_JE);
    metric.push_back(SIM_MI);
    metric.push_back(SIM_NMI);
    metric.push_back(SIM_CR_XY);
    metric.push_back(SIM_CR_YX);
  }

  ofstream fout;
  ostream *out = &cout;
  if (table_name) {
    verbose = -1;
    fout.open(table_name, ios::out);
    if (!fout) {
      FatalError("Failed to open output file: " << table_name);
    }
    out = &fout;
  }

  int width = 22;
  if (verbose > 0) {
    delim = '\n';
    for (size_t i = 0; i < metric.size(); ++i) {
      width = max(width, static_cast<int>(ToPrettyString(metric[i]).length()));
    }
  } else if (delim.empty()) {
    delim = ',';
    if (verbose >= 0) delim += ' ';
  }

  // Initialize I/O module
  InitializeIOLibrary();

  // Input and (transformed) image instances
  RegisteredImage::InputImageType target_image, source_image;
  RegisteredImage target, source;

  target.InputImage(&target_image);
  source.InputImage(&source_image);
  source.InterpolationMode(interp);

  // Registered images must be updated each time because ImageSimilarity subclasses
  // are allowed to change the intensities of the registered images. For example,
  // the NormalizedIntensityCrossCorrelation normalizes the values to [0, 1].
  target.SelfUpdate(true);
  source.SelfUpdate(true);

  // Read target image
  if (verbose > 1) cout << "Reading target image from " << target_name << "...", cout.flush();
  target_image.Read(target_name);
  target_image.PutBackgroundValueAsDouble(target_padding, true);
  if (target_image.T() > 1) {
    if (verbose > 1) cout << " failed" << endl;
    FatalError("Target image must be a 2D or 3D image!");
  }
  if (verbose > 1) cout << " done" << endl;

  // Target region of interest
  if (i2 < 0) i2 = target_image.X();
  if (j2 < 0) j2 = target_image.Y();
  if (k2 < 0) k2 = target_image.Z();

  ImageAttributes target_region = target_image.Attributes();
  target_region._x       = i2 - i1;
  target_region._y       = j2 - j1;
  target_region._z       = k2 - k1;
  target_region._xorigin = .0;
  target_region._yorigin = .0;
  target_region._zorigin = .0;

  // Adjust origin
  double x1 = i1, y1 = j1, z1 = k1;
  target_image.ImageToWorld(x1, y1, z1);
  double x2 = .0, y2 = .0, z2 = .0;
  target_region.LatticeToWorld(x2, y2, z2);
  target_region._xorigin = x1 - x2;
  target_region._yorigin = y1 - y2;
  target_region._zorigin = z1 - z2;

  const int nvox = target_region.NumberOfPoints();

  // Read input (and evaluate) source transformation
  // ATTENTION: MUST be done **before** source.Initialize() is called!
  UniquePtr<Transformation> dofin;
  UniquePtr<DisplacementField> disp;
  if (dofin_name) {
    if (verbose > 1) cout << "Reading source transformation from " << dofin_name << "...", cout.flush();
    dofin.reset(Transformation::New(dofin_name));
    source.Transformation(dofin.get());
    HomogeneousTransformation *lin;
    lin = dynamic_cast<HomogeneousTransformation *>(dofin.get());
    if (invert && lin) {
      lin->Invert();
    } else if (!lin && (invert || dofin->RequiresCachingOfDisplacements())) {
      disp.reset(new DisplacementField(target_region, 3));
      const double t  = source.GetTOrigin();
      const double t0 = target.GetTOrigin();
      if (invert) {
        dofin->InverseDisplacement(*disp, t, t0);
      } else {
        dofin->Displacement(*disp, t, t0);
      }
      source.ExternalDisplacement(disp.get());
    }
    if (verbose > 1) cout << " done" << endl;
  }

  // Initialize (transformed) images
  if (verbose > 1) cout << "Initializing registered images...", cout.flush();
  target.Initialize(target_region);
  source.Initialize(target_region);
  if (IsNaN(target_padding)) {
    target.ClearBackgroundValue();
  } else {
    target.PutBackgroundValueAsDouble(target_padding);
  }
  if (IsNaN(source_padding)) {
    source.ClearBackgroundValue();
  } else {
    source.PutBackgroundValueAsDouble(source_padding);
  }
  target.Recompute();
  if (verbose > 1) cout << " done" << endl;

  if (debug) {
    target.Write("debug_target.nii.gz");
  }

  // Target region of interest and image overlap masks
  BinaryImage mask(target_region), overlap(target_region);

  if (mask_name) {
    if (verbose > 1) cout << "Reading target mask from " << mask_name << "...", cout.flush();
    BinaryImage input_mask(mask_name);
    if (input_mask.T() > 1) {
      FatalError("Mask image must be a 2D or 3D scalar image!");
    }
    if (input_mask.HasSpatialAttributesOf(&target)) {
      mask.CopyFrom(input_mask);
    } else {
      GenericNearestNeighborInterpolateImageFunction<BinaryImage> nn;
      nn.Input(&input_mask);
      nn.Initialize();
      double x, y, z;
      for (int k = 0; k < target.Z(); ++k)
      for (int j = 0; j < target.Y(); ++j)
      for (int i = 0; i < target.X(); ++i) {
        x = i, y = j, z = k;
        target.ImageToWorld(x, y, z);
        mask(i, j, k) = static_cast<BinaryPixel>(nn.Evaluate(x, y, z));
      }
    }
    if (verbose > 1) cout << " done" << endl;
  } else {
    const VoxelType *tgt = target.Data();
    if (IsNaN(target_padding)) {
      for (int idx = 0; idx < nvox; ++idx, ++tgt) {
        mask(idx) = (IsNaN(*tgt) ? 0 : 1);
      }
    } else {
      for (int idx = 0; idx < nvox; ++idx, ++tgt) {
        mask(idx) = (*tgt > target_padding ? 1 : 0);
      }
    }
  }

  target.PutMask(&mask);

  // Target intensity range in region of interest
  double tmin, tmax;
  target_image.GetMinMaxAsDouble(tmin, tmax);
  if (fequal(tmin, tmax)) {
    FatalError("Target image has homogeneous intensity = " << tmin);
  }

  // Initialize similarity measures
  if (verbose > 1) cout << "Initializing similarity measures...", cout.flush();
  Array<UniquePtr<ImageSimilarity> > sim(metric.size());
  bool use_shared_histogram = false;
  JointHistogram samples;

  for (size_t i = 0; i < metric.size(); ++i) {
    sim[i].reset(ImageSimilarity::New(metric[i], ToString(metric[i]).c_str(), 1.0));
    sim[i]->Target(&target);
    sim[i]->Source(&source);
    sim[i]->Mask(&overlap);
    sim[i]->Domain(target_region);
    sim[i]->DivideByInitialValue(false);
    sim[i]->SkipTargetInitialization(true);
    sim[i]->SkipSourceInitialization(true);
    HistogramImageSimilarity *p;
    p = dynamic_cast<HistogramImageSimilarity *>(sim[i].get());
    if (p != nullptr) {
      p->Samples(&samples);
      p->UseParzenWindow(parzen_window);
      use_shared_histogram = true;
    }
  }
  if (verbose > 1) cout << " done" << endl;

  // Print table header
  const string target_id = (full_path_id ? target_name : FileName(target_name));

  if (table_header && verbose < 0) {
    if (image_id_col) {
      *out << "Target" << delim << "Source" << delim;
    }
    for (size_t i = 0; i < metric.size(); ++i) {
      if (i > 0) *out << delim;
      *out << ToString(metric[i]);
    }
    *out << endl;
  }

  // Evaluate similarity of (transformed) source image(s)
  for (size_t n = 0; n < source_name.size(); ++n) {

    Indent indent;
    if (verbose > 0 && n > 0) cout << "\n";

    // Read source image
    if (verbose > 1) {
      cout << "Reading source image from " << source_name[n] << "...";
      cout.flush();
    }
    source_image.Read(source_name[n]);
    source_image.PutBackgroundValueAsDouble(source_padding, true);
    if (verbose > 1) {
      cout << " done" << endl;
      if (source.Transformation()) {
        cout << "Transforming source image...";
        cout.flush();
      }
    }
    source.Recompute();
    if (verbose > 1 && source.Transformation()) {
      cout << " done\n" << endl;
    }
    if (debug) {
      if (source_name.size() == 1) {
        source.Write("debug_source.nii.gz");
      } else {
        source.Write(("debug_source_" + ToString(n + 1) + ".nii.gz").c_str());
      }
    }

    // Print image IDs
    const string source_id = (full_path_id ? source_name[n] : FileName(source_name[n]));

    if (verbose > 0) {
      cout << "Similarity of target image " << target_id << " and source image " << source_id << ":\n";
      ++indent;
    } else if (image_id_col) {
      if (verbose == 0) cout << "Target = ";
      *out << target_id << delim;
      if (verbose == 0) cout << "Source = ";
      *out << source_id;
    }
    cout.flush();

    // Compute overlap mask
    overlap = mask;
    for (int idx = 0; idx < nvox; ++idx) {
      overlap(idx) = (overlap(idx) && source.IsForeground(idx) ? 1 : 0);
    }

    // Fill joint histogram
    if (use_shared_histogram) {
      double smin, smax;
      if (n > 0) target.Recompute();
      source_image.GetMinMaxAsDouble(smin, smax);
      if (fequal(smin, smax)) {
        FatalError("Source image has homogeneous intensity = " << smin);
      }
      int tbins = target_bins, sbins = source_bins;
      if (tbins <= 0) tbins = DefaultNumberOfBins(&target_image, tmin, tmax);
      if (sbins <= 0) sbins = DefaultNumberOfBins(&source_image, smin, smax);
      double twidth = (tmax - tmin) / tbins;
      double swidth = (smax - smin) / sbins;
      samples.Initialize(tmin, tmax, twidth, smin, smax, swidth);
      const VoxelType *tgt = target.Data();
      const VoxelType *src = source.Data();
      for (int idx = 0; idx < nvox; ++idx, ++tgt, ++src) {
        if (mask(idx)) {
          samples.Add(samples.ValToBinX(*tgt), samples.ValToBinY(*src));
        }
      }
      if (verbose > 0) {
        cout.precision(digits);
        #define Print(name, value) \
          cout << delim << indent << ToString(name, width, ' ', true) << " = " << value
        Print("No. of samples",         samples.NumberOfSamples());
        Print("No. of histogram bins",  samples.NumberOfBinsX() << " x " << samples.NumberOfBinsY());
        Print("Size of histogram bins", samples.WidthX() << " x " << samples.WidthY());
        Print("Mean target intensity",  samples.MeanX());
        Print("Mean source intensity",  samples.MeanY());
        Print("Target image range", "[" << samples.MinX() << ", " << samples.MaxX() << "]");
        Print("Source image range", "[" << samples.MinY() << ", " << samples.MaxY() << "]");
        Print("Target image variance",  samples.VarianceX());
        Print("Source image variance",  samples.VarianceY());
        Print("Target image entropy",   samples.EntropyX());
        Print("Source image entropy",   samples.EntropyY());
        #undef Print
        cout.flush();
      }
    }

    // Evaluate similarity measure(s)
    for (size_t i = 0; i < sim.size(); ++i) {
      if (image_id_col || i > 0) {
        *out << delim;
      }
      if (verbose > 0) {
        cout << indent << ToPrettyString(metric[i], width, ' ', true) << " = ";
      } else if (verbose == 0) {
        cout << sim[i]->Name() << " = ";
      }
      if (verbose > 0) cout.flush();
      // Initialize similarity measure, possibly changes [min, max] range
      sim[i]->Initialize();
      sim[i]->Update();
      // Evaluate image similarity
      *out << setprecision(digits) << sim[i]->RawValue();
      // Reset registered image properties which may be modified by similarity measures
      target.MinIntensity(NaN);
      target.MaxIntensity(NaN);
      source.MinIntensity(NaN);
      source.MaxIntensity(NaN);
      if (IsNaN(target_padding)) {
        target.ClearBackgroundValue();
      } else {
        target.PutBackgroundValueAsDouble(target_padding);
      }
      if (IsNaN(source_padding)) {
        source.ClearBackgroundValue();
      } else {
        source.PutBackgroundValueAsDouble(source_padding);
      }
      if (verbose > 0) cout.flush();
    }

    *out << endl;
  }

  // Close output file
  if (fout.is_open()) {
    fout.close();
  }
  return 0;
}
