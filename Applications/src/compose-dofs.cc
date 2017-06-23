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
#include "mirtk/Transformations.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
// Print help screen
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <T1> <T2>... <T> [-target <image>]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes the composition T of the given input transformations such that\n";
  cout << "\n";
  cout << "  .. math::\n";
  cout << "\n";
  cout << "     T(x) = Tn o ... o T2 o T1(x)\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -target <image>\n";
  cout << "      Target image on which images will be resampled using the composed transformation.\n";
  cout << "      The finite grid of the target image is used to determine an appropriate domain\n";
  cout << "      on which to approximate the displacement field of the composite transformation\n";
  cout << "      when :option:`-approximate` is given.\n";
  cout << "  -[no]rotation\n";
  cout << "      Whether to allow rotation    when composite transformation is affine. (default: on)\n";
  cout << "  -[no]translation\n";
  cout << "      Whether to allow translation when composite transformation is affine. (default: on)\n";
  cout << "  -[no]scaling\n";
  cout << "      Whether to allow scaling when composite transformation is affine. (default: on)\n";
  cout << "  -[no]shearing\n";
  cout << "      Whether to allow shearing when composite transformation is affine. (default: on)\n";
  cout << "  -approximate\n";
  cout << "      Approximate the composed transformation using a single FFD. (default: off)\n";
  cout << "  -scale <s1> [s2...]\n";
  cout << "      Scaling factors for each input (SV) FFD in the same order as the positional\n";
  cout << "      input file name arguments. These factors can only be applied to single FFDs,\n";
  cout << "      and are mainly useful for velocity based transformations. For rigid, similarity,\n";
  cout << "      or affine transformations, use -1 to invert the input transformation. (default: 1)\n";
  cout << "  -bch [<n> [yes|no]]\n";
  cout << "      Use Baker-Campbell-Hausdorff (BCH) formula to approximate composition of SV FFDs.\n";
  cout << "      All input transformations must be of type cubic B-spline SV FFD. Arguments\n";
  cout << "      are optional. The first argument is the number of BCH terms to use. The minimum\n";
  cout << "      is 2 terms, i.e., the sum of left and right velocity fields. The second argument\n";
  cout << "      is a boolean flag indicating whether or not the Lie brackets should be computed\n";
  cout << "      using the Jacobian of the vector fields (yes) or if it should be approximated\n";
  cout << "      as the difference of the compositions in either order. (default: 6 yes)\n";
  cout << "  -nobch\n";
  cout << "      Do not use BCH formula to compose SV FFDs. Instead, evaluate composite displacements\n";
  cout << "      and approximate these even when all input transformations are of type SV FFD.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(3);
  int N = NUM_POSARGS - 1;

  ImageAttributes attr;
  FluidFreeFormTransformation t;
  bool translation    = true;
  bool rotation       = true;
  bool scaling        = true;
  bool shearing       = true;
  bool approximate    = false;
  int  no_bch_terms   = -1;
  bool lie_derivative = true;
  Array<double> scales;
  double dx = 0., dy = 0., dz = 0.;

  for (ALL_OPTIONS) {
    if (OPTION("-target")) {
      InitializeIOLibrary();
      GreyImage target(ARGUMENT);
      attr = target.Attributes();
    }
    else if (OPTION("-scale")) {
      do {
        double s;
        PARSE_ARGUMENT(s);
        if (AreEqual(s, 0.) || IsNaN(s) || IsInf(s)) {
          FatalError("Scaling factor must be neither zero, +/- inf, or NaN!");
        }
        scales.push_back(s);
        if (scales.size() > static_cast<size_t>(N)) {
          FatalError("Too many arguments for option -scale!");
        }
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-bch")) {
      approximate = true;
      no_bch_terms = 4;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(no_bch_terms);
        if (HAS_ARGUMENT) {
          PARSE_ARGUMENT(lie_derivative);
        }
      }
    }
    else if (OPTION("-nobch")) {
      no_bch_terms = 0;
    }
    else if (OPTION("-rigid")) {
      bool bval = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(bval);
      rotation = translation = bval;
    }
    else if (OPTION("-norigid")) {
      rotation = translation = false;
    }
    else if (OPTION("-affine") || OPTION("-global")) {
      bool bval = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(bval);
      rotation = translation = scaling = shearing = bval;
    }
    else if (OPTION("-noaffine") || OPTION("-noglobal")) {
      rotation = translation = scaling = shearing = false;
    }
    else if (OPTION("-spacing")) {
      PARSE_ARGUMENT(dx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(dy);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(dz);
        else dz = 0.;
      } else {
        dy = dz = dx;
      }
    }
    else HANDLE_BOOL_OPTION(translation);
    else HANDLE_BOOL_OPTION(rotation);
    else HANDLE_BOOL_OPTION(scaling);
    else HANDLE_BOOL_OPTION(shearing);
    else HANDLE_BOOL_OPTION(approximate);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (verbose) {
    cout << "Compose remaining transformations..." << endl;
  }
  ImageAttributes disp_attr = attr;
  for (int n = 1; n <= N; ++n) {
    // Read n-th transformation
    if (verbose) {
      cout << "Reading transformation " << n << " from " << POSARG(n) << endl;
    }
    UniquePtr<Transformation> dof(Transformation::New(POSARG(n)));
    // Scale velocities / invert transformation
    if (n <= static_cast<int>(scales.size()) && !AreEqual(scales[n - 1], 1.)) {
      auto aff  = dynamic_cast<HomogeneousTransformation *>(dof.get());
      auto ffd  = dynamic_cast<FreeFormTransformation    *>(dof.get());
      auto mffd = dynamic_cast<MultiLevelTransformation  *>(dof.get());
      if (mffd) {
        if (mffd->NumberOfLevels() == 0) {
          aff = mffd->GetGlobalTransformation();
        } else if (mffd->GetGlobalTransformation()->IsIdentity()) {
          if (mffd->NumberOfLevels() == 1) {
            ffd = mffd->GetLocalTransformation(0);
          }
        }
      }
      if (aff) {
        if (AreEqual(scales[n - 1], -1.)) {
          cout << "  Inverting transformation" << endl;
          aff->Invert();
        } else {
          FatalError("Cannot -scale affine transformation, use '-1' to invert it!");
        }
      } else if (ffd) {
        const double scale = scales[n - 1];
        cout << "  Scaling FFD coefficients by factor " << scale << endl;
        for (int i = 0; i < ffd->NumberOfDOFs(); ++i) {
          ffd->Put(i, scale * ffd->Get(i));
        }
      } else {
        FatalError("Cannot -scale non-FFD transformation! Option mainly used for SV FFDs.");
      }
    }
    // Compose with current chain
    t.PushTransformation(dof.get(), &disp_attr);
  }
  if (verbose) {
    cout << "Compose transformations... done\n" << endl;
  }

  // Approximate the composed transformation using a single free-form deformation
  if (approximate && t.NumberOfLevels() > 1) {
    double rms_error;

    // Get copy of global transformation
    const auto global = t.GetGlobalTransformation()->GetMatrix();

    // Use the most dense control point lattice in the local transformation stack
    ImageAttributes domain = attr;
    if (domain) {
      domain.PutAffineMatrix(global, true);
    } else {
      for (int i = 0; i < t.NumberOfLevels(); ++i) {
        const auto &ffd_attr = t.GetLocalTransformation(i)->Attributes();
        if (ffd_attr.NumberOfLatticePoints() > domain.NumberOfLatticePoints()) {
          domain = ffd_attr;
        }
      }
    }
    domain._t  = 1;
    domain._dt = 0.;

    // Reset global transformation
    t.GetGlobalTransformation()->Reset();

    // Determine if only SV FFDs are being composed and use BCH formula (if not -nobch option given)
    bool common_type_is_svffd = true;
    for (int i = 0; i < t.NumberOfLevels(); ++i) {
      if (t.GetLocalTransformation(i)->TypeOfClass() != TRANSFORMATION_BSPLINE_FFD_SV) {
        common_type_is_svffd = false;
        break;
      }
    }
    if (common_type_is_svffd && no_bch_terms < 0) {
      no_bch_terms = 4;
    } else if (!common_type_is_svffd && no_bch_terms > 0) {
      FatalError("Cannot use Baker-Campbell-Hausdorff (BCH) formula for non-SV FFD transformations");
    }

    // Approximate composite transformation by single (SV) FFD
    UniquePtr<FreeFormTransformation3D> ffd;

    if (no_bch_terms > 0) {

      if (verbose) {
        cout << "Approximate composed SV FFD using Baker-Campbell-Hausdorff (BCH) formula...";
        cout.flush();
      }
      UniquePtr<BSplineFreeFormTransformationSV> svffd;
      auto *first = dynamic_cast<BSplineFreeFormTransformationSV *>(t.GetLocalTransformation(0));
      if (attr) {
        svffd.reset(new BSplineFreeFormTransformationSV());
        svffd->Initialize(domain, dx, dy, dz, first);
        svffd->Parameter(first->Parameter());
      } else {
        svffd.reset(new BSplineFreeFormTransformationSV(*first));
      }
      const auto orig_lie_derivative = svffd->LieDerivative();
      const auto orig_no_bch_terms = svffd->NumberOfBCHTerms();
      svffd->NumberOfBCHTerms(no_bch_terms);
      svffd->LieDerivative(lie_derivative);
      for (int i = 1; i < t.NumberOfLevels(); ++i) {
        auto *next = dynamic_cast<BSplineFreeFormTransformationSV *>(t.GetLocalTransformation(i));
        svffd->CombineWith(next);
      }
      svffd->NumberOfBCHTerms(orig_no_bch_terms);
      svffd->LieDerivative(orig_lie_derivative);
      rms_error = svffd->EvaluateRMSError(domain, &t);
      ffd.reset(svffd.release());

    } else {

      if (verbose) {
        cout << "Approximate composed transformation using a single FFD...";
        cout.flush();
      }
      ffd.reset(new BSplineFreeFormTransformation3D(domain, dx, dy, dz));
      rms_error = ffd->ApproximateAsNew(domain, &t);

    }

    // Set new composite global and local transformations
    t.Clear();
    t.GetGlobalTransformation()->PutMatrix(global);
    t.PushLocalTransformation(ffd.release());

    // Report approximation error
    if (verbose) {
      cout << " done";
      if (verbose > 1) {
        cout << ": RMSE = " << rms_error;
      }
      cout << endl;
    }
  }

  // Merge global and (first) local transformation
  if (t.NumberOfLevels() > 0 && !translation && !rotation && !scaling && !shearing) {
    if (verbose) cout << "Merging global into first local transformation...";
    t.MergeGlobalIntoLocalDisplacement();
    if (verbose) cout << " done" << endl;
  }

  // Write composite transformation
  const char *output_name = POSARG(NUM_POSARGS);
  if (verbose) cout << "Writing composite transformation to " << output_name << endl;
  if (t.NumberOfLevels() == 0) {
    AffineTransformation aff(*t.GetGlobalTransformation());
    aff.PutMatrix(t.GetAffineTransformation()->GetMatrix() * aff.GetMatrix());
    if (!translation) {
      aff.PutTranslationX(0.);
      aff.PutTranslationY(0.);
      aff.PutTranslationZ(0.);
    }
    if (!rotation) {
      aff.PutRotationX(0.);
      aff.PutRotationY(0.);
      aff.PutRotationZ(0.);
    }
    if (!scaling) {
      aff.PutScale(100.);
    }
    if (!shearing) {
      aff.PutShearXY(0.);
      aff.PutShearYZ(0.);
      aff.PutShearXZ(0.);
    }
    aff.Write(output_name);
    if (verbose) aff.Print(2);
  } else if (t.NumberOfLevels() == 1 && !translation && !rotation && !scaling && !shearing) {
    t.GetLocalTransformation(0)->Write(output_name);
    if (verbose > 2) t.GetLocalTransformation(0)->Print(2);
  } else {
    t.Write(output_name);
    if (verbose > 2) t.Print(2);
  }
}
