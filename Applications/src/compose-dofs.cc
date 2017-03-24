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
#include "mirtk/GenericImage.h"
#include "mirtk/Transformations.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <T1> <T2>... <T> [-target <image>]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes the composition T of the given input transformations such that" << endl;
  cout << endl;
  cout << "  .. math::" << endl;
  cout << endl;
  cout << "     T(x) = Tn o ... o T2 o T1(x)" << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -target <image>    Target image on which images will be resampled using" << endl;
  cout << "                     the composed transformation." << endl;
  cout << "  -[no]rotation      Whether to allow rotation    when composite transformation is affine. (default: on)" << endl;
  cout << "  -[no]translation   Whether to allow translation when composite transformation is affine. (default: on)" << endl;
  cout << "  -[no]scaling       Whether to allow scaling     when composite transformation is affine. (default: on)" << endl;
  cout << "  -[no]shearing      Whether to allow shearing    when composite transformation is affine. (default: on)" << endl;
  cout << "  -approximate       Approximate the composed transformation using a single FFD. (default: off)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Compose current composite transformation with any displacement field
void PushDisplacement(FluidFreeFormTransformation &t1,
                      Transformation              *t2,
                      ImageAttributes             &attr)
{
  if (!attr) {
    if (t1.NumberOfLevels() == 0) {
      MultiLevelTransformation *mffd;
      mffd = dynamic_cast<MultiLevelTransformation *>(t2);
      if (mffd && mffd->NumberOfLevels() > 0) {
        attr = mffd->GetLocalTransformation(-1)->Attributes();
      } else {
        // Note: As long as no local transformation was encountered before,
        //       no global transformation or single-level FFD has to be
        //       approximated by a displacement field. This error should
        //       thus never be encountered...
        cerr << "Internal error: Specify -target image to circumvent it" << endl;
        cerr << "                and please consider reporting the issue" << endl;
        exit(1);
      }
    } else {
      attr = t1.GetLocalTransformation(-1)->Attributes();
    }
  }
  GenericImage<double> disp(attr, 3);
  t2->Displacement(disp);
  t1.PushLocalTransformation(new LinearFreeFormTransformation3D(disp));
}

// -----------------------------------------------------------------------------
/// Compose current composite transformation with rigid/affine transformation
void PushTransformation(FluidFreeFormTransformation &t1,
                        HomogeneousTransformation   *t2,
                        ImageAttributes             &attr,
                        bool last = false)
{
  if (t2->IsIdentity()) return;
  if (t1.NumberOfLevels() == 0) {
    HomogeneousTransformation *global = t1.GetGlobalTransformation();
    global->PutMatrix(t2->GetMatrix() * global->GetMatrix());
  } else {
    if (last) {
      AffineTransformation *post = t1.GetAffineTransformation();
      post->PutMatrix(post->GetMatrix() * t2->GetMatrix());
    } else {
      PushDisplacement(t1, t2, attr);
    }
  }
}

// -----------------------------------------------------------------------------
/// Compose current composite transformation with single-level FFD
void PushTransformation(FluidFreeFormTransformation &t1,
                        FreeFormTransformation      *t2,
                        ImageAttributes             &attr,
                        bool = false)
{
  t1.PushLocalTransformation(t2);
}

// -----------------------------------------------------------------------------
/// Compose current composite transformation with multi-level FFD
void PushTransformation(FluidFreeFormTransformation      &t1,
                        MultiLevelFreeFormTransformation *t2,
                        ImageAttributes                  &attr,
                        bool = false)
{
  if (t2->NumberOfLevels() == 0) {
    PushTransformation(t1, t2->GetGlobalTransformation(), attr);
  } else if (t2->NumberOfLevels() == 1) {
    t2->MergeGlobalIntoLocalDisplacement();
    t1.PushLocalTransformation(t2->RemoveLocalTransformation(0));
  } else {
    PushDisplacement(t1, t2, attr);
  }
}

// -----------------------------------------------------------------------------
/// Compose current composite transformation with multi-level SV FFD
void PushTransformation(FluidFreeFormTransformation                &t1,
                        MultiLevelStationaryVelocityTransformation *t2,
                        ImageAttributes                            &attr,
                        bool = false)
{
  if (t2->NumberOfLevels() == 0) {
    PushTransformation(t1, t2->GetGlobalTransformation(), attr);
  } else if (t2->NumberOfLevels() == 1) {
    t2->MergeGlobalIntoLocalDisplacement();
    t1.PushLocalTransformation(t2->RemoveLocalTransformation(0));
  } else {
    PushDisplacement(t1, t2, attr);
  }
}

// -----------------------------------------------------------------------------
/// Compose current composite transformation with another composite transformation
void PushTransformation(FluidFreeFormTransformation &t1,
                        FluidFreeFormTransformation *t2,
                        ImageAttributes             &attr,
                        bool last = false)
{
  PushTransformation(t1, t2->GetGlobalTransformation(), attr);
  while (t2->NumberOfLevels() > 0) {
    PushTransformation(t1, t2->RemoveLocalTransformation(0), attr);
  }
  if (last) {
    AffineTransformation *post = t1.GetAffineTransformation();
    post->PutMatrix(post->GetMatrix() * t2->GetAffineTransformation()->GetMatrix());
  } else {
    PushTransformation(t1, t2->GetAffineTransformation(), attr);
  }
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
  double dx = 0., dy = 0., dz = 0.;

  for (ALL_OPTIONS) {
    if (OPTION("-target")) {
      InitializeIOLibrary();
      GreyImage target(ARGUMENT);
      attr = target.Attributes();
    }
    else if (OPTION("-bch")) {
      approximate = true;
      no_bch_terms = 4;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(no_bch_terms);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(lie_derivative);
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

  HomogeneousTransformation                  *aff;
  FreeFormTransformation                     *ffd;
  MultiLevelFreeFormTransformation           *mffd;
  MultiLevelStationaryVelocityTransformation *svmffd;
  FluidFreeFormTransformation                *fluid;

  // Compose affine transformations at the end
  if (verbose) {
    cout << "Compose affine transformations at end..." << endl;
  }
  AffineTransformation *post = t.GetAffineTransformation();
  for (int n = N; n >= 1; --n) {
    if (verbose) {
      cout << "Reading transformation " << n << ": " << POSARG(n) << endl;
    }
    UniquePtr<Transformation> dof(Transformation::New(POSARG(n)));
    if ((aff = dynamic_cast<HomogeneousTransformation *>(dof.get()))) {
      post->PutMatrix(post->GetMatrix() * aff->GetMatrix());
      --N; // remove transformation from composition chain
    } else {
      break;
    }
  }
  if (verbose) {
    cout << "Compose affine transformations at end... done" << endl;
  }

  // Compose remaining transformations from the start
  if (verbose) {
    cout << "\nCompose remaining transformations..." << endl;
  }
  for (int n = 1; n <= N; ++n) {

    // Read n-th transformation
    if (verbose) {
      cout << "Reading transformation " << n << " from " << POSARG(n) << endl;
    }
    UniquePtr<Transformation> dof(Transformation::New(POSARG(n)));
    aff    = dynamic_cast<HomogeneousTransformation                  *>(dof.get());
    ffd    = dynamic_cast<FreeFormTransformation                     *>(dof.get());
    mffd   = dynamic_cast<MultiLevelFreeFormTransformation           *>(dof.get());
    svmffd = dynamic_cast<MultiLevelStationaryVelocityTransformation *>(dof.get());
    fluid  = dynamic_cast<FluidFreeFormTransformation                *>(dof.get());

    // Compose current composite transformation with n-th transformation
    const bool last = (n == N);
    if      (aff   ) PushTransformation(t, aff,    attr, last);
    else if (mffd  ) PushTransformation(t, mffd,   attr, last);
    else if (svmffd) PushTransformation(t, svmffd, attr, last);
    else if (fluid ) PushTransformation(t, fluid,  attr, last);
    else if (ffd   ) {
      PushTransformation(t, ffd, attr, last);
      dof.release();
    }
    else {
      cerr << "Unsupported transformation file " << POSARG(n) << endl;
      cerr << "  Type name = " << dof->NameOfClass() << endl;
      exit(1);
    }

    // Transform target image attributes by affine transformation such that
    // attributes of consecutively approximated displacment fields overlap
    // with the thus far transformed image grid
    if (aff && attr) {
      attr.PutAffineMatrix(aff->GetMatrix(), true);
    }
  }
  if (verbose) {
    cout << "Compose remaining transformations... done\n" << endl;
  }

  // Approximate the composed transformation using a single free-form deformation
  if (approximate && t.NumberOfLevels() > 1) {
    double rms_error;

    // Use the most dense control point lattice in the local transformation stack
    ImageAttributes domain = attr;
    if (!domain) {
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
    const auto global = t.GetGlobalTransformation()->GetMatrix();
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
      auto * const first = dynamic_cast<const BSplineFreeFormTransformationSV *>(t.GetLocalTransformation(0));
      if (attr) {
        svffd.reset(new BSplineFreeFormTransformationSV());
        svffd->Initialize(domain, dx, dy, dz, first);
      } else {
        svffd.reset(new BSplineFreeFormTransformationSV(*first));
      }
      svffd->NumberOfBCHTerms(no_bch_terms);
      svffd->LieDerivative(lie_derivative);
      for (int i = 1; i < t.NumberOfLevels(); ++i) {
        svffd->CombineWith(t.GetLocalTransformation(i));
      }
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
