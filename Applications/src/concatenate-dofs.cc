/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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
  cout << "  -target <image>   Target image on which images will be resampled using" << endl;
  cout << "                    the composed transformation." << endl;
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
                        bool                             last = false)
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
                        bool                             last = false)
{
  PushTransformation(t1, t2->GetGlobalTransformation(), attr);
  for (int l = 0; l < t2->NumberOfLevels(); ++l) {
    PushTransformation(t1, t2->RemoveLocalTransformation(l), attr);
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
  int N = NUM_POSARGS;

  ImageAttributes attr;
  FluidFreeFormTransformation t;

  for (ALL_OPTIONS) {
    if (OPTION("-target")) {
      InitializeIOLibrary();
      GreyImage target(ARGUMENT);
      attr = target.Attributes();
    }
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
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
  for (int n = N - 1; n >= 1; --n) {
    if (verbose) {
      cout << "Reading transformation " << n << ": " << POSARG(n) << endl;
    }
    unique_ptr<Transformation> dof(Transformation::New(POSARG(n)));
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
  for (int n = 1; n < N; ++n) {

    // Read n-th transformation
    if (verbose) {
      cout << "Reading transformation " << n << " from " << POSARG(n) << endl;
    }
    unique_ptr<Transformation> dof(Transformation::New(POSARG(n)));
    aff    = dynamic_cast<HomogeneousTransformation                  *>(dof.get());
    ffd    = dynamic_cast<FreeFormTransformation                     *>(dof.get());
    mffd   = dynamic_cast<MultiLevelFreeFormTransformation           *>(dof.get());
    svmffd = dynamic_cast<MultiLevelStationaryVelocityTransformation *>(dof.get());
    fluid  = dynamic_cast<FluidFreeFormTransformation                *>(dof.get());

    // Compose current composite transformation with n-th transformation
    const bool last = (n == N - 1);
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
    cout << "Compose remaining transformations... done" << endl;
  }

  // Write composite transformation
  const char *output_name = POSARG(NUM_POSARGS);
  if (verbose) cout << "Writing composite transformation to " << output_name << endl;
  t.Write(output_name);
  if (verbose) t.Print(2);
}
