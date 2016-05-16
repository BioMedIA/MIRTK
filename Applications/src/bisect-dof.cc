/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Christian Ledig
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

#include "mirtk/Matrix.h"
#include "mirtk/HomogeneousTransformation.h"

using namespace mirtk;


// ===========================================================================
// Help
// ===========================================================================

// ---------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofin> <dofout>" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This command bisects a rigid or affine transformation by calculating the" << endl;
  cout << "  matrix square root of the transformation matrix." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);
  UniquePtr<Transformation> dof(Transformation::New(POSARG(1)));
  HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(dof.get());
  if (!lin) {
    FatalError("Input transformation must be either Rigid, Similarity, or Affine");
  }
  lin->PutMatrix(lin->GetMatrix().Sqrt());
  lin->Write(POSARG(2));
  return 0;
}
