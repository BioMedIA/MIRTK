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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkImageIOConfig.h>

#include <mirtkTransformations.h>

using namespace mirtk;


// ===========================================================================
// Help
// ===========================================================================

// ---------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofout> <dofin>" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This command bisects a rigid or affine transformation by calculating the" << endl;
  cout << "  matrix square root of the transformation matrix." << endl;
  cout << endl;
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InitializeImageIOLibrary();

  REQUIRES_POSARGS(2);

  string        dofout;
  string	dofin;

  // Parse arguments
  dofout = POSARG(1);
  dofin  = POSARG(2);

  // Read transformation from file
  unique_ptr<Transformation> t(Transformation::New(dofin.c_str()));
  // Determine actual type of transformation
  RigidTransformation         *rigid      = NULL;
  AffineTransformation        *affine     = NULL;
  if (!((rigid      = dynamic_cast<RigidTransformation *>(t.get())) ||
        (affine     = dynamic_cast<AffineTransformation        *>(t.get())))) {
    cerr << EXECNAME << ": Cannot process transformation \"" << dofin << "\" of type " << t->NameOfClass() << endl;
    exit(1);
  }
  if (affine) {
  	Matrix bisectMat = affine->GetMatrix();
  	bisectMat = bisectMat.Sqrt();
  	affine->PutMatrix(bisectMat);
  	affine->Write(dofout.c_str());
  }
  else if (rigid) {
        Matrix bisectMat = rigid->GetMatrix();
	bisectMat = bisectMat.Sqrt();
	rigid->PutMatrix(bisectMat);
	rigid->Write(dofout.c_str());
  }
  return 0;
}
