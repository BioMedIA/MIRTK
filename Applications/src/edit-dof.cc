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

#include "mirtk/FreeFormTransformation.h"
#include "mirtk/MultiLevelTransformation.h"

#include "mirtk/BSplineFreeFormTransformationSV.h"
#include "mirtk/BSplineFreeFormTransformationTD.h"

#include <locale>
#include <algorithm>


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofin> <dofout> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This tool provides a generic way of modifying the parameters of a" << endl;
  cout << "  transformation (linear as well as non-linear). It uses the generic" << endl;
  cout << "  interface member function Transformation::Set whose arguments" << endl;
  cout << "  are the name of the parameter and the value to be set as string." << endl;
  cout << endl;
  cout << "  To allow this tool to be used like any other command-line tool," << endl;
  cout << "  it performs a generic conversion from option names to parameter names." << endl;
  cout << "  While option names are all lowercase and start with one or two hyphens" << endl;
  cout << "  and use hyphens also as word separators, parameter names start with a" << endl;
  cout << "  captial letter and use spaces as word separators but also hyphens for" << endl;
  cout << "  some compound words. See examples below." << endl;
  PrintStandardOptions(cout);
  cout << endl;
  cout << "Examples:" << endl;
  cout << "  " << name << " svffd.dof.gz svffd.dof.gz -number-of-integration-steps 128" << endl;
  cout << endl;
  cout << "      The option \"-number-of-integration-steps\" used here corresponds to" << endl;
  cout << "      the parameter \"Number of integration steps\" of the SV FFD transformation." << endl;
  cout << "      The new value for the parameter will be 128." << endl;
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Predicate used by ci_find
/// \see http://stackoverflow.com/questions/3152241/case-insensitive-stdstring-find
template <typename charT>
struct ci_equal {
  ci_equal(const std::locale& loc) : _loc(loc) {}

  bool operator()(charT ch1, charT ch2) {
    return toupper(ch1, _loc) == toupper(ch2, _loc);
  }

private:
  const std::locale& _loc;
};

// -----------------------------------------------------------------------------
/// Find substring (case insensitive)
/// \see http://stackoverflow.com/questions/3152241/case-insensitive-stdstring-find
template <typename T>
size_t ci_find(const T& str1, const T& str2, const std::locale& loc = std::locale())
{
  typename T::const_iterator it;
  it = std::search(str1.begin(), str1.end(),
                   str2.begin(), str2.end(),
                   ci_equal<typename T::value_type>(loc));
  if (it != str1.end()) return it - str1.begin();
  return string::npos;
}

// -----------------------------------------------------------------------------
/// Convert command-line option name to parameter name where the first word
/// starts with a captial letter and separating hyphens are replaced by spaces.
/// This function explicitly has to know about hyphenated compound words used
/// for parameters names of transformations such as "cross-sectional".
string ToParameterName(const char *opt)
{
  size_t pos;
  // skip leading hyphens
  if (opt[0] != '-') {
    cerr << "Invalid argument: " << opt << "! Expected an option starting with a leading hyphen." << endl;
    exit(1);
  }
  ++opt;
  if (opt[0] == '-') ++opt; // optional second hyphen
  // allocate output string
  string name;
  name.reserve(strlen(opt)-1);
  // add starting letter and capitalize it if necessary
  if      (isupper(opt[0])) name += opt[0];
  else if (islower(opt[0])) name += 'A' + (opt[0] - 'a');
  else {
    cerr << "Invalid option: -[-]" << opt << "! Expected a letter after the leading hyphen(s)." << endl;
    exit(1);
  }
  ++opt;
  // copy remaining option characters, replacing hyphens by spaces
  for (; opt[0]; ++opt) {
    switch (opt[0]) {
      case '-': name += ' '; break;
      default:  name += opt[0];
    }
  }
  // hyphenate compound words again
  pos = ci_find(name, string("cross sectional"));
  if (pos != string::npos) name[pos+5] = '-';
  // capitalize acronyms (again)
  pos = ci_find(name, string("bch"));
  if (pos != string::npos) { name[pos] = 'B'; name[pos+1] = 'C'; name[pos+2] = 'H'; }
  // end abbreviations with a period
  pos = ci_find(name, string("no "));
  if (pos != string::npos) name.insert(pos+2, ".");
  return name;
}

// -----------------------------------------------------------------------------
FreeFormTransformation *ConvertToDisplacements(FreeFormTransformation *ffd)
{
  cerr << "ConvertToDisplacements: Not implemented" << endl;
  exit(1);
  return ffd;
}

// -----------------------------------------------------------------------------
FreeFormTransformation *CastToDisplacements(FreeFormTransformation *ffd)
{
  BSplineFreeFormTransformationSV *svffd = dynamic_cast<BSplineFreeFormTransformationSV *>(ffd);
  if (svffd) return new BSplineFreeFormTransformation3D(*svffd);
  BSplineFreeFormTransformationTD *tdffd = dynamic_cast<BSplineFreeFormTransformationTD *>(ffd);
  if (tdffd) return new BSplineFreeFormTransformation4D(*tdffd);
  return ffd;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(2);

  // Parse command-line options
  ParameterList params;
  double        scale   = 1.0;
  bool          v2d     = false;
  bool          convert = true;

  for (ALL_OPTIONS) {
    if      (OPTION("-disp"))    v2d     = true;
    else if (OPTION("-cast"))    convert = false;
    else if (OPTION("-convert")) convert = true;
    else if (OPTION("-scale"))   PARSE_ARGUMENT(scale);
    else if (STANDARD_OPTION)    PARSE_STANDARD_OPTION();
    else if (IS_OPTION) Insert(params, ToParameterName(OPTNAME), ARGUMENT);
    else {
      cerr << "Invalid argument: " << argv[ARGIDX] << endl;
      exit(1);
    }
  }

  // Read transformation
  Transformation *dof = Transformation::New(POSARG(1));

  // Modify (non-DoF) parameters
  int exit_code = 0;
  for (ParameterConstIterator it = params.begin(); it != params.end(); ++it) {
    if (!dof->Set(it->first.c_str(), it->second.c_str())) {
      cerr << "Transformation of type " << dof->NameOfClass() << " either has no parameter named \"" << it->first << "\"," << endl;
      cerr << "does not allow this parameter to be modified through the generic interface, or does not implement this interface." << endl;
      exit_code = 1;
    }
  }

  // Scale DoFs
  if (scale != 1.0) {
    for (int i = 0; i < dof->NumberOfDOFs(); ++i) {
      dof->Put(i, scale * dof->Get(i));
    }
  }

  // Convert velocities to displacements
  if (v2d) {
    MultiLevelTransformation *mffd = dynamic_cast<MultiLevelTransformation *>(dof);
    FreeFormTransformation   *affd = dynamic_cast<FreeFormTransformation   *>(dof);
    FreeFormTransformation   *ffd;

    if (mffd) {
      for (int i = 0; i < mffd->NumberOfLevels(); ++i) {
        if (convert) {
          ffd = ConvertToDisplacements(mffd->GetLocalTransformation(i));
        } else {
          ffd = CastToDisplacements(mffd->GetLocalTransformation(i));
        }
        if (ffd != mffd->GetLocalTransformation(i)) {
          delete mffd->GetLocalTransformation(i);
          mffd->PutLocalTransformation(ffd, i);
        }
      }
    } else if (affd) {
      if (convert) ffd = ConvertToDisplacements(affd);
      else         ffd = CastToDisplacements   (affd);
      if (ffd != affd) delete affd;
      dof = ffd;
    }
  }

  // Write modified transformation
  if (exit_code == 0) {
    if (verbose) dof->Print();
    dof->Write(POSARG(2));
  }

  // Clean up
  delete dof;

  return exit_code;
}
