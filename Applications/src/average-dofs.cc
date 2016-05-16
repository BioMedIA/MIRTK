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

#include "mirtk/Transformations.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/DisplacementToVelocityFieldBCH.h"
#include "mirtk/VelocityToDisplacementFieldSS.h"
#include "mirtk/ImageToInterpolationCoefficients.h"

using namespace mirtk;


// ===========================================================================
// Constants
// ===========================================================================

const double NaN     = numeric_limits<double>::quiet_NaN();
const double EPSILON = 0.001; // default -epsilon value

// ===========================================================================
// Help
// ===========================================================================

// ---------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofout> <dofin>...       [options]" << endl;
  cout << "       " << name << " <dofout> -dofnames <file> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  This command averages a number of input transformations." << endl;
  cout << "  It can be used to construct a brain image atlas and has been" << endl;
  cout << "  utilized for the construction of the spatio-temporal fetal/neonatal" << endl;
  cout << "  brain atlases available at http://brain-development.org/brain-atlases/." << endl;
  cout << endl;
  cout << "  The transformations which can be averaged by this program are listed below," << endl;
  cout << "  commonly expressed as \"sum of global and local transformations\", i.e.," << endl;
  cout << endl;
  cout << "  .. math::" << endl;
  cout << endl;
  cout << "     T(x) = T_{global}(x) + T_{local}(x)" << endl;
  cout << endl;
  cout << "  ======  ===============================================================================" << endl;
  cout << "  Linear  :math:`T(x) = Ax,         T_{global}(x) = A x, T_{local}(x) = 0`" << endl;
  cout << "  FFD     :math:`T(x) =  x + d(x),  T_{global}(x) = I x, T_{local}(x) = d(x)`" << endl;
  cout << "  MFFD    :math:`T(x) = Ax + d(x),  T_{global}(x) = A x, T_{local}(x) = d(x)`" << endl;
  cout << "  Fluid   :math:`T(x) =  x + d(Ax), T_{global}(x) = A x, T_{local}(x) = d(Ax) - (Ax - x)`" << endl;
  cout << "  ======  ===============================================================================" << endl;
  cout << endl;
  cout << "  where A is a linear 4x4 transformation matrix, and d(x) is a non-linear displacement." << endl;
  cout << "  Note that in case of the fluid transformation, d(x) is the total displacement" << endl;
  cout << "  instead of the local displacement only. By subtracting the displacement induced" << endl;
  cout << "  by the global transformation, the displacement corresponding to the local" << endl;
  cout << "  deformation only is obtained. Note that in case of MFFDs, the total displacement" << endl;
  cout << "  introduced by all local transformation levels is used to compute the average local" << endl;
  cout << "  transformation and thus multiple levels are allowed. Representing all transformations" << endl;
  cout << "  as sum of global and local transformations enables the computation of an average" << endl;
  cout << "  transformation for a given target domain from an arbitrary mix of input transformations." << endl;
  //
  //         Note: This computation is in fact implemented by the overridden virtual member
  //               function FluidFreeFormTransformation::LocalDisplacement and hence
  //               the virtual function Transformation::LocalDisplacement can be used
  //               in all cases to get the local component of the total displacement field.
  //
  cout << endl;
  cout << "  Accordingly, the resulting average transformation is set to be equal the sum of the" << endl;
  cout << "  Log-Euclidean means of the global transformation matrices, A, of the input transformations" << endl;
  cout << "  and the corresponding local transformations. The bi-invariant global mean is computed" << endl;
  cout << "  instead using the barycentric fixed point iteration when :option:`-bi-invariant` is given." << endl;
  cout << "  Before the computation of the local mean transformation, any dependency of the local" << endl;
  cout << "  displacement on the global component is removed by pre-multiplying the displacement" << endl;
  cout << "  vectors by the inverse matrix of the global transformation [1]_. Moreover, the (optional)" << endl;
  cout << "  local component of the final average transformation has to be taken into consideration" << endl;
  cout << "  as well. Therefore, all local input displacements are made relative to the same position" << endl;
  cout << "  after global transformation, i.e.," << endl;
  cout << endl;
  cout << "  .. math::" << endl;
  cout << endl;
  cout << "     d(x) = \\bar{A} \\circ A^{-1} \\circ T_{local}(x)" << endl;
  cout << endl;
  cout << "  The average of all transformation parameters is computed in the \"log-space\" by default." << endl;
  cout << "  Thus, in case of the local displacement fields, d(x), the corresponding stationary" << endl;
  cout << "  velocity fields, :math:`v(x) = log(d(x))`, are averaged instead. If the input transformation" << endl;
  cout << "  itself is parameterized by stationary velocities, these velocities are used directly" << endl;
  cout << "  when possible to avoid the redundant :math:`log(exp)` computation and its numerical error." << endl;
  cout << "  The resulting average velocity field is then either combined with a suitable" << endl;
  cout << "  velocity-based transformation model or exponentiated again to obtain the final" << endl;
  cout << "  average displacement field, :math:`\\bar{d}(x)`." << endl;
  cout << endl;
  cout << "  The output transformation will be an affine transformation if only the global transformation" << endl;
  cout << "  parameters were averaged, or a MFFD with one level corresponding to the average FFD otherwise." << endl;
  cout << "  The type of the output FFD depends on the type of the input transformations. If all transformations" << endl;
  cout << "  are either a linear transformation or a FFD with cubic B-spline or linear interpolation kernel" << endl;
  cout << "  (may also be contained in a one-level MFFD), the ouput FFD will contain the average parameters" << endl;
  cout << "  of the input FFDs. If the average is computed in the log-space, this only applies if all input" << endl;
  cout << "  FFDs are parameterized by stationary velocity fields. In all other cases, the output FFD will" << endl;
  cout << "  be either defined in the image domain of the input target image or the first encountered input" << endl;
  cout << "  FFD with a linear interpolation kernel and a control point spacing corresponding to the voxel size." << endl;
  cout << "  In other words, the output FFD is in general a dense (\"non-parametric\") displacement field." << endl;
  cout << endl;
  cout << "  .. [1] Rueckert et al., Automatic construction of 3-D statistical deformation" << endl;
  cout << "         models of the brain using nonrigid registration, IEEE TMI, 22(8), 1014–25 (2003)" << endl;
  cout << endl;
  cout << "Input options:" << endl;
  cout << "  -target <file>             Common reference image space of input transformations." << endl;
  cout << "                             If not specified and local transformations are to be averaged," << endl;
  cout << "                             the attributes of the first local transformation are used." << endl;
  cout << "  -dofnames <file>           Text file listing input transformations and associated values." << endl;
  cout << "  -dofdir <dir>              Directory used to make relative paths in :option:`-dofnames` text file absolute. (default: cwd)" << endl;
  cout << "  -prefix <string>           Prefix for transformation name entries in :option:`-dofnames` text file. (default: none)" << endl;
  cout << "  -suffix <string>           Suffix for transformation name entries in :option:`-dofnames` text file. (default: none)" << endl;
  cout << "  -gaussian <mean> <sigma>   Use Gaussian kernel weights. Requires :option:`-dofnames` to specify transformation values." << endl;
  cout << "                             By default, if :option:`-dofnames` is used or <sigma> is 0, the values specified in the" << endl;
  cout << "                             text file are directly used as kernel weights for the averaging instead of" << endl;
  cout << "                             using these values as arguments for the Gaussian kernel function." << endl;
  cout << "  -epsilon <value>           Weight threshold. (default: " << EPSILON << ")" << endl;
  cout << "  -add-identity              Assume additional identity transformation as part of input transformations. (default: none)" << endl;
  cout << "  -add-identity-with-weight <value>   Assume additional identity transformation with given value/weight. (default: none)" << endl;
  cout << "  -add-identity-for-dofname <name>    Assume identity transformation for named input transformation." << endl;
  cout << "                                      Note that if this option is used, the named input transformation file" << endl;
  cout << "                                      which is listed in the :option:`-dofnames` list does not need to exist." << endl;
  cout << "                                      (default: read input transformation from file)" << endl;
  cout << endl;
  cout << "Average transformation options:" << endl;
  cout << "  -[no]rotation              Average rotation    or assume none to be present. (default: off)" << endl;
  cout << "  -[no]translation           Average translation or assume none to be present. (default: off)" << endl;
  cout << "  -[no]scaling               Average scaling     or assume none to be present. (default: on)" << endl;
  cout << "  -[no]shearing              Average shearing    or assume none to be present. (default: on)" << endl;
  cout << "  -[no]deformation           Average deformation or assume none to be present. (default: on)" << endl;
  cout << "  -[no]log                   Whether to average local transformations in log-space." << endl;
  cout << "                             (default: yes, unless the :option:`-dofs` are averaged directly)" << endl;
  cout << "  -log-euclidean             Compute Log-Euclidean means. (default: off)" << endl;
  cout << "  -bi-invariant              Compute global bi-invariant mean, i.e., exponential barycenter. (default: on)" << endl;
  cout << "  -dofs                      Average the local transformation parameters directly." << endl;
  cout << "                             (default: corresponding (dense) displacement fields are averaged)" << endl;
  cout << "  -inverse-dofs              Average inverse input transformations and invert the average again. (default: off)" << endl;
  cout << endl;
  cout << "  -rigid                     Enables averaging of :option:`-translation` and :option:`-rotation` components," << endl;
  cout << "                             and disables averaging of :option:`-scaling`, :option:`-shearing`, and :option:`-deformation`." << endl;
  cout << "  -norigid                   Disables averaging of :option:`-translation` and :option:`-rotation` components." << endl;
  cout << "  -affine                    Enables averaging of :option:`-translation`, :option:`-rotation`," << endl;
  cout << "                             :option:`-scaling`, and :option:`-shearing`, components and disables output of" << endl;
  cout << "                             average :option:`-deformation`." << endl;
  cout << "  -noaffine                  Disables averaging of :option:`-rotation`, :option:`-translation`," << endl;
  cout << "                             :option:`-scaling`, and :option:`-shearing` components." << endl;
  cout << "  -all                       Enables averaging of :option:`-rotation`, :option:`-translation`," << endl;
  cout << "                             :option:`-scaling`, :option:`-shearing`, and :option:`-deformation`." << endl;
  cout << endl;
  cout << "  -linear                    Force average output transformation to use linear interpolation." << endl;
  cout << "  -cubic                     Force average output transformation to use cubic B-spline interpolation." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// ===========================================================================
// Auxiliaries
// ===========================================================================

// ---------------------------------------------------------------------------
bool isRelativePath(const string &path)
{
  return path.size() > 0 && path[0] != '/';
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char **argv)
{
  InitializeIOLibrary();

  REQUIRES_POSARGS(1);

  string        dofout;
  Array<string> dofin;
  Array<double> value;

  // Parse arguments
  dofout = POSARG(1);
  for (OPTIONAL_POSARGS) dofin.push_back(ARGUMENT);

  // Parse options
  const char *target_name    = NULL;
  const char *doflist        = NULL;
  const char *dofdir         = NULL;
  const char *prefix         = NULL;
  const char *suffix         = NULL;
  bool        translation    = false;
  bool        rotation       = false;
  bool        scaling        = true;
  bool        shearing       = true;
  bool        deformation    = true;
  bool        invert         = false;      // invert input transformations
  bool        invavg         = false;      // invert output transformation
  int         logspace       = -1;         // i.e., not set explicitly to "true"
  bool        avgdofs        = false;      // average DoFs directly if possible
  bool        approxglobal   = true;       // approximate vs. truncate average
                                           // global transformation if #DoFs < 12
  bool        bsplineffd     = false;      // output B-spline FFD when avgdofs == false
  bool        biinvariant    = true;       // bi-invariant mean if possible
  const char *identity_name  = "identity"; // name of (input) transformation which is to be replaced by identity
  double      identity_value = numeric_limits<double>::infinity(); // i.e., no implicit identity transformation added by default
  double      mean           = .0;
  double      sigma          = .0;
  double      epsilon        = EPSILON;
  int         frechet_iter   = 20;

  for (ALL_OPTIONS) {
    if      (OPTION("-target"))        target_name = ARGUMENT;
    else if (OPTION("-dofnames"))      doflist     = ARGUMENT;
    else if (OPTION("-dofdir"))        dofdir      = ARGUMENT;
    else if (OPTION("-prefix"))        prefix      = ARGUMENT;
    else if (OPTION("-suffix"))        suffix      = ARGUMENT;
    else if (OPTION("-invert"))        invert      = true;
    else if (OPTION("-inverse"))       invavg      = true;
    else if (OPTION("-translation"))   translation = true;
    else if (OPTION("-rotation"))      rotation    = true;
    else if (OPTION("-scaling"))       scaling     = true;
    else if (OPTION("-shearing"))      shearing    = true;
    else if (OPTION("-deformation"))   deformation = true;
    else if (OPTION("-notranslation")) translation = false;
    else if (OPTION("-norotation"))    rotation    = false;
    else if (OPTION("-noscaling"))     scaling     = false;
    else if (OPTION("-noshearing"))    shearing    = false;
    else if (OPTION("-nodeformation")) deformation = false;
    else if (OPTION("-bi-invariant"))  biinvariant = true;
    else if (OPTION("-log-euclidean")) biinvariant = false, logspace = true;
    else if (OPTION("-log"))           logspace    = true;
    else if (OPTION("-nolog"))         logspace    = false;
    else if (OPTION("-linear"))        bsplineffd  = false;
    else if (OPTION("-cubic"))         bsplineffd  = true;
    else if (OPTION("-dofs"))          avgdofs     = true;
    else if (OPTION("-inverse-dofs"))  avgdofs     = true, invert = invavg = true;
    else if (OPTION("-epsilon"))       epsilon     = atof(ARGUMENT);
    else if (OPTION("-gaussian")) { mean = atof(ARGUMENT); sigma = atof(ARGUMENT); }
    else if (OPTION("-rigid"))    { translation = rotation = true;  scaling = shearing = deformation = false; }
    else if (OPTION("-norigid"))  { translation = rotation = false; }
    else if (OPTION("-affine"))   { translation = rotation = scaling = shearing = true;  deformation = false; }
    else if (OPTION("-noaffine")) { translation = rotation = scaling = shearing = false; }
    else if (OPTION("-all"))      { translation = rotation = scaling = shearing = deformation = true; }
    else if (OPTION("-add-identity"))             { identity_value = NaN; }
    else if (OPTION("-add-identity-with-weight")) { identity_value = atof(ARGUMENT); }
    else if (OPTION("-add-identity-for-dofname")) { identity_name  = ARGUMENT; }
    else if (OPTION("-max-frechet-iterations")) frechet_iter = atoi(ARGUMENT);
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  if (sigma < .0) {
    cerr << EXECNAME << ": Gaussian sigma value must be positive" << endl;
    exit(1);
  }

  // Parse input text file
  if (doflist) {
    if (dofin.size() > 0) {
      cerr << EXECNAME << ": Text file specified for transformations but " << dofin.size() << " have already been read." << endl;
      exit(1);
    }
    ifstream in(doflist);
    if (!in) {
      cout << EXECNAME << ": Unable to open file " << doflist << endl;
      exit(1);
    }
    size_t lineno = 0;
    while (in) {
      lineno++;
      dofin.push_back("");
      value.push_back(NaN);
      in >> dofin.back() >> value.back();
      if (in.eof()) {
        dofin.pop_back();
        value.pop_back();
        break;
      }
      if (dofin.back().empty() || dofin.back()[0] == '#') {
        dofin.pop_back();
        value.pop_back();
      } else if (dofin.back() != identity_name) {
        if (prefix) dofin.back().insert(0, prefix);
        if (suffix) dofin.back() += suffix;
      }
      if (in.fail()) {
        cout << EXECNAME << ": Failed to parse line " << lineno << " of file " << doflist << endl;
        exit(1);
      }
    }
    in.close();
    if (dofin.empty()) {
      cerr << EXECNAME << ": No transformations listed in file " << doflist << endl;
      exit(1);
    }
    if (verbose) {
      cout << "Read " << dofin.size() << " transformations from file " << doflist << endl;
    }
  } else if (sigma != .0) {
      cerr << EXECNAME << ": Input transformations must be specified using -dofnames when using -gaussian kernel." << endl;
      exit(1);
  }

  // Number of input transformations, including optional additional implicit identity transformation
  const int N = static_cast<int>(IsInf(identity_value) ? dofin.size() : dofin.size()+1u);

  // Compute actual weights
  double *w = NULL; // actual weights
  double  W = .0;   // sum of all weights

  if (value.size() > 0) {
    // Add value for implicit identity transformation
    if (!IsInf(identity_value)) value.push_back(identity_value);
    w = new double[value.size()];
    for (size_t i = 0; i < value.size(); ++i) {
      if (sigma > .0) {
        w[i] = 1.0 / sqrt(2.0) / sigma;
        if (!IsNaN(value[i])) w[i] *= exp(-pow((mean - value[i])/sigma, 2.0)/2.0);
      } else {
        if (IsNaN(value[i])) w[i] = 1.0;
        else                 w[i] = value[i];
      }
      if (w[i] < epsilon) w[i] = .0;
      W += w[i];
    }
  }

  // If no weights specified, all transformations contribute equally
  if (W == .0) W = N;

  // Some verbose reporting
  if (verbose > 1) {
    int    nzero = 0;
    size_t width = 0;
    double minw  = (w ? w[0] : 1.0);
    double maxw  = (w ? w[0] : 1.0);
    double dw    = 1.0;
    for (size_t i = 0; i < value.size(); ++i) {
      if (w[i] < minw) minw = w[i];
      if (w[i] > maxw) maxw = w[i];
      if (w[i] == .0) { ++nzero; if (verbose < 2) continue; }
      if (dofin[i].length() > width) width = dofin[i].length();
      for (size_t j = i+1; j < value.size(); ++j) {
        double d = fabs(w[i] - w[j]);
        if (d > 1e-6 && d < dw) dw = d;
      }
    }
    if (value.size() > 0) {
      bool *reported = new bool[value.size()];
      for (size_t i = 0; i < value.size(); ++i) reported[i] = (verbose < 2 && w[i] == .0);
      for (double level = maxw; level >= minw; level -= dw) {
        for (size_t i = 0; i < value.size(); ++i) {
          if (reported[i] || w[i] < level) continue;
          printf("%3d. Weight of %-*s = %.5f\n", int(i+1), int(width), dofin[i].c_str(), w[i]);
          reported[i] = true;
        }
      }
      delete[] reported;
      printf("\n");
    }
    printf("Number of transformations with zero     weight: %d\n", nzero);
    printf("Number of transformations with non-zero weight: %d\n", N - nzero);
    printf("\n");
    printf("Total number of transformations:                %d\n",   N);
    printf("Total weight of transformations:                %.5f\n", W);
  }

  // Make input file paths absolute
  if (dofdir) {
    for (size_t i = 0; i < dofin.size(); ++i) {
      if (dofin[i] != identity_name && isRelativePath(dofin[i])) {
        dofin[i] = string(dofdir) + '/' + dofin[i];
      }
    }
  }

  // (Common) type and attributes of local input/output transformations
  string          mtype; // mulit-level transformation type
  string          type;  // common transformation type
  ImageAttributes attr;  // common FFD lattice attributes

  // Average displacement fields sampled in specified target image domain
  if (!avgdofs && target_name) {
    UniquePtr<BaseImage> target(BaseImage::New(target_name));
    attr = target->GetImageAttributes();
  // Otherwise, determine common type of local input transformations
  // if a mix of transformations is given, use attributes of first FFD
  // to define the image domain of the average displacement field
  } else {
    if (verbose) cout << "Checking type of input transformations...", cout.flush();
    for (size_t i = 0; i < dofin.size(); ++i) {
      if (dofin[i] == identity_name) continue;
      UniquePtr<Transformation> t(Transformation::New(dofin[i].c_str()));
      Transformation           *p    = t.get();
      MultiLevelTransformation *mffd = dynamic_cast<MultiLevelTransformation *>(t.get());
      if (mffd) {
        if      (mffd->NumberOfLevels() == 0) p = mffd->GetGlobalTransformation();
        else if (mffd->NumberOfLevels() == 1) p = mffd->GetLocalTransformation(0);
        else {
          attr = mffd->GetLocalTransformation(-1)->Attributes();
          type.clear();
          break;
        }
        if      (mtype.empty())                mtype = mffd->NameOfClass();
        else if (mtype != mffd->NameOfClass()) mtype.clear();
      }
      FreeFormTransformation *ffd = dynamic_cast<FreeFormTransformation *>(p);
      if (ffd) {
        if (type.empty()) {
          type = ffd->NameOfClass();
          attr = ffd->Attributes();
        } else if (type != ffd->NameOfClass() || ffd->Attributes() != attr) {
          mtype.clear();
          type .clear();
          break;
        }
      }
    }
    if (verbose) cout << " done" << endl;
  }

  // Cases in which local transformation parameters are averaged directly and stored
  // using the very same transformation model as given as input which does not require
  // a conversion of the input transformations to dense displacement fields
  if (avgdofs && deformation) {
    if (attr.NumberOfLatticePoints() > 0 && // input contains FFDs to be averaged
        type != "BSplineFreeFormTransformationSV"   &&
        type != "LinearFreeFormTransformation3D"    &&
        type != "BSplineFreeFormTransformation3D") {
      if (verbose) cout << endl;
      cerr << EXECNAME << ": Cannot average local input transformation parameters directly (-dofs option)." << endl;
      cerr << endl;
      cerr << "  Transformations must therefore either all be linear transformations or transformations with local" << endl;
      cerr << "  component of type \"Linear FFD\", \"B-spline FFD\", or \"SV B-spline FFD\" (see dofprint output)." << endl;
      exit(1);
    }
    if (logspace == 1) {
      if (type != "BSplineFreeFormTransformationSV") {
        if (verbose) cout << endl;
        cerr << EXECNAME << ": Combining options -log and -dofs only possible for SV B-spline FFDs! Use only one of them." << endl;
        exit(1);
      }
    } else if (logspace == 0) {
      if (type != "LinearFreeFormTransformation" && type != "BSplineFreeFormTransformation3D") {
        if (verbose) cout << endl;
        cerr << EXECNAME << ": Combining options -nolog and -dofs only possible for Linear/B-spline FFDs! Use only one of them." << endl;
        exit(1);
      }
    }
  }

  if (verbose) {
    cout << "Average deformation:   " << ToString(deformation) << endl;
    cout << "Average DoFs directly: " << ToString(avgdofs) << endl;
    cout << "Average in log-space:  " << ToString(logspace) << endl;
  }

  // Initialize intermediate data structures
  Matrix              *A = NULL;
  GenericImage<double> d, avgD;

  if (translation || rotation || scaling || shearing) {
    A = new Matrix[N];
    for (int i = 0; i < N; ++i) {
      A[i].Initialize(4, 4);
      A[i].Ident();
    }
  }
  if (deformation && attr.NumberOfLatticePoints() > 0) {
    d   .Initialize(attr, 3);
    avgD.Initialize(attr, 3);
  }

  // Read parameters of input transformations
  for (size_t i = 0; i < dofin.size(); ++i) {
    if (dofin[i] == identity_name) continue;
    // Read transformation from file
    UniquePtr<Transformation> t(Transformation::New(dofin[i].c_str()));
    // Determine actual type of transformation
    HomogeneousTransformation   *global     = NULL;
    RigidTransformation         *rigid      = NULL;
    SimilarityTransformation    *similarity = NULL;
    AffineTransformation        *affine     = NULL;
    FreeFormTransformation3D    *ffd        = NULL;
    MultiLevelTransformation    *mffd       = NULL;
    FluidFreeFormTransformation *fluid      = NULL;
    if (!((fluid      = dynamic_cast<FluidFreeFormTransformation *>(t.get())) ||
          (mffd       = dynamic_cast<MultiLevelTransformation    *>(t.get())) ||
          (ffd        = dynamic_cast<FreeFormTransformation3D    *>(t.get())) ||
          (affine     = dynamic_cast<AffineTransformation        *>(t.get())) ||
          (similarity = dynamic_cast<SimilarityTransformation    *>(t.get())) ||
          (rigid      = dynamic_cast<RigidTransformation         *>(t.get())))) {
      cerr << EXECNAME << ": Cannot process transformation \"" << dofin[i] << "\" of type " << t->NameOfClass() << endl;
      exit(1);
    }
    // Set base class pointers for common interface access,
    // e.g., rigid components can be modified using the common
    // RigidTransformation interface even if the actual transformation
    // is of type SimilarityTransformation, AffineTransformation,
    // or even a subclass of MultiLevelTransformation.
    if (fluid) mffd = fluid;
    if (mffd) {
      affine = mffd->GetGlobalTransformation();
      if (mffd->NumberOfLevels() > 0) {
        FreeFormTransformation *affd = mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);
        ffd = dynamic_cast<FreeFormTransformation3D *>(affd);
        if (!ffd) {
          cerr << EXECNAME << ": Cannot process MFFD \"" << dofin[i] << "\" with level of type " << affd->NameOfClass() << endl;
          exit(1);
        }
      }
    }
    if (affine    ) similarity = affine;
    if (similarity) rigid      = similarity;
    if (rigid     ) global     = rigid;
    // Get local transformation parameters
    if (deformation && ffd) {
      if (invert) {
        cerr << EXECNAME << ": -invert option only supported for rigid/affine transformations" << endl;
        exit(1);
      }
      // Get inverse of the linear transformation matrix
      Matrix invA(4, 4);
      if (global) {
        invA = global->GetMatrix();
        invA.Invert();
      } else invA.Ident();
      // Average parameters at control points directly when possible
      double x, y, z;
      if (avgdofs) {
        for (int k = 0; k < attr._z; ++k) {
          for (int j = 0; j < attr._y; ++j) {
            for (int i = 0; i < attr._x; ++i) {
              // Get control point parameters
              ffd->Get(i, j, k, x, y, z);
              // Remove dependency on global transformation
              // Note: Applies also when the FFD parameters are stationary velocities, see below.
              d(i, j, k, 0) = invA(0, 0) * x + invA(0, 1) * y + invA(0, 2) * z;
              d(i, j, k, 1) = invA(1, 0) * x + invA(1, 1) * y + invA(1, 2) * z;
              d(i, j, k, 2) = invA(2, 0) * x + invA(2, 1) * y + invA(2, 2) * z;
            }
          }
        }
      // Otherwise,...
      } else {
        // Get local displacement field
        Matrix i2w = attr.GetImageToWorldMatrix();
        for (int k = 0; k < attr._z; ++k) {
          for (int j = 0; j < attr._y; ++j) {
            for (int i = 0; i < attr._x; ++i) {
              // Convert voxel indices to world coordinates
              x = i2w(0, 0) * i + i2w(0, 1) * j + i2w(0, 2) * k + i2w(0, 3);
              y = i2w(1, 0) * i + i2w(1, 1) * j + i2w(1, 2) * k + i2w(1, 3);
              z = i2w(2, 0) * i + i2w(2, 1) * j + i2w(2, 2) * k + i2w(2, 3);
              // Evaluate (total) local voxel displacement
              t->LocalDisplacement(x, y, z);
              // Remove dependency on global transformation
              d(i, j, k, 0) = invA(0, 0) * x + invA(0, 1) * y + invA(0, 2) * z;
              d(i, j, k, 1) = invA(1, 0) * x + invA(1, 1) * y + invA(1, 2) * z;
              d(i, j, k, 2) = invA(2, 0) * x + invA(2, 1) * y + invA(2, 2) * z;
            }
          }
        }
        // Compute stationary velocity field
        if (logspace) {
          DisplacementToVelocityFieldBCH<double> dtov;
          dtov.Input (&d);
          dtov.Output(&d);
          dtov.Run();
          // Smooth velocities if only few transformations are being averaged,
          // otherwise rely on the average velocity field to be sufficiently smooth
          if (N < 5) {
            GaussianBlurring<double> blur(max(attr._dx, max(attr._dy, attr._dz)));
            blur.Input (&d);
            blur.Output(&d);
            blur.Run();
          }
        }
      }
      // Add to sum of local transformations
      if (w) d *= w[i];
      avgD += d;
    }
    // Get global transformation parameters
    // (**after** local parameters as the following modifies the transformation)
    if (A && global) {
      if (invert) global->Invert();
      if (!translation) {
        rigid->PutTranslationX(0);
        rigid->PutTranslationY(0);
        rigid->PutTranslationZ(0);
      }
      if (!rotation) {
        rigid->PutRotationX(0);
        rigid->PutRotationY(0);
        rigid->PutRotationZ(0);
      }
      if (!scaling) {
        if (affine) {
          affine->PutScaleX(100);
          affine->PutScaleY(100);
          affine->PutScaleZ(100);
        } else if (similarity) {
          similarity->PutScale(100);
        }
      }
      if (!shearing && affine) {
        affine->PutShearXY(0);
        affine->PutShearXZ(0);
        affine->PutShearYZ(0);
      }
      A[i] = global->GetMatrix();
    }
  }

  // Compute average global transformation
  AffineTransformation globalAvg;
  globalAvg.AllowTranslations(translation);
  globalAvg.AllowRotations(rotation);
  globalAvg.AllowScaling(scaling);
  globalAvg.AllowShearing(shearing);

  if (A) {
    if (verbose) cout << "\nComputing average global transformation..." << endl;
    Matrix avgA;
    if (biinvariant) avgA = BiInvariantMean (N, A, w, frechet_iter);
    else             avgA = LogEuclideanMean(N, A, w);
    if (invavg) avgA.Invert();
    if (approxglobal) globalAvg.ApproximateAsNew(avgA);
    else              globalAvg.PutMatrix(avgA);
    if (verbose) {
      Indent indent(1);
      if (verbose > 1) {
        cout << "\n";
        avgA.Print(indent);
      }
      const Matrix avgAInv = avgA.Inverse();
      Matrix residual(4, 4);
      for (int i = 0; i < N; ++i) {
        if (w[i] == .0) continue;
        residual += (avgAInv * A[i]).Log() * w[i];
      }
      residual /= W;
      cout << "\n" << indent << "Deviation from barycentric equation   = " << residual.Norm() << endl;
      residual.Ident();
      residual -= avgA * globalAvg.GetMatrix().Inverse();
      cout << indent << "Norm of (A * inv(global._matrix) - I) = " << residual.Norm() << "\n" << endl;
      globalAvg.Print(indent + 1);
    }
    if (verbose) cout << "\nComputing average global transformation... done" << endl;
  }

  // If global average only...
  if (avgD.IsEmpty()) {

    // Write average affine transformation
    globalAvg.Write(dofout.c_str());

  // Otherwise,...
  } else {

    if (invavg) {
      cerr << EXECNAME << ": -inverse[-dofs] option only supported for global transformations" << endl;
      exit(1);
    }

    // Divide sum of local transformations by total weight to get average
    avgD /= W;
    // Convert velocities back to displacements, unless the output transformation
    // is parameterized by a stationary velocity field
    if (logspace && !avgdofs) {
      VelocityToDisplacementFieldSS<double> vtod;
      vtod.Input (&avgD);
      vtod.Output(&avgD);
      vtod.Run();
    }
    // Add dependency on average global transformation such that
    // avgT = avgA (x + avgD) = avgA x + avgA avgD = avgT_global + avgT_local
    //
    // Note that also in case of avgD being a stationary velocity field (avgdofs == true),
    // we can pre-multiply it the same way with the global transformation matrix as in
    // case of a displacement field. This can be illustrated by looking at the formula
    // of the family of explicit Runge-Kutta methods for the exponentiation step which
    // consists of a sum of the velocities at different time steps, i.e.,
    // avgT(x) = avgA x + avgA sum_i b_i h v(x_i) = avgA x + sum_i b_i h (avgA v(x_i)).
    double x, y, z;
    const Matrix &avgA = globalAvg.GetMatrix();
    for (int k = 0; k < attr._z; ++k) {
      for (int j = 0; j < attr._y; ++j) {
        for (int i = 0; i < attr._x; ++i) {
          x = avgD(i, j, k, 0), y = avgD(i, j, k, 1), z = avgD(i, j, k, 2);
          avgD(i, j, k, 0) = avgA(0, 0) * x + avgA(0, 1) * y + avgA(0, 2) * z;
          avgD(i, j, k, 1) = avgA(1, 0) * x + avgA(1, 1) * y + avgA(1, 2) * z;
          avgD(i, j, k, 2) = avgA(2, 0) * x + avgA(2, 1) * y + avgA(2, 2) * z;
        }
      }
    }
    // Construct FFD from average deformation
    FreeFormTransformation *ffd = NULL;
    if (avgdofs) {
      if (type == "LinearFreeFormTransformation3D") {
        ffd = new LinearFreeFormTransformation3D(avgD, true);
      } else if (type == "BSplineFreeFormTransformation3D") {
        ffd = new BSplineFreeFormTransformation3D(avgD, true);
      } else if (type == "BSplineFreeFormTransformationSV") {
        ffd = new BSplineFreeFormTransformationSV(avgD, true);
      }
    }
    if (!ffd) {
      if (bsplineffd) {
        ConvertToCubicBSplineCoefficients(avgD, 0);
        ConvertToCubicBSplineCoefficients(avgD, 1);
        ConvertToCubicBSplineCoefficients(avgD, 2);
        ffd = new BSplineFreeFormTransformation3D(avgD, true);
      } else {
        ffd = new LinearFreeFormTransformation3D(avgD, true);
      }
    }
    // Write average non-linear transformation
    MultiLevelTransformation *mffd = NULL;
    if (mtype == "MultiLevelStationaryVelocityTransformation") {
      mffd = new MultiLevelStationaryVelocityTransformation(globalAvg);
    } else {
      mffd = new MultiLevelFreeFormTransformation(globalAvg);
    }
    mffd->PushLocalTransformation(ffd);
    mffd->Write(dofout.c_str());
    delete mffd;
  }

  // Clean up
  delete[] w;
  delete[] A;

  return 0;
}
