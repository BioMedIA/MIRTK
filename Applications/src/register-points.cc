/*
 * Medical Image Registration ToolKit (MIRTK)
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

#include "mirtk/PointCorrespondence.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/RegisteredPointSet.h"

#include "mirtk/Transformation.h"
#include "mirtk/TransformationModel.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/MultiLevelTransformation.h"

using namespace mirtk;


// TODO: Support also non-rigid transformations. In case of a FFD, support especially also a MFFD
//       which implements the multi-level scatter data approximation with cubic B-splines.


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Register point sets or surfaces by iteratively approximating the residual target registration error\n";
  cout << "  given current point correspondences and transformation estimate. By default, correspondences are defined\n";
  cout << "  by closest points. The default settings implement the iterative closest points (ICP) algorithm. When the\n";
  cout << "  input point sets are surface meshes, :option:`-closest-surface-points` can be used as correspondences\n";
  cout << "  instead of discrete points in order to increase the accuracy of the surface registration at the expense\n";
  cout << "  of a slightly more costly point correspondence update at each iteration.\n";
  cout << "\n";
  cout << "Required options:\n";
  cout << "  -t, -target <path>...\n";
  cout << "    One or more target point sets or surfaces. If multiple targets are specified,\n";
  cout << "    the number of target point sets must match the number of :option:`-source` point sets.\n";
  cout << "    In this case, the respective pairs of target and source point sets are registered\n";
  cout << "    to one another using the same target to source transformation. This option can\n";
  cout << "    be given multiple times to append target point sets.\n";
  cout << "  -s, -source <path>...\n";
  cout << "    One or more source point sets or surfaces. If a single :option:`-target` point set is given,\n";
  cout << "    all source point sets are registered to this common target point set. Otherwise,\n";
  cout << "    the number of target and source point sets must be the same. This option can be\n";
  cout << "    can be given multiple times to append source point sets.\n";
  cout << "  -o, -dofout <path>\n";
  cout << "    Output transformation file (.dof, .dof.gz).\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -i, -dofin <path>|Id|identity\n";
  cout << "    File path of input transformation from which initial -model parameters are derived.\n";
  cout << "  -m, -model Rigid|Similarity|Affine\n";
  cout << "    Transformation model to use. (default: Rigid)\n";
  cout << "  -n, -iterations <n>\n";
  cout << "    Number of iterations. (default: 100)\n";
  cout << "  -c, -cor, -corr, -correspondence <name>\n";
  cout << "    Name of correspondence type, e.g., 'closest point', 'closest cell'. (default: closest point)\n";
  cout << "  -cp, -closest-point\n";
  cout << "    Alias for :option:`-correspondence 'closest point'`.\n";
  cout << "  -csp, -closest-surface-point, -closest-cell\n";
  cout << "    Alias for :option:`-correspondence 'closest surface point'`.\n";
  cout << "  -p, -par, -corpar, -corrpar <name> <value>\n";
  cout << "    Set parameter of chosen correspondence type.\n";
  cout << "  -f, -feature <name> [<weight>]\n";
  cout << "    Point or cell features based on which correspondences are determined.\n";
  cout << "    By default, the spatial 3D coordinates of the input points are used.\n";
  cout << "  -[no]symmetric [on|off]\n";
  cout << "    Approximate symmetric registration error. When this flag is set, the distance\n";
  cout << "    from every source point to its corresponding target point is considered in addition\n";
  cout << "    to the distance of each target point to its corresponding source point. (default: off)\n";
  cout << "  -[no]inverse [on|off]\n";
  cout << "    By default, the output transformation is applied to the target points.\n";
  cout << "    When this flag is set, the output transformation is applied to the source\n";
  cout << "    points instead. (default: off)\n";
  PrintStandardOptions(cout);
  cout << "\n";
  cout.flush();
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
inline int TargetIndex(int m, int n, int i)
{
  return n == 1 ? 0 : i;
}

// -----------------------------------------------------------------------------
void Update(Array<UniquePtr<PointCorrespondence>> &cmaps)
{
  for (auto &cmap : cmaps) {
    cmap->Update();
  }
}

// -----------------------------------------------------------------------------
void Update(Array<RegisteredPointSet> &psets)
{
  for (auto &pset : psets) {
    pset.Update();
  }
}

// -----------------------------------------------------------------------------
void Update(Array<RegisteredPointSet> &targets,
            Array<RegisteredPointSet> &sources,
            Array<UniquePtr<PointCorrespondence>> &cmaps)
{
  Update(targets);
  Update(sources);
  Update(cmaps);
}

// -----------------------------------------------------------------------------
double EvaluateRMSError(const Array<RegisteredPointSet> &targets,
                        const Array<RegisteredPointSet> &sources,
                        const Array<UniquePtr<PointCorrespondence>> &cmaps,
                        bool symmetric = false)
{
  Point p, q;

  const int m = static_cast<int>(sources.size());
  const int n = static_cast<int>(targets.size());

  double error = 0;
  int count = 0;

  for (int i = 0; i < m; ++i) {
    const int j = TargetIndex(m, n, i);
    auto &cmap = cmaps[i];
    auto &source = sources[i];
    auto &target = targets[j];
    for (int t = 0; t < target.NumberOfPoints(); ++t) {
      target.GetPoint(t, p);
      if (cmap->GetSourcePoint(t, q)) {
        error += pow(q._x - p._x, 2) + pow(q._y - p._y, 2) + pow(q._z - p._z, 2);
        ++count;
      }
    }
    if (symmetric) {
      for (int s = 0; s < source.NumberOfPoints(); ++s) {
        source.GetPoint(s, q);
        if (cmap->GetTargetPoint(s, p)) {
          error += pow(q._x - p._x, 2) + pow(q._y - p._y, 2) + pow(q._z - p._z, 2);
          ++count;
        }
      }
    }
  }

  if (count == 0) FatalError("No corresponding points found");
  return sqrt(error / count);
}

// -----------------------------------------------------------------------------
void Fit(Transformation *dof,
         const Array<RegisteredPointSet> &targets,
         const Array<RegisteredPointSet> &sources,
         const Array<UniquePtr<PointCorrespondence>> &cmaps,
         bool symmetric = false, bool inverse = false)
{
  const int m = static_cast<int>(sources.size());
  const int n = static_cast<int>(targets.size());

  int no = 0;
  for (int i = 0; i < m; ++i) {
    const int j = TargetIndex(m, n, i);
    no += targets[j].NumberOfPoints();
    if (symmetric) no += sources[i].NumberOfPoints();
  }

  Array<double> x(no), y(no), z(no), dx(no), dy(no), dz(no);
  Point p, q;

  int k = 0;
  for (int i = 0; i < m; ++i) {
    const int j = TargetIndex(m, n, i);
    auto &cmap = cmaps[i];
    auto &source = sources[i];
    auto &target = targets[j];
    for (int t = 0; t < target.NumberOfPoints(); ++t) {
      target.GetInputPoint(t, p);
      if (cmap->GetInputSourcePoint(t, q)) {
        if (inverse) swap(p, q);
        x[k] = p._x;
        y[k] = p._y;
        z[k] = p._z;
        dx[k] = q._x - p._x;
        dy[k] = q._y - p._y;
        dz[k] = q._z - p._z;
        ++k;
      }
    }
    if (symmetric) {
      for (int s = 0; s < source.NumberOfPoints(); ++s) {
        source.GetInputPoint(s, q);
        if (cmap->GetInputTargetPoint(s, p)) {
          if (inverse) swap(p, q);
          x[k] = p._x;
          y[k] = p._y;
          z[k] = p._z;
          dx[k] = q._x - p._x;
          dy[k] = q._y - p._y;
          dz[k] = q._z - p._z;
          ++k;
        }
      }
    }
  }
  dof->ApproximateAsNew(x.data(), y.data(), z.data(), dx.data(), dy.data(), dz.data(), k);
}

// -----------------------------------------------------------------------------
void PrintProgress(ostream &os, int i, double error, bool flush = true)
{
  const streamsize w = os.width(0);
  const streamsize p = os.precision(5);
  const ios::fmtflags f = os.flags();

  os << setw(3) << i << ". RMS error = " << setprecision(5) << error << "\n";
  if (flush) os.flush();

  os.width(w);
  os.precision(p);
  os.flags(f);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(0);

  Array<string> target_names;
  Array<string> source_names;

  TransformationModel model = TM_Rigid;
  string dofin_name;
  string dofout_name;

  PointCorrespondence::TypeId ctype = PointCorrespondence::ClosestPoint;
  Array<string> feature_name;
  Array<double> feature_weight;
  ParameterList param;

  int iterations = 100;
  double epsilon = 0.01;
  bool inverse = false;
  bool symmetric = false;

  verbose = 1;

  for (ALL_OPTIONS) {
    if (OPTION("-t") || OPTION("-target")) {
      do {
        target_names.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-s") || OPTION("-source")) {
      do {
        source_names.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-i") || OPTION("-dofin")) {
      dofin_name = ARGUMENT;
    }
    else if (OPTION("-o") || OPTION("-dofout")) {
      dofout_name = ARGUMENT;
    }
    else if (OPTION("-m") || OPTION("-model")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, model)) FatalError("Invalid -model argument: " << arg);
    }
    else if (OPTION("-n") || OPTION("-iterations")) {
      PARSE_ARGUMENT(iterations);
    }
    else if (OPTION("-c") || OPTION("-cor") || OPTION("-corr") || OPTION("-correspondence")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, ctype)) FatalError("Invalid -correspondence argument: " << arg);
    }
    else if (OPTION("-p") || OPTION("-par") || OPTION("-corpar") || OPTION("-corrpar")) {
      Insert(param, ARGUMENT, ARGUMENT);
    }
    else if (OPTION("-f") || OPTION("-feature")) {
      feature_name.push_back(ARGUMENT);
      if (HAS_ARGUMENT) feature_weight.push_back(atof(ARGUMENT));
      else              feature_weight.push_back(1.0);
    }
    else if (OPTION("-cp") || OPTION("-closest-point")) {
      ctype = PointCorrespondence::ClosestPoint;
    }
    else if (OPTION("-csp") || OPTION("-closest-surface-point") || OPTION("-closest-cell")) {
      ctype = PointCorrespondence::ClosestCell;
    }
    else if (OPTION("-epsilon")) {
      PARSE_ARGUMENT(epsilon);
    }
    else HANDLE_BOOL_OPTION(inverse);
    else HANDLE_BOOL_OPTION(symmetric);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  const int m = static_cast<int>(source_names.size());
  const int n = static_cast<int>(target_names.size());

  // Check required arguments
  if (dofout_name.empty()) {
    FatalError("Option -dofout is required!");
  }
  if (m == 0 || n == 0) {
    FatalError("Options -target and -source are required!");
  }
  if (n > 1 && n != m) {
    FatalError("Either specify a single -target or one target for each -source point set!");
  }
  if (!IsLinear(model)) {
    FatalError("Currently only Rigid, Similarity, and Affine transformation -model supported!");
  }

  // By default, use point coordinates as featurs for point matching
  if (feature_name.empty()) {
    feature_name.push_back("spatial coordinates");
    feature_weight.push_back(1.0);
  }

  // Initialize transformation
  UniquePtr<Transformation> dof(Transformation::New(ToTransformationType(model)));
  if (!dofin_name.empty() && dofin_name != "Id" && dofin_name != "Identity" && dofin_name != "identity") {
    if (verbose) cout << "Reading transformation...", cout.flush();
    UniquePtr<Transformation> dofin(Transformation::New(dofin_name.c_str()));
    const MultiLevelTransformation *mffd = dynamic_cast<const MultiLevelTransformation *>(dofin.get());
    const HomogeneousTransformation *ilin = nullptr;
    if (mffd) {
      ilin = mffd->GetGlobalTransformation();
    } else {
      ilin = dynamic_cast<const HomogeneousTransformation *>(dofin.get());
      if (!ilin) {
        FatalError("Input -dofin must be a linear or multi-level transformation");
      }
    }
    HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(dof.get());
    mirtkAssert(lin != nullptr, "expected transformation for model=" << model << " to be of type HomogeneousTransformation");
    lin->CopyFrom(ilin);
    if (verbose) cout << " done" << endl;
  }

  // Initialize point sets
  if (verbose) cout << "Initialize point sets...", cout.flush();
  Array<RegisteredPointSet> sources(m);
  Array<RegisteredPointSet> targets(n);
  for (int i = 0; i < m; ++i) {
    auto pset = ReadPointSet(source_names[i].c_str());
    if (pset->GetNumberOfPoints() == 0) {
      FatalError("Failed to open source point set or point set contains no points: " << source_names[i]);
    }
    sources[i].InputPointSet(pset);
    if (inverse) {
      sources[i].Transformation(dof.get());
    }
    sources[i].Initialize();
  }
  for (int j = 0; j < n; ++j) {
    auto pset = ReadPointSet(target_names[j].c_str());
    if (pset->GetNumberOfPoints() == 0) {
      FatalError("Failed to open target point set or point set contains no points: " << target_names[j]);
    }
    targets[j].InputPointSet(pset);
    if (!inverse) {
      targets[j].Transformation(dof.get());
    }
    targets[j].Initialize();
  }
  if (verbose) cout << " done" << endl;

  // Initialize correspondence maps
  if (verbose) cout << "Initialize correspondence maps...", cout.flush();
  Array<UniquePtr<PointCorrespondence>> cmaps(sources.size());
  for (int i = 0; i < m; ++i) {
    const int j = TargetIndex(m, n, i);
    auto &cmap = cmaps[i];
    cmap.reset(PointCorrespondence::New(ctype));
    cmap->FromTargetToSource(true);
    cmap->FromSourceToTarget(symmetric);
    cmap->Parameter(param);
    cmap->Source(&sources[i]);
    cmap->Target(&targets[j]);
    for (size_t f = 0; f < feature_name.size(); ++f) {
      cmap->AddFeature(feature_name[f].c_str(), feature_weight[f]);
    }
    cmap->Initialize();
  }
  if (verbose) cout << " done" << endl;

  // Iterate least squares fitting
  double error, last_error = numeric_limits<double>::infinity();
  for (int iter = 0; true; ++iter) {
    // Update correspondences
    Update(targets, sources, cmaps);
    // Check for convergence
    error = EvaluateRMSError(targets, sources, cmaps, symmetric);
    if (verbose) PrintProgress(cout, iter, error);
    if (last_error - error < epsilon) {
      if (verbose) {
        cout << "Converged after " << iter << " iterations." << endl;
      }
      break;
    }
    last_error = error;
    if (iter >= iterations) {
      if (verbose) {
        cout << "Terminated after " << iter << " iterations." << endl;
      }
      break;
    }
    // Update transformation
    Fit(dof.get(), targets, sources, cmaps, symmetric, inverse);
  }

  // Write resulting transformation
  dof->Write(dofout_name.c_str());

  return 0;
}
