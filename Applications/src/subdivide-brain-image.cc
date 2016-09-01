/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/Vector3.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/PointSet.h"
#include "mirtk/GenericImage.h"
#include "mirtk/ConnectedComponents.h"
#include "mirtk/Dilation.h"
#include "mirtk/Erosion.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " [<input>] <output> [options]\n";
  cout << "Usage: " << name << " [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  This program reads a structural brain segmentation and derives from it\n";
  cout << "  a segmentation of the brain volume into the following output labels.\n";
  cout << "  This output segmentation can then be used to reconstruct topologically\n";
  cout << "  correct (i.e., closed genus-0) surfaces of the cGM/WM interface for the\n";
  cout << "  left and right hemisphere, where subcortical and deep brain structures\n";
  cout << "  are enclosed by these so-called white surfaces. Additionally, the\n";
  cout << "  brainstem+cerebellum segment can be respresented by another closed\n";
  cout << "  surface mesh. The union of these reconstructured surfaces encloses\n";
  cout << "  the entire brain volume, yet excluding cortical grey matter.\n";
  cout << "  By deforming the joint brain surface towards the cGM/CSF interface,\n";
  cout << "  the pial surface which encloses the entire brain volume including\n";
  cout << "  subcortical structures can be obtained. The right/left hemisphere\n";
  cout << "  assignment of cortical grey matter follows from the point correspondences\n";
  cout << "  between white and pial surfaces, respectively, the RH/LH label may\n";
  cout << "  be assigned to white surface mesh nodes upon merging the right/left\n";
  cout << "  white surface meshes. See merge-surfaces -source-array option.\n";
  cout << "\n";
  cout << "  Output labels:\n";
  cout << "  - 0: Background\n";
  cout << "  - 1: Cortical grey matter\n";
  cout << "  - 2: Cerebral white matter and deep brain structures of right hemisphere\n";
  cout << "  - 3: Cerebral white matter and deep brain structures of left  hemisphere\n";
  cout << "  - 4: Brainstem, including cerebellum when :option:`-brainstem+cerebellum` is on\n";
  cout << "  - 5: Cerebellum, when :option:`-brainstem+cerebellum` is off\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  label image. See :option:`-input-labels`.\n";
  cout << "  output   Output label image. See :option:`-output-labels`.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -input-labels, -input <file>\n";
  cout << "      Input segmentation label image, use either this option or <input> argument.\n";
  cout << "      For each of the following options which take a <labels> argument,\n";
  cout << "      the respective segmentation mask is derived from this input label\n";
  cout << "      image by merging the specified segments of this brain segmentation.\n";
  cout << "      Multiple integer labels can be specified as the <labels> argument,\n";
  cout << "      each separated by at least one space character. A range of labels\n";
  cout << "      can be specified as <first>..<last>. When the <labels> argument is\n";
  cout << "      a single argument that is neither an integer nor a label range,\n";
  cout << "      the argument is expected to be the name of a binary mask image file.\n";
  cout << "      This option is required when any of the other options has one or more\n";
  cout << "      <labels> arguments. It is ignored otherwise. (default: none)\n";
  cout << "  -output-labels, -output <file>\n";
  cout << "      Output label image, use either this option or <output> argument.\n";
  cout << "  -hemispheres <file>\n";
  cout << "      Hemispheres mask (0: outside, 1: right, 2: left) from\n";
  cout << "      which cutting plane is computed. (default: none)\n";
  cout << "  -right-hemisphere, -rh <file>|<labels>\n";
  cout << "      Mask or labels of structures of right hemisphere used to compute cutting plane.\n";
  cout << "      This option is required when no :option:`-hemispheres` mask is given.\n";
  cout << "  -left-hemisphere, -lh <file>|<labels>\n";
  cout << "      Mask or labels of structures of left  hemisphere used to compute cutting plane.\n";
  cout << "      This option is required when no :option:`-hemispheres` mask is given.\n";
  cout << "  -subcortical, -sb <file>|<labels>\n";
  cout << "      Subcortical / deep brain structures segmentation mask or labels.\n";
  cout << "      When specified, these structures are cut into RH and LH. (default: none)\n";
  cout << "  -white-matter, -wm <file>|<labels>\n";
  cout << "      White matter segmentation mask or labels. (default: none)\n";
  cout << "  -grey-matter, -gm <file>|<labels>\n";
  cout << "      Grey matter segmentation mask or labels. (default: none)\n";
  cout << "  -brainstem, -bs <file>|<labels>\n";
  cout << "      Brainstem segmentation mask or labels. (default: none)\n";
  cout << "  -cerebellum, -cb <file>|<labels>\n";
  cout << "      Cerebellum segmentation mask or labels. (default: none)\n";
  cout << "  -closing-iterations <n>\n";
  cout << "      No. of iterations used to close holes between right/left subcortical,\n";
  cout << "      brainstem, and cerebellum segmentations. (default: 5)\n";
  cout << "  -brainstem+cerebellum, -cerebellum+brainstem, -bs+cb, -cb+bs [on|off]\n";
  cout << "      Whether to merge brainstem and cerebellum. (default: off)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Output brain segmentation labels
enum BrainRegion
{
  BG  = 0, ///< Background/undefined
  GM  = 1, ///< Cortical grey matter
  RH  = 2, ///< Cerebrum (excl. grey matter) of right brain hemisphere
  LH  = 3, ///< Cerebrum (excl. grey matter) of left  brain hemisphere
  CB  = 4, ///< Cerebellum
  BS  = 5, ///< Brainstem / midbrain
  UH  = 6, ///< Unspecified hemisphere
  NUM = 7  ///< Number of brain region labels (i.e., max label + 1)
};

// -----------------------------------------------------------------------------
/// Parameters of cutting plane in normal form
struct Plane
{
  Point   c; ///< Center
  Vector3 n; ///< Normal
  double  b; ///< Offset

  Plane()
  :
    n(0., 0., 0.), b(0.)
  {}

  Plane(const Vector3 &n, const Point &p)
  :
    c(p), n(n), b(-n.Dot(Vector3(c.x, c.y, c.z)))
  {}

  void Center(const Point &p)
  {
    c = p;
    b = - n.Dot(Vector3(c.x, c.y, c.z));
  }

  void Normal(const Vector3 &v)
  {
    n = v;
    b = - n.Dot(Vector3(c.x, c.y, c.z));
  }

  operator bool() const
  {
    return fequal(n.SquaredLength(), 1., 1e-3);
  }

  double SignedDistance(const Point &p) const
  {
    return n.Dot(Vector3(p.x, p.y, p.z)) + b;
  }

  static bool Tangents(const Vector3 &n, Vector3 &e1, Vector3 &e2)
  {
    Vector3 v(n[1], n[2], n[0]);
    e1 = n.Cross(v);
    if (e1.Dot(e1) < 1e-6) {
      v[1] *= -1.;
      e1 = n.Cross(v);
      if (e1.Dot(e1) < 1e-6) return false;
    }
    e2 = n.Cross(e1);
    e1.Normalize();
    e2.Normalize();
    return true;
  }

  bool Tangents(Vector3 &e1, Vector3 &e2) const
  {
    return Tangents(n, e1, e2);
  }
};

// -----------------------------------------------------------------------------
/// Print plane equation
ostream &operator <<(ostream &os, const Plane &plane)
{
  os << "x^T [" << plane.n.x << " " << plane.n.y << " " << plane.n.z << "] + " << plane.b << " = 0";
  return os;
}

// -----------------------------------------------------------------------------
/// Add label range to set of labels
void AddLabels(UnorderedSet<int> &labels, int a, int b)
{
  for (int i = a; i <= b; ++i) labels.insert(i);
}

// -----------------------------------------------------------------------------
/// Parse single label or label range argument
bool ParseLabelRange(const char *arg, int &a, int &b)
{
  const Array<string> parts = Split(arg, "..");
  if (parts.size() == 1) {
    if (FromString(parts[0], a)) {
      b = a;
      return true;
    }
  } else if (parts.size() == 2) {
    if (FromString(parts[0], a) && FromString(parts[1], b)) {
      if (a > b) swap(a, b);
      return true;
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Parse <file>|<labels> command option argument
const char *GetFileNameOrAddLabels(int &OPTIDX, int &argc, char *argv[], UnorderedSet<int> &labels)
{
  int a, b;
  Array<const char *> args;
  do {
    args.push_back(ARGUMENT);
  } while (HAS_ARGUMENT);
  for (const auto &arg : args) {
    if (ParseLabelRange(arg, a, b)) {
      AddLabels(labels, a, b);
    } else {
      if (args.size() == 1) return args[0];
      FatalError("Invalid <labels> argument: " << arg);
    }
  }
  return nullptr;
}

// -----------------------------------------------------------------------------
/// Add segments with specified labels to segmentation mask
void AddLabels(ByteImage &mask, const GreyImage &label_image, const UnorderedSet<int> &labels)
{
  const int nvox = mask.NumberOfVoxels();
  for (int vox = 0; vox < nvox; ++vox) {
    if (labels.find(label_image(vox)) != labels.end()) {
      mask(vox) = 1;
    }
  }
}

// -----------------------------------------------------------------------------
/// Check whether brain regions map contains a certain label
bool Contains(const ByteImage &regions, BrainRegion region)
{
  const int nvox = regions.NumberOfVoxels();
  for (int vox = 0; vox < nvox; ++vox) {
    if (regions(vox) == region) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Initialize boolean array of whether or not a brain region is selected given set
void InitializeSelectionLUT(const UnorderedSet<int> &labels, bool selected[NUM])
{
  for (int i = 0; i < NUM; ++i) selected[i] = false;
  for (auto label : labels) selected[label] = true;
}

// -----------------------------------------------------------------------------
/// Get set of boundary points
void AddBoundaryPoints(PointSet                &points,
                       const ByteImage         &regions,
                       const UnorderedSet<int> &region1,
                       const UnorderedSet<int> &region2)
{
  Point p;
  bool selection1[NUM], selection2[NUM], boundary;
  InitializeSelectionLUT(region1, selection1);
  InitializeSelectionLUT(region2, selection2);
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    if (selection1[regions(i, j, k)]) {
      boundary = false;
      for (int nk = k-1; nk <= k+1; ++nk) {
        if (nk < 0 || nk >= regions.Z()) continue;
        for (int nj = j-1; nj <= j+1; ++nj) {
          if (nj < 0 || nj >= regions.Y()) continue;
          for (int ni = i-1; ni <= i+1; ++ni) {
            if (ni < 0 || ni >= regions.X()) continue;
            if (selection2[regions(ni, nj, nk)]) {
              boundary = true;
              break;
            }
          }
          if (boundary) break;
        }
        if (boundary) break;
      }
      if (boundary) {
        p = Point(i, j, k);
        regions.ImageToWorld(p);
        points.Add(p);
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Get set of boundary points
void AddBoundaryPoints(PointSet &points, const ByteImage &regions,
                       BrainRegion region1, BrainRegion region2)
{
  UnorderedSet<int> set1, set2;
  set1.insert(region1);
  set2.insert(region2);
  AddBoundaryPoints(points, regions, set1, set2);
}

// -----------------------------------------------------------------------------
/// Get set of boundary points
PointSet BoundaryPoints(const ByteImage         &regions,
                        const UnorderedSet<int> &region1,
                        const UnorderedSet<int> &region2)
{
  PointSet points;
  AddBoundaryPoints(points, regions, region1, region2);
  return points;
}

// -----------------------------------------------------------------------------
/// Add world coordinates of brain region points to given point set
void AddPoints(PointSet &points, const ByteImage &regions, BrainRegion region)
{
  Point p;
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    if (regions(i, j, k) == region) {
      p = Point(i, j, k);
      regions.ImageToWorld(p);
      points.Add(p);
    }
  }
}

// -----------------------------------------------------------------------------
/// Principal component analysis of point set
bool PrincipalDirections(const PointSet &points, Vector3 dir[3])
{
  if (points.Size() == 0) return false;

  Point p;
  Matrix3x3 covar(0.);
  const Point c = points.Centroid();
  for (int i = 0; i < points.Size(); ++i) {
    p = points(i) - c;
    for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) {
      covar[r][c] += p[r] * p[c];
    }
  }

  Vector3       axis[3];
  Array<double> eigval(3);
  covar.EigenSolveSymmetric(eigval.data(), axis);

  eigval[0] = abs(eigval[0]);
  eigval[1] = abs(eigval[1]);
  eigval[2] = abs(eigval[2]);
  Array<int> order = DecreasingOrder(eigval);
  if (eigval[order[2]] < 1e-6) return false;

  dir[0] = axis[order[0]], dir[0].Normalize();
  dir[1] = axis[order[1]], dir[1].Normalize();
  dir[2] = axis[order[2]], dir[2].Normalize();
  return true;
}

// -----------------------------------------------------------------------------
/// Principal component analysis of structure shape
bool PrincipalDirections(const ByteImage &regions, BrainRegion region, Vector3 dir[3])
{
  PointSet points;
  AddPoints(points, regions, region);
  return PrincipalDirections(points, dir);
}

// -----------------------------------------------------------------------------
/// Get principal direction with most spread
Vector3 PrincipalDirection(const ByteImage &regions, BrainRegion region)
{
  Vector3 dir[3] = {0.};
  PrincipalDirections(regions, region, dir);
  return dir[0];
}

// -----------------------------------------------------------------------------
/// Mask of brain regions cut by specified plane
ByteImage Cut(const ByteImage &regions, const Plane &plane, double r)
{
  ImageAttributes attr;

  Vector3 e1, e2;
  plane.Tangents(e1, e2);

  const double ds = min(min(regions.XSize(), regions.YSize()), regions.ZSize());
  attr._dx = attr._dy = attr._dz = ds;

  attr._x = attr._y = 2 * iceil(r / ds) + 1;
  attr._z = 1;

  attr._xorigin = plane.c.x;
  attr._yorigin = plane.c.y;
  attr._zorigin = plane.c.z;

  attr._xaxis[0] = e1.x;
  attr._xaxis[1] = e1.y;
  attr._xaxis[2] = e1.z;

  attr._yaxis[0] = e2.x;
  attr._yaxis[1] = e2.y;
  attr._yaxis[2] = e2.z;

  attr._zaxis[0] = plane.n.x;
  attr._zaxis[1] = plane.n.y;
  attr._zaxis[2] = plane.n.z;

  const Matrix w2i = attr.GetWorldToImageMatrix();
  const Matrix i2w = attr.GetImageToWorldMatrix();

  attr._w2i = &w2i;
  attr._i2w = &i2w;

  GenericNearestNeighborInterpolateImageFunction<ByteImage> nn;
  nn.Input(&regions);
  nn.DefaultValue(BG);
  nn.Initialize();

  double x, y, z;
  ByteImage cut(attr);
  for (int j = 0; j < cut.Y(); ++j)
  for (int i = 0; i < cut.X(); ++i) {
    x = i, y = j, z = 0.;
    attr.LatticeToWorld(x, y, z);
    regions.WorldToImage(x, y, z);
    cut(i, j) = nn.Evaluate(x, y, z);
  }

  return cut;
}

// -----------------------------------------------------------------------------
/// Calculate diameter of brainstem oriented bounding box
double MaxBrainstemDiameter(const ByteImage &regions, const Plane &plane)
{
  Vector3 e1, e2, p;
  plane.Tangents(e1, e2);
  double x, x1 = inf, x2 = -inf;
  double y, y1 = inf, y2 = -inf;
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    if (regions(i, j, k) == BS) {
      p = Vector3(i, j, k);
      regions.ImageToWorld(p.x, p.y, p.z);
      x = e1.Dot(p);
      y = e2.Dot(p);
      x1 = min(x1, x);
      x2 = max(x2, x);
      y1 = min(y1, y);
      y2 = max(y2, y);
    }
  }
  return max(y2 - y1, x2 - x1);
}

// -----------------------------------------------------------------------------
/// Test if specified cut is a desired brainstem cut
bool ValidBrainstemCut(const ByteImage &cut)
{
  int nbs = 0;
  for (int j = 0; j < cut.Y(); ++j)
  for (int i = 0; i < cut.X(); ++i) {
    if (cut(i, j) == BS) {
      ++nbs;
      for (int nj = j-1; nj <= j+1; ++nj)
      for (int ni = i-1; ni <= i+1; ++ni) {
        if (cut(ni, nj) != BS && cut(ni, nj) != BG) {
          return false;
        }
      }
    }
  }
  return nbs > 0;
}

// -----------------------------------------------------------------------------
/// Cutting plane separating brainstem from cerebrum
Plane BrainstemCuttingPlane(const ByteImage &regions)
{
  Plane plane;

  // Sets of brainstem brain regions labels
  UnorderedSet<int> set1, set2;
  set1.insert(BS);

  set2.insert(RH);
  set2.insert(LH);
  set2.insert(UH);

  // Plane normal direction is longitudinal axis of brainstem
  plane.Normal(PrincipalDirection(regions, BS));
  if (!plane) return plane;

  // Choose plane center point which is closest to the boundary with the cerebrum
  plane.Center(BoundaryPoints(regions, set1, set2).Centroid());

  // Adjust normal direction such that most of brainstem is in positive half space
  Point p;
  int npos = 0, nneg = 0;
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    if (regions(i, j, k) == BS) {
      p = Point(i, j, k);
      regions.ImageToWorld(p);
      if (plane.SignedDistance(p) >= 0.) ++npos;
      else                               ++nneg;
    }
  }
  if (nneg > npos) plane.Normal(plane.n * -1.);

  // Move center point away from cerebrum until a clean cutting point is found,
  // where no brainstem voxel is connected to a non-brainstem voxel.
  // The other side of the cut will contain part of the midbrain.
  // This is ok for our purposes, where we only want to divide the brain volume
  // into clean disjoint regions whose surface can initially be reconstructed
  // by shrinking a genus-0 surface (e.g., sphere, convex hull) towards
  // each hemisphere of the (WM) segmentation individually such that these
  // shrink well into the sulci in the medial fissure. The left and right white
  // surfaces are afterwards combined with the brainstem+cerebellum surface
  // to form a single genus-0 surface enclosing the entire brain volume.
  const int max_iter = 100;

  const Point  c = plane.c;
  const double r = MaxBrainstemDiameter(regions, plane);
  const double h = min(min(regions.XSize(), regions.YSize()), regions.ZSize());

  int iter;
  for (iter = 0; iter < max_iter; ++iter) {
    plane.Center(plane.c + h * plane.n);
    ByteImage cut = Cut(regions, plane, r);
    if (Contains(cut, CB)) {
      plane.Center(plane.c - h * plane.n);
      break;
    }
    if (ValidBrainstemCut(cut)) break;
  }
  if (iter == max_iter) {
    plane.Center(c);
  }

  if (debug) {
    Cut(regions, plane, r).Write("debug_brainstem_cut.nii.gz");
  }
  return plane;
}

// -----------------------------------------------------------------------------
/// Naive cutting plane calculation assuming symmetric brain hemispheres mask
Plane SymmetricCuttingPlane(const ByteImage &regions, BrainRegion a, BrainRegion b)
{
  int   comp;
  Point center[2], p;
  int   numvox[2] = {0};
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    const auto &region = regions(i, j, k);
    if      (region == a) comp = 0;
    else if (region == b) comp = 1;
    else continue;
    p = Point(i, j, k);
    regions.ImageToWorld(p);
    center[comp] += p;
    numvox[comp] += 1;
  }
  if (numvox[0] == 0 || numvox[1] == 0) {
    return Plane();
  }
  center[0] /= numvox[0];
  center[1] /= numvox[1];
  Plane plane;
  plane.n[0] = center[0].x - center[1].x;
  plane.n[1] = center[0].y - center[1].y;
  plane.n[2] = center[0].z - center[1].z;
  plane.n.Normalize();
  plane.Center((center[0] + center[1]) / 2.0);
  return plane;
}

// -----------------------------------------------------------------------------
/// Cutting plane separating left/right brain hemisphere clusters
Plane MedialCuttingPlane(const ByteImage &regions)
{
  Plane plane;
  PointSet points;
  AddBoundaryPoints(points, regions, RH, LH);
  AddBoundaryPoints(points, regions, LH, RH);
  if (verbose) {
    cout << "No. of medial points = " << points.Size() << endl;
  }
  if (points.Size() < 1000) {
    plane = SymmetricCuttingPlane(regions, RH, LH);
  } else {
    Vector3 dir[3];
    if (PrincipalDirections(points, dir)) {
      plane.Normal(dir[2]);
      plane.Center(points.Centroid());
    }
  }
  // TODO: Minimize misclassification of RH/LH based on cutting plane by
  //       trying different rotation angles and slight offsets of the
  //       cutting plane. See also merge-surfaces.cc: FindCuttingPlane.
  return plane;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(0);

  const char *labels_name = nullptr; // Input  labels image
  const char *output_name = nullptr; // Output labels image

  const char *hemis_name  = nullptr; // RH/LH cerebrum label image
  const char *rhmask_name = nullptr; // RH cerebrum mask
  const char *lhmask_name = nullptr; // LH cerebrum mask
  const char *wmmask_name = nullptr; // Cerebral WM mask
  const char *gmmask_name = nullptr; // Cerebral GM mask
  const char *sbmask_name = nullptr; // Subcortical structures mask
  const char *bsmask_name = nullptr; // Brainstem mask
  const char *cbmask_name = nullptr; // Cerebellum mask

  UnorderedSet<int> rhmask_labels;
  UnorderedSet<int> lhmask_labels;
  UnorderedSet<int> wmmask_labels;
  UnorderedSet<int> gmmask_labels;
  UnorderedSet<int> sbmask_labels;
  UnorderedSet<int> bsmask_labels;
  UnorderedSet<int> cbmask_labels;

  int closing_iter = 5;
  bool merge_bs_cb = false;

  if (NUM_POSARGS == 1) {
    output_name = POSARG(1);
  } else if (NUM_POSARGS == 2) {
    labels_name = POSARG(1);
    output_name = POSARG(2);
  } else if (NUM_POSARGS > 2) {
    FatalError("Too many positional arguments!");
  }

  for (ALL_OPTIONS) {
    if (OPTION("-input-labels") || OPTION("-input")) {
      if (labels_name) {
        FatalError("Use either <input> argument or -input-labels option, not both!");
      }
      labels_name = ARGUMENT;
    }
    else if (OPTION("-output-labels") || OPTION("-output")) {
      if (output_name) {
        FatalError("Use either <output> argument or -output-labels option, not both!");
      }
      output_name = ARGUMENT;
    }
    else if (OPTION("-hemispheres")) {
      hemis_name = ARGUMENT;
    }
    else if (OPTION("-right-hemisphere") || OPTION("-rh")) {
      rhmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, rhmask_labels);
    }
    else if (OPTION("-left-hemisphere") || OPTION("-lh")) {
      lhmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, lhmask_labels);
    }
    else if (OPTION("-white-matter") || OPTION("-wm")) {
      wmmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, wmmask_labels);
    }
    else if (OPTION("-grey-matter") || OPTION("-gm")) {
      gmmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, gmmask_labels);
    }
    else if (OPTION("-subcortical") || OPTION("-sb")) {
      sbmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, sbmask_labels);
    }
    else if (OPTION("-brainstem") || OPTION("-bs")) {
      bsmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, bsmask_labels);
    }
    else if (OPTION("-cerebellum") || OPTION("-cb")) {
      cbmask_name = GetFileNameOrAddLabels(OPTIDX, argc, argv, cbmask_labels);
    }
    else if (OPTION("-brainstem+cerebellum") || OPTION("-bs+cb") ||
             OPTION("-cerebellum+brainstem") || OPTION("-cb+bs")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(merge_bs_cb);
      else merge_bs_cb = true;
    }
    else if (OPTION("-closing-iterations")) {
      PARSE_ARGUMENT(closing_iter);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Output labels image name required
  if (output_name == nullptr) {
    FatalError("Missing -output-labels file name!");
  }

  // Initialize I/O library
  InitializeIOLibrary();

  // Read input labels image
  GreyImage labels;
  if (!rhmask_labels.empty() ||
      !lhmask_labels.empty() ||
      !wmmask_labels.empty() ||
      !gmmask_labels.empty() ||
      !sbmask_labels.empty() ||
      !bsmask_labels.empty() ||
      !cbmask_labels.empty()) {
    if (labels_name == nullptr) {
      FatalError("Input -labels image required when using <labels> option argument(s)!");
    }
    labels.Read(labels_name);
  } else {
    if ((verbose > 0 || debug > 0) && labels_name != nullptr) {
      Warning("Input -labels image unused, using provided segmentation masks instead.");
    }
  }

  ImageAttributes attr = labels.Attributes();
  int             nvox = attr.NumberOfSpatialPoints();

  // Read input segmentation masks
  ByteImage rhmask, lhmask, sbmask, bsmask, cbmask;

  if (rhmask_name) {
    rhmask.Read(rhmask_name);
    if (!labels.IsEmpty() && rhmask.Attributes() != attr) {
      FatalError("Right hemispheres mask must have same attributes as input labels image!");
    }
  }
  if (!rhmask_labels.empty()) {
    if (rhmask.IsEmpty()) rhmask.Initialize(attr);
    AddLabels(rhmask, labels, rhmask_labels);
  }
  if (lhmask_name) {
    lhmask.Read(lhmask_name);
    if (!labels.IsEmpty() && lhmask.Attributes() != attr) {
      FatalError("Left hemispheres mask must have same attributes as input labels image!");
    }
  }
  if (!lhmask_labels.empty()) {
    if (lhmask.IsEmpty()) lhmask.Initialize(attr);
    AddLabels(lhmask, labels, lhmask_labels);
  }
  if (hemis_name && (rhmask.IsEmpty() || lhmask.IsEmpty())) {
    ByteImage hemis(hemis_name);
    if (!labels.IsEmpty() && hemis.Attributes() != attr) {
      FatalError("Hemispheres mask must have same attributes as input labels image!");
    }
    if (rhmask.IsEmpty()) {
      rhmask.Initialize(hemis.Attributes());
      const int nvox = hemis.NumberOfVoxels();
      for (int vox = 0; vox < nvox; ++vox) {
        rhmask(vox) = (hemis(vox) == 1);
      }
    }
    if (lhmask.IsEmpty()) {
      lhmask.Initialize(hemis.Attributes());
      const int nvox = hemis.NumberOfVoxels();
      for (int vox = 0; vox < nvox; ++vox) {
        lhmask(vox) = (hemis(vox) == 2);
      }
    }
  }
  if (rhmask.IsEmpty() || lhmask.IsEmpty()) {
    FatalError("Right and left hemisphere mask <file> or <labels> required!");
  }
  if (labels.IsEmpty()) {
    attr = rhmask.Attributes();
    nvox = attr.NumberOfSpatialPoints();
  }
  if (lhmask.Attributes() != attr) {
    FatalError("Attributes of right and left hemispheres masks must be identical!");
  }
  if (sbmask_name) {
    sbmask.Read(sbmask_name);
    if (sbmask.Attributes() != attr) {
      FatalError("Attributes of subcortical mask must match those of hemisphere masks!");
    }
  }
  if (!sbmask_labels.empty()) {
    if (sbmask.IsEmpty()) sbmask.Initialize(attr);
    AddLabels(sbmask, labels, sbmask_labels);
  }
  if (bsmask_name) {
    bsmask.Read(bsmask_name);
    if (bsmask.Attributes() != attr) {
      FatalError("Attributes of brainstem mask must match those of hemisphere masks!");
    }
  }
  if (!bsmask_labels.empty()) {
    if (bsmask.IsEmpty()) bsmask.Initialize(attr);
    AddLabels(bsmask, labels, bsmask_labels);
  }
  if (cbmask_name) {
    cbmask.Read(cbmask_name);
    if (cbmask.Attributes() != attr) {
      FatalError("Attributes of cerebellum mask must match those of hemisphere masks!");
    }
  }
  if (!cbmask_labels.empty()) {
    if (cbmask.IsEmpty()) cbmask.Initialize(attr);
    AddLabels(cbmask, labels, cbmask_labels);
  }

  // Initial brain regions image used to find cutting planes
  ByteImage regions(attr, 1);
  for (int vox = 0; vox < nvox; ++vox) {
    if      (rhmask(vox)) regions(vox) = RH;
    else if (lhmask(vox)) regions(vox) = LH;
  }
  if (!bsmask.IsEmpty()) {
    for (int vox = 0; vox < nvox; ++vox) {
      if (bsmask(vox)) regions(vox) = BS;
    }
  }
  if (!cbmask.IsEmpty()) {
    for (int vox = 0; vox < nvox; ++vox) {
      if (cbmask(vox)) regions(vox) = CB;
    }
  }

  if (closing_iter > 0) {
    // Close holes in subcortical mask
    if (!sbmask.IsEmpty()) {
      BinaryImage closed(attr, 1);
      for (int vox = 0; vox < nvox; ++vox) {
        closed(vox) = sbmask(vox);
      }
      Dilate<BinaryPixel>(&closed, closing_iter,     CONNECTIVITY_18);
      Erode <BinaryPixel>(&closed, closing_iter + 1, CONNECTIVITY_18);
      for (int vox = 0; vox < nvox; ++vox) {
        if (regions(vox) == BG && (closed(vox) != 0 || sbmask(vox) != 0)) {
          regions(vox) = UH;
        }
      }
    }

    // Close holes in brainstem+cerebellum region
    if (!bsmask.IsEmpty()) {
      BinaryImage closed(attr, 1);
      for (int vox = 0; vox < nvox; ++vox) {
        closed(vox) = BinaryPixel(regions(vox) == BS ? 1 : 0);
      }
      Dilate<BinaryPixel>(&closed, closing_iter,     CONNECTIVITY_18);
      Erode <BinaryPixel>(&closed, closing_iter + 1, CONNECTIVITY_18);
      int ncb = 0;
      for (int vox = 0; vox < nvox; ++vox) {
        const auto &region = regions(vox);
        if (region == BG) {
          if (closed(vox) != 0) {
            regions(vox) = BS;
          }
        } else if (region == CB) {
          closed(vox) = 1;
          ++ncb;
        } else if (region != BS) {
          closed(vox) = 0;
        }
      }
      if (ncb > 0) {
        Dilate<BinaryPixel>(&closed, closing_iter,     CONNECTIVITY_18);
        Erode <BinaryPixel>(&closed, closing_iter + 1, CONNECTIVITY_18);
      }
      for (int vox = 0; vox < nvox; ++vox) {
        const auto &region = regions(vox);
        if (region == BG) {
          if (closed(vox) != 0) {
            regions(vox) = CB;
          }
        } else if (region == UH) {
          closed(vox) = 1;
        } else if (region != BS && region != CB) {
          closed(vox) = 0;
        }
      }
      Dilate<BinaryPixel>(&closed, closing_iter,     CONNECTIVITY_18);
      Erode <BinaryPixel>(&closed, closing_iter + 1, CONNECTIVITY_18);
      for (int vox = 0; vox < nvox; ++vox) {
        if (regions(vox) == BG && closed(vox) != 0) {
          regions(vox) = UH;
        }
      }
    }
  }

  if (debug) {
    regions.Write("debug_regions.nii.gz");
  }

  // Determine cutting planes
  Plane bs_plane = BrainstemCuttingPlane(regions);
  Plane rl_plane = MedialCuttingPlane(regions);
  if (verbose) {
    if (bs_plane) cout << "Brainstem cutting plane equation: " << bs_plane << endl;
    if (rl_plane) cout << "Medial cutting plane equation:    " << rl_plane << endl;
  }

  // Cut WM segmentation into RH and LH
  if (wmmask_name || !wmmask_labels.empty()) {
    ByteImage wm;
    if (wmmask_name) wm.Read(wmmask_name);
    if (!wmmask_labels.empty()) {
      if (wm.IsEmpty()) wm.Initialize(attr);
      AddLabels(wm, labels, wmmask_labels);
    }

    ConnectedComponents<GreyPixel> cc;
    GreyImage wm_lh(wm.Attributes());
    GreyImage wm_rh(wm.Attributes());
    GreyImage wm_lh_cc, wm_rh_cc;

    Point p;
    for (int k = 0; k < wm.Z(); ++k)
    for (int j = 0; j < wm.Y(); ++j)
    for (int i = 0; i < wm.X(); ++i) {
      if (wm(i, j, k) != 0) {
        p = Point(i, j, k);
        wm.ImageToWorld(p);
        if (rl_plane.SignedDistance(p) < 0.) wm_lh(i, j, k) = 1;
        else                                 wm_rh(i, j, k) = 1;
      }
    }

    cc.Input (&wm_rh);
    cc.Output(&wm_rh_cc);
    cc.Run();

    cc.Input (&wm_lh);
    cc.Output(&wm_lh_cc);
    cc.Run();

    for (int vox = 0; vox < nvox; ++vox) {
      wm_rh(vox) = (wm_rh_cc(vox) == 1 || wm_lh_cc(vox) > 1 ? 1 : 0);
      wm_lh(vox) = (wm_lh_cc(vox) == 1 || wm_rh_cc(vox) > 1 ? 1 : 0);
    }

    cc.Input (&wm_rh);
    cc.Output(&wm_rh_cc);
    cc.Run();

    cc.Input (&wm_lh);
    cc.Output(&wm_lh_cc);
    cc.Run();

    for (int vox = 0; vox < nvox; ++vox) {
      if      (wm_rh_cc(vox) == 1) wm(vox) = RH;
      else if (wm_lh_cc(vox) == 1) wm(vox) = LH;
      else                         wm(vox) = BG;
    }

    // Replace initial RH/LH labels by new WM segmentation based RH/LH labels
    for (int vox = 0; vox < nvox; ++vox) {
      if (wm(vox) != BG) {
        regions(vox) = wm(vox);
      } else if (regions(vox) == RH || regions(vox) == LH) {
        regions(vox) = BG;
      }
    }
  }

  // Add cortical GM segmentation
  if (gmmask_name || !gmmask_labels.empty()) {
    ByteImage gm;
    if (gmmask_name) gm.Read(gmmask_name);
    if (!gmmask_labels.empty()) {
      if (gm.IsEmpty()) gm.Initialize(attr);
      AddLabels(gm, labels, gmmask_labels);
    }
    ConnectedComponents<BytePixel> cc;
    cc.Input (&gm);
    cc.Output(&gm);
    for (int vox = 0; vox < nvox; ++vox) {
      if (gm(vox) == 1) {
        regions(vox) = GM;
      }
    }
  }

  // Cut subcortical structures and brainstem into RH, LH, and BS+CB
  Point p;
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    p = Point(i, j, k);
    regions.ImageToWorld(p);
    auto &region = regions(i, j, k);
    if (region == BS || region == CB) {
      if (bs_plane.SignedDistance(p) < 0.) {
        if (rl_plane.SignedDistance(p) < 0.) {
          region = LH;
        } else {
          region = RH;
        }
      } else if (merge_bs_cb) {
        region = min(BS, CB);
      }
    } else if (region == UH || (!sbmask.IsEmpty() && sbmask(i, j, k) != 0)) {
      if (rl_plane.SignedDistance(p) < 0.) {
        region = LH;
      } else {
        region = RH;
      }
    }
  }

  // Write output labels
  regions.Write(output_name);
  return 0;
}
