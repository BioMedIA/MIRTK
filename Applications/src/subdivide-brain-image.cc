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
#include "mirtk/Closing.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/EuclideanDistanceTransform.h"

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
  cout << "  - 4: Brainstem, including cerebellum when :option:`-brainstem-and-cerebellum` is on\n";
  cout << "  - 5: Cerebellum, when :option:`-brainstem-and-cerebellum` is off\n";
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
  cout << "  -closing, -closing-iterations <n>\n";
  cout << "      No. of iterations used to close holes between right/left subcortical,\n";
  cout << "      brainstem, and cerebellum segmentations. (default: 0)\n";
  cout << "  -subcortical-closing <n>\n";
  cout << "      No. of iterations used to close holes between right/left subcortical segmentation. (default: 0)\n";
  cout << "  -brainstem-closing <n>\n";
  cout << "      No. of iterations used to close holes in brainstem segmentation. (default: 0)\n";
  cout << "  -cerebellum-closing <n>\n";
  cout << "      No. of iterations used to close holes in cerebellum segmentation. (default: 0)\n";
  cout << "  -brainstem-and-cerebellum, -cerebellum-and-brainstem, -bscb, -cbbs [on|off]\n";
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
    c(p), n(n)
  {
    this->n.Normalize();
    b = - this->n.Dot(Vector3(c._x, c._y, c._z));
  }

  void Center(const Point &p)
  {
    c = p;
    b = - n.Dot(Vector3(c._x, c._y, c._z));
  }

  void Normal(const Vector3 &v)
  {
    n = v;
    n.Normalize();
    b = - n.Dot(Vector3(c._x, c._y, c._z));
  }

  operator bool() const
  {
    return fequal(n.SquaredLength(), 1., 1e-3);
  }

  double SignedDistance(const Point &p) const
  {
    return n.Dot(Vector3(p._x, p._y, p._z)) + b;
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
  os << "x^T [" << plane.n._x << " " << plane.n._y << " " << plane.n._z << "] + " << plane.b << " = 0";
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
/// Get set of boundary voxel
void AddBoundaryPoints(PointSet                &points,
                       const ByteImage         &regions,
                       const UnorderedSet<int> &region1,
                       const UnorderedSet<int> &region2,
                       bool wc = true)
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
        if (wc) regions.ImageToWorld(p);
        points.Add(p);
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Get set of boundary points
void AddBoundaryPoints(PointSet &points, const ByteImage &regions,
                       BrainRegion region1, BrainRegion region2,
                       bool wc = true)
{
  UnorderedSet<int> set1, set2;
  set1.insert(region1);
  set2.insert(region2);
  AddBoundaryPoints(points, regions, set1, set2, wc);
}

// -----------------------------------------------------------------------------
/// Get set of boundary points
PointSet BoundaryPoints(const ByteImage         &regions,
                        const UnorderedSet<int> &region1,
                        const UnorderedSet<int> &region2,
                        bool wc = true)
{
  PointSet points;
  AddBoundaryPoints(points, regions, region1, region2, wc);
  return points;
}

// -----------------------------------------------------------------------------
/// Get set of boundary points
PointSet BoundaryPoints(const ByteImage &regions, BrainRegion region1, BrainRegion region2, bool wc = true)
{
  PointSet points;
  AddBoundaryPoints(points, regions, region1, region2, wc);
  return points;
}

// -----------------------------------------------------------------------------
/// Add world coordinates of brain region points to given point set
void AddPoints(PointSet &points, const ByteImage &regions, BrainRegion region, bool wc = true)
{
  Point p;
  for (int k = 0; k < regions.Z(); ++k)
  for (int j = 0; j < regions.Y(); ++j)
  for (int i = 0; i < regions.X(); ++i) {
    if (regions(i, j, k) == region) {
      p = Point(i, j, k);
      if (wc) regions.ImageToWorld(p);
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

  attr._xorigin = plane.c._x;
  attr._yorigin = plane.c._y;
  attr._zorigin = plane.c._z;

  attr._xaxis[0] = e1._x;
  attr._xaxis[1] = e1._y;
  attr._xaxis[2] = e1._z;

  attr._yaxis[0] = e2._x;
  attr._yaxis[1] = e2._y;
  attr._yaxis[2] = e2._z;

  attr._zaxis[0] = plane.n._x;
  attr._zaxis[1] = plane.n._y;
  attr._zaxis[2] = plane.n._z;

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
    cut(i, j) = static_cast<BytePixel>(nn.Evaluate(x, y, z));
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
      regions.ImageToWorld(p._x, p._y, p._z);
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
Plane BrainstemCuttingPlane(const ByteImage &regions, const Vector3 &n)
{
  // Plane normal direction is longitudinal axis of brainstem
  Plane plane;
  plane.Normal(n);
  if (!plane) return plane;

  // Sets of brainstem brain regions labels
  UnorderedSet<int> set1, set2;
  set1.insert(BS);

  set2.insert(RH);
  set2.insert(LH);
  set2.insert(UH);

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
/// Cutting plane separating brainstem from cerebrum
Plane BrainstemCuttingPlane(const ByteImage &regions)
{
  // Plane normal direction is longitudinal axis of brainstem
  return BrainstemCuttingPlane(regions, PrincipalDirection(regions, BS));
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
  plane.n[0] = center[0]._x - center[1]._x;
  plane.n[1] = center[0]._y - center[1]._y;
  plane.n[2] = center[0]._z - center[1]._z;
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

// -----------------------------------------------------------------------------
bool TouchesOtherHemisphere(const GreyImage &cc, const GreyImage &other, int vox)
{
  const GreyPixel label = cc(vox);
  int i, j, k;
  Stack<int> active;
  UnorderedSet<int> visited;
  active.push(vox);
  while (!active.empty()) {
    vox = active.top();
    active.pop();
    visited.insert(vox);
    cc.IndexToVoxel(vox, i, j, k);
    for (int nk = k-1; nk <= k+1; ++nk)
    for (int nj = j-1; nj <= j+1; ++nj)
    for (int ni = i-1; ni <= i+1; ++ni) {
      if (cc.IsInside(ni, nj, nk)) {
        if (other(ni, nj, nk) == 1) {
          return true;
        } else if (cc(ni, nj, nk) == label) {
          vox = cc.VoxelToIndex(ni, nj, nk);
          if (visited.find(vox) == visited.end()) {
            active.push(vox);
          }
        }
      }
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
void FixHemisphereLabels(ByteImage &regions)
{
  const auto attr = regions.Attributes();
  const int  nvox = regions.NumberOfVoxels();

  ConnectedComponents<GreyPixel> cc;
  GreyImage wm_lh(attr), wm_rh(attr);
  GreyImage wm_lh_cc, wm_rh_cc;

  for (int vox = 0; vox < nvox; ++vox) {
    if      (regions(vox) == RH) wm_rh(vox) = 1;
    else if (regions(vox) == LH) wm_lh(vox) = 1;
  }

  cc.Input (&wm_rh);
  cc.Output(&wm_rh_cc);
  cc.Run();

  cc.Input (&wm_lh);
  cc.Output(&wm_lh_cc);
  cc.Run();

  for (int vox = 0; vox < nvox; ++vox) {
    if (regions(vox) == RH && wm_rh_cc(vox) != 1) {
      if (TouchesOtherHemisphere(wm_rh_cc, wm_lh_cc, vox)) {
        regions(vox) = LH;
      } else {
        regions(vox) = BG;
      }
    } else if (regions(vox) == LH && wm_lh_cc(vox) != 1) {
      if (TouchesOtherHemisphere(wm_lh_cc, wm_rh_cc, vox)) {
        regions(vox) = RH;
      } else {
        regions(vox) = BG;
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Resample regions mask such that image axes are aligned with orthogonal cutting planes
ByteImage Resample(const ByteImage &regions, const Plane &rl_plane, const Plane &bs_plane,
                   int xmargin = 0, int ymargin = 0, int zmargin = 0,
                   const char *fgmask_name = nullptr)
{
  // Determine RAS attributes
  const double ds = min(min(regions.XSize(), regions.YSize()), regions.ZSize());

  Vector3 xaxis = rl_plane.n;              // L -> R
  Vector3 yaxis = xaxis.Cross(bs_plane.n); // P -> A
  Vector3 zaxis = xaxis.Cross(yaxis);      // I -> S

  xaxis.Normalize();
  yaxis.Normalize();
  zaxis.Normalize();

  int bi[2], bj[2], bk[2];
  regions.BoundingBox(bi[0], bj[0], bk[0], bi[1], bj[1], bk[1]);

  Point  p, q;
  double bounds[6] = {+inf, -inf, +inf, -inf, +inf, -inf};
  double limits[6] = {+inf, -inf, +inf, -inf, +inf, -inf};
  for (int ck = 0; ck < 2; ++ck)
  for (int cj = 0; cj < 2; ++cj)
  for (int ci = 0; ci < 2; ++ci) {
    p = Point(bi[ci], bj[cj], bk[ck]);
    regions.ImageToWorld(p);
    bounds[0] = min(bounds[0], p._x);
    bounds[1] = max(bounds[1], p._x);
    bounds[2] = min(bounds[2], p._y);
    bounds[3] = max(bounds[3], p._y);
    bounds[4] = min(bounds[4], p._z);
    bounds[5] = max(bounds[5], p._z);
    q._x = xaxis.Dot(p);
    q._y = yaxis.Dot(p);
    q._z = zaxis.Dot(p);
    limits[0] = min(limits[0], q._x);
    limits[1] = max(limits[1], q._x);
    limits[2] = min(limits[2], q._y);
    limits[3] = max(limits[3], q._y);
    limits[4] = min(limits[4], q._z);
    limits[5] = max(limits[5], q._z);
  }

  ImageAttributes attr;
  attr._xorigin = bounds[0] + .5 * (bounds[1] - bounds[0]);
  attr._yorigin = bounds[2] + .5 * (bounds[3] - bounds[2]);
  attr._zorigin = bounds[4] + .5 * (bounds[5] - bounds[4]);
  attr._x       = iceil((limits[1] - limits[0]) / ds);
  attr._y       = iceil((limits[3] - limits[2]) / ds);
  attr._z       = iceil((limits[5] - limits[4]) / ds);
  attr._dx      = ds;
  attr._dy      = ds;
  attr._dz      = ds;
  memcpy(attr._xaxis, xaxis, 3 * sizeof(double));
  memcpy(attr._yaxis, yaxis, 3 * sizeof(double));
  memcpy(attr._zaxis, zaxis, 3 * sizeof(double));

  // Resample region labels image
  ByteImage output(attr, 1);
  RealImage pbmaps[NUM], input(regions.Attributes(), 1);
  GenericLinearInterpolateImageFunction<RealImage> pbfunc;
  pbfunc.Input(&input);
  pbfunc.Initialize();

  for (BytePixel c = 0; c < NUM; ++c) {
    for (int k = 0; k < regions.Z(); ++k)
    for (int j = 0; j < regions.Y(); ++j)
    for (int i = 0; i < regions.X(); ++i) {
      input(i, j, k) = (regions(i, j, k) == c ? 1. : 0.);
    }
    auto &pbmap = pbmaps[c];
    pbmap.Initialize(attr, 1);
    for (int k = 0; k < pbmap.Z(); ++k)
    for (int j = 0; j < pbmap.Y(); ++j)
    for (int i = 0; i < pbmap.X(); ++i) {
      p = Point(i, j, k);
      output.ImageToWorld(p);
      pbfunc.WorldToImage(p);
      pbmap(i, j, k) = pbfunc.Evaluate(p);
    }
  }

  RealPixel prob;
  for (int k = 0; k < output.Z(); ++k)
  for (int j = 0; j < output.Y(); ++j)
  for (int i = 0; i < output.X(); ++i) {
    prob = 0.;
    for (BytePixel c = 0; c < NUM; ++c) {
      const auto &pbmap = pbmaps[c];
      if (pbmap(i, j, k) > prob) {
        prob = pbmap(i, j, k);
        output(i, j, k) = c;
      }
    }
  }

  // Crop image
  output.PutBackgroundValueAsDouble(0.);
  output.BoundingBox(bi[0], bj[0], bk[0], bi[1], bj[1], bk[1]);

  bi[0] -= xmargin, bi[1] += xmargin;
  bj[0] -= ymargin, bj[1] += ymargin;
  bk[0] -= zmargin, bk[1] += zmargin;

  if (fgmask_name) {
    ByteImage fgmask(fgmask_name);
    for (int k = 0; k < fgmask.Z(); ++k)
    for (int j = 0; j < fgmask.Y(); ++j)
    for (int i = 0; i < fgmask.X(); ++i) {
      if (fgmask(i, j, k) != 0) {
        p = Point(i, j, k);
        fgmask.ImageToWorld(p);
        output.WorldToImage(p);
        bi[0] = min(bi[0], ifloor(p._x));
        bi[1] = max(bi[1], iceil (p._x));
        bj[0] = min(bj[0], ifloor(p._y));
        bj[1] = max(bj[1], iceil (p._y));
        bk[0] = min(bk[0], ifloor(p._z));
        bk[1] = max(bk[1], iceil (p._z));
      }
    }
  }

  if (bi[0] < 0) bi[0] = 0;
  if (bj[0] < 0) bj[0] = 0;
  if (bk[0] < 0) bk[0] = 0;

  if (bi[1] >= output.X()) bi[1] = output.X() - 1;
  if (bj[1] >= output.Y()) bj[1] = output.Y() - 1;
  if (bk[1] >= output.Z()) bk[1] = output.Z() - 1;

  output = output.GetRegion(bi[0], bj[0], bk[0], bi[1], bj[1], bk[1]);

  // Dilate BS 2 voxels into WM at BS/WM boundary to cut it cleanly again
  for (int k = output.Z()-1; k > 1; --k) {
    for (int j = 0; j < output.Y(); ++j)
    for (int i = 0; i < output.X(); ++i) {
      auto &region = output(i, j, k);
      if (region == LH || region == RH || region == UH) {
        if (output(i, j, k-1) == BS || output(i, j, k-2) == BS) {
          region = BS;
        }
      }
    }
  }

  // Cut WM and BS(+CB) with cutting planes / output image axes
  for (int k = 0; k < output.Z(); ++k)
  for (int j = 0; j < output.Y(); ++j)
  for (int i = 0; i < output.X(); ++i) {
    auto &region = output(i, j, k);
    if (region == BS || region == CB) {
      p = Point(i, j, k);
      output.ImageToWorld(p);
      if (bs_plane.SignedDistance(p) < 0.) {
        if (rl_plane.SignedDistance(p) < 0.) {
          region = LH;
        } else {
          region = RH;
        }
      }
    } else if (region == UH || region == RH || region == LH) {
      p = Point(i, j, k);
      output.ImageToWorld(p);
      if (rl_plane.SignedDistance(p) < 0.) {
        region = LH;
      } else {
        region = RH;
      }
    }
  }
  FixHemisphereLabels(output);

  return output;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  Point p;

  REQUIRES_POSARGS(0);

  const char *labels_name = nullptr; // Input  labels image
  const char *output_name = nullptr; // Output labels image

  const char *fgmask_name  = nullptr; // brain mask
  const char *hemis_name   = nullptr; // RH/LH cerebrum label image
  const char *rhmask_name  = nullptr; // RH cerebrum mask
  const char *lhmask_name  = nullptr; // LH cerebrum mask
  const char *wmmask_name  = nullptr; // Cerebral WM mask
  const char *gmmask_name  = nullptr; // Cerebral GM mask
  const char *sbmask_name  = nullptr; // Subcortical structures mask
  const char *bsmask_name  = nullptr; // Brainstem mask
  const char *cbmask_name  = nullptr; // Cerebellum mask
  const char *depth_name   = nullptr; // Output cortical depth map

  UnorderedSet<int> rhmask_labels;
  UnorderedSet<int> lhmask_labels;
  UnorderedSet<int> wmmask_labels;
  UnorderedSet<int> gmmask_labels;
  UnorderedSet<int> sbmask_labels;
  UnorderedSet<int> bsmask_labels;
  UnorderedSet<int> cbmask_labels;

  int  sb_closing  = 0;     // Subcortical segmentation closing iterations
  int  bs_closing  = 0;     // Brainstem   segmentation closing iterations
  int  cb_closing  = 0;     // Cerebellum  segmentation closing iterations
  int  xmargin     = 0;     // Margin in x direction after resampling in no. of voxels
  int  ymargin     = 0;     // Margin in y direction after resampling in no. of voxels
  int  zmargin     = 0;     // Margin in z direction after resampling in no. of voxels
  bool merge_bs_cb = false; // Whether to merge brainstem and cerebellum segments

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
    else if (OPTION("-output-inner-cortical-distance")) {
      depth_name = ARGUMENT;
    }
    else if (OPTION("-fg") || OPTION("-foreground") || OPTION("-brain")) {
      fgmask_name = ARGUMENT;
    }
    else if (OPTION("-margin")) {
      PARSE_ARGUMENT(xmargin);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(ymargin);
        if (HAS_ARGUMENT) {
          PARSE_ARGUMENT(zmargin);
        } else {
          zmargin = 0;
        }
      } else {
        ymargin = zmargin = xmargin;
      }
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
    else if (OPTION("-brainstem-and-cerebellum") || OPTION("-bscb") ||
             OPTION("-cerebellum-and-brainstem") || OPTION("-cbbs") ||
             OPTION("-brainstem+cerebellum") || OPTION("-bs+cb") ||
             OPTION("-cerebellum+brainstem") || OPTION("-cb+bs")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(merge_bs_cb);
      else merge_bs_cb = true;
    }
    else if (OPTION("-closing") || OPTION("-closing-iterations")) {
      PARSE_ARGUMENT(sb_closing);
      bs_closing = cb_closing = sb_closing;
    }
    else if (OPTION("-subcortical-closing")) {
      PARSE_ARGUMENT(sb_closing);
    }
    else if (OPTION("-brainstem-closing")) {
      PARSE_ARGUMENT(bs_closing);
    }
    else if (OPTION("-cerebellum-closing")) {
      PARSE_ARGUMENT(cb_closing);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Output labels image name required
  if (output_name == nullptr) {
    FatalError("Missing -output-labels file name!");
  }

  // ---------------------------------------------------------------------------
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

  const double ds = max(max(attr._dx, attr._dy), attr._dz);

  // ---------------------------------------------------------------------------
  // Initial brain regions image used to find cutting planes
  ByteImage regions(attr, 1);
  regions.PutBackgroundValueAsDouble(0.);
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

  if (debug) {
    regions.Write("debug_regions.nii.gz");
  }

  if (sb_closing > 0 || bs_closing > 0 || cb_closing > 0) {
    // Close holes in subcortical mask
    if (!sbmask.IsEmpty() && sb_closing > 0) {
      BinaryImage closed(attr, 1);
      for (int vox = 0; vox < nvox; ++vox) {
        closed(vox) = sbmask(vox);
      }
      Dilate<BinaryPixel>(&closed, sb_closing,     CONNECTIVITY_18);
      Erode <BinaryPixel>(&closed, sb_closing + 1, CONNECTIVITY_18);
      for (int vox = 0; vox < nvox; ++vox) {
        if (regions(vox) == BG && (closed(vox) != 0 || sbmask(vox) != 0)) {
          regions(vox) = UH;
        }
      }
      if (debug) {
        regions.Write("debug_regions+sbmask_closed.nii.gz");
      }
    }

    // Close holes in brainstem+cerebellum region
    if (!bsmask.IsEmpty() && (bs_closing > 0 || cb_closing > 0)) {
      int ncb = 0;
      BinaryImage closed(attr, 1);
      if (bs_closing > 0) {
        for (int vox = 0; vox < nvox; ++vox) {
          closed(vox) = BinaryPixel(regions(vox) == BS ? 1 : 0);
        }
        Dilate<BinaryPixel>(&closed, bs_closing,     CONNECTIVITY_18);
        Erode <BinaryPixel>(&closed, bs_closing + 1, CONNECTIVITY_18);
        for (int vox = 0; vox < nvox; ++vox) {
          const auto &region = regions(vox);
          if (region == BG) {
            if (closed(vox) != 0) {
              regions(vox) = BS;
            }
          } else if (region == CB) {
            closed(vox) = 1, ++ncb;
          } else if (region != BS) {
            closed(vox) = 0;
          }
        }
        if (debug) {
          regions.Write("debug_regions+brainstem_closed.nii.gz");
        }
      } else if (cb_closing > 0) {
        for (int vox = 0; vox < nvox; ++vox) {
          if (regions(vox) == CB) {
            closed(vox) = BinaryPixel(1), ++ncb;
          } else {
            closed(vox) = BinaryPixel(0);
          }
        }
      }
      if (ncb > 0 && cb_closing > 0) {
        Dilate<BinaryPixel>(&closed, cb_closing,     CONNECTIVITY_18);
        Erode <BinaryPixel>(&closed, cb_closing + 1, CONNECTIVITY_18);
        for (int vox = 0; vox < nvox; ++vox) {
          const auto &region = regions(vox);
          if (region == BG && closed(vox) != 0) {
            regions(vox) = CB;
          }
        }
        if (debug) {
          regions.Write("debug_regions+cerebellum_closed.nii.gz");
        }
      }
    }

    // Final closed regions image
    if (debug) {
      regions.Write("debug_regions+all_closed.nii.gz");
    }
  }

  // ---------------------------------------------------------------------------
  // Determine cutting planes
  Plane rl_plane, bs_plane;
  {
    rl_plane = MedialCuttingPlane(regions);
    Vector3 xaxis = rl_plane.n;
    Vector3 yaxis = xaxis.Cross(PrincipalDirection(regions, BS));
    Vector3 zaxis = xaxis.Cross(yaxis);
    bs_plane = BrainstemCuttingPlane(regions, zaxis);
  }
  if (verbose) {
    if (bs_plane) cout << "Brainstem cutting plane equation: " << bs_plane << endl;
    if (rl_plane) cout << "Medial cutting plane equation:    " << rl_plane << endl;
  }

  // ---------------------------------------------------------------------------
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
    for (int vox = 0; vox < nvox; ++vox) {
      if (gm(vox) == 1) {
        regions(vox) = GM;
      }
    }
  }

  // Cut subcortical structures and brainstem into RH, LH, and BS+CB
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
      }
    } else if (region == UH || (!sbmask.IsEmpty() && sbmask(i, j, k) != 0)) {
      if (rl_plane.SignedDistance(p) < 0.) {
        region = LH;
      } else {
        region = RH;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Resample regions mask such that image axes are aligned with determined
  // orthogonal cutting planes x: L -> R, y: P -> A, z: I -> S
  regions = Resample(regions, rl_plane, bs_plane, xmargin, ymargin, zmargin, fgmask_name);
  attr = regions.Attributes();
  nvox = regions.NumberOfSpatialVoxels();
  if (debug) {
    regions.Write("debug_output.nii.gz");
  }

  // ---------------------------------------------------------------------------
  {
    // Determine bounding box of interhemisphere WM
    const bool world = false;
    PointSet voxels = BoundaryPoints(regions, RH, LH, world);
    int bounds[6] = {static_cast<int>(voxels(0)._x), static_cast<int>(voxels(0)._x),
                     static_cast<int>(voxels(0)._y), static_cast<int>(voxels(0)._y),
                     static_cast<int>(voxels(0)._z), static_cast<int>(voxels(0)._z)};
    for (int n = 1, i, j, k; n < voxels.Size(); ++n) {
      const auto &voxel = voxels(n);
      i = static_cast<int>(voxel._x);
      j = static_cast<int>(voxel._y);
      k = static_cast<int>(voxel._z);
      bounds[0] = min(bounds[0], i);
      bounds[1] = max(bounds[1], i);
      bounds[2] = min(bounds[2], j);
      bounds[3] = max(bounds[3], j);
      bounds[4] = min(bounds[4], k);
      bounds[5] = max(bounds[5], k);
    }
    voxels.Clear();

    bounds[0] = max(0,               bounds[0] - 20);
    bounds[1] = min(regions.X() - 1, bounds[1] + 20);

    // -------------------------------------------------------------------------
    // 1. Fill small holes in xz slices nearby interhemispheric bounding box

    // Get first xz slice
    GreyImage holes = regions.GetRegion(bounds[0],     bounds[2],     bounds[4],
                                        bounds[1] + 1, bounds[2] + 1, bounds[5] + 1);

    // Mark holes in current xz slice, update y origin after each iteration
    NeighborhoodOffsets offsets(&regions, CONNECTIVITY_6);
    for (int j = bounds[2]; j <= bounds[3]; ++j) {

      // Extract non-WM mask in xz slice nearby RH/LH boundary
      for (int k = bounds[4], kk = 0; k <= bounds[5]; ++k, ++kk)
      for (int i = bounds[0], ii = 0; i <= bounds[1]; ++i, ++ii) {
        if (regions(i, j, k) != RH && regions(i, j, k) != LH) {
          holes(ii, 0, kk) = 1;
        } else {
          holes(ii, 0, kk) = 0;
        }
      }

      // Label holes
      ConnectedComponents<GreyPixel> cc;
      cc.Input (&holes);
      cc.Output(&holes);
      cc.Ordering(CC_SmallestFirst);
      cc.Run();

      // Fill small interhemisphere holes
      for (GreyPixel hid = 1; hid <= cc.NumberOfComponents(); ++hid) {

        // Hole may not exceed a preset size
        bool discard = (cc.ComponentSize(hid) > 20);

        // Hole may not extend beyond boundary of extracted mask
        if (!discard) {
          for (int ii = 0, kk = holes.Z()-1; ii < holes.X(); ++ii) {
            if (holes(ii, 0, 0) == hid || holes(ii, 0, kk) == hid) {
              discard = true;
              break;
            }
          }
        }
        if (!discard) {
          for (int ii = holes.X()-1, kk = 0; kk < holes.Z(); ++kk) {
            if (holes(0, 0, kk) == hid || holes(ii, 0, kk) == hid) {
              discard = true;
              break;
            }
          }
        }

        // Hole must touch both RH and LH segment
        if (!discard) {
          bool touches_rh = false;
          bool touches_lh = false;
          for (int kk = 1; kk < holes.Z()-1; ++kk)
          for (int ii = 1; ii < holes.X()-1; ++ii) {
            if (holes(ii, 0, kk) == hid) {
              const auto data = regions.Data(bounds[0] + ii, j, bounds[4] + kk);
              for (int n = 0; n < offsets.Size(); ++n) {
                const auto region = data + offsets(n);
                if (*region == RH) {
                  touches_rh = true;
                } else if (*region == LH) {
                  touches_lh = true;
                }
              }
            }
            if (touches_rh && touches_lh) break;
          }
          discard = (!touches_rh || !touches_lh);
        }

        // Fill hole if not to be discarded
        if (!discard) {
          int i, k;
          for (int kk = 0; kk < holes.Z(); ++kk)
          for (int ii = 0; ii < holes.X(); ++ii) {
            if (holes(ii, 0, kk) == hid) {
              i = bounds[0] + ii;
              k = bounds[4] + kk;
              p = Point(i, j, k);
              regions.ImageToWorld(p);
              if (rl_plane.SignedDistance(p) < 0.) {
                regions(i, j, k) = LH;
              } else {
                regions(i, j, k) = RH;
              }
            }
          }
        }
      }
    }

    // -------------------------------------------------------------------------
    // 2. Previous step may cut off some BG because of xz hole size threshold.
    //    Thus, fill in this BG with RH/LH, i.e., when it does not reach the
    //    boundary of the interhemispheric bounding box in xy slice

    // Get first xy slice
    holes = regions.GetRegion(bounds[0],     bounds[2],     bounds[4],
                              bounds[1] + 1, bounds[3] + 1, bounds[5] + 1);

    // Extract non-WM mask in xy slice nearby RH/LH boundary
    for (int k = bounds[4], kk = 0; k <= bounds[5]; ++k, ++kk)
    for (int j = bounds[2], jj = 0; j <= bounds[3]; ++j, ++jj)
    for (int i = bounds[0], ii = 0; i <= bounds[1]; ++i, ++ii) {
      if (regions(i, j, k) != RH && regions(i, j, k) != LH) {
        holes(ii, jj, kk) = 1;
      } else {
        holes(ii, jj, kk) = 0;
      }
    }

    // Label holes
    ConnectedComponents<GreyPixel> cc;
    cc.Input (&holes);
    cc.Output(&holes);
    cc.Ordering(CC_SmallestFirst);
    cc.Run();

    // Fill small interhemisphere holes
    for (GreyPixel hid = 1; hid <= cc.NumberOfComponents(); ++hid) {

      bool discard = false;

      // Hole may not extend beyond boundary of extracted mas
      if (!discard) {
        const int ii = holes.X() - 1;
        for (int kk = 0; kk < holes.Z(); ++kk)
        for (int jj = 0; jj < holes.Y(); ++jj) {
          if (holes(0, jj, kk) == hid || holes(ii, jj, kk) == hid) {
            discard = true;
            break;
          }
        }
      }
      if (!discard) {
        const int jj = holes.Y() - 1;
        for (int kk = 0; kk < holes.Z(); ++kk)
        for (int ii = 0; ii < holes.X(); ++ii) {
          if (holes(ii, 0, kk) == hid || holes(ii, jj, kk) == hid) {
            discard = true;
            break;
          }
        }
      }
      if (!discard) {
        const int kk = holes.Z() - 1;
        for (int jj = 0; jj < holes.Y(); ++jj)
        for (int ii = 0; ii < holes.X(); ++ii) {
          if (holes(ii, jj, 0) == hid || holes(ii, jj, kk) == hid) {
            discard = true;
            break;
          }
        }
      }

      // Fill hole if not to be discarded
      if (!discard) {
        int i, j, k;
        for (int kk = 0; kk < holes.Z(); ++kk)
        for (int jj = 0; jj < holes.Y(); ++jj)
        for (int ii = 0; ii < holes.X(); ++ii) {
          if (holes(ii, jj, kk) == hid) {
            i = bounds[0] + ii;
            j = bounds[2] + jj;
            k = bounds[4] + kk;
            p = Point(i, j, k);
            regions.ImageToWorld(p);
            if (rl_plane.SignedDistance(p) < 0.) {
              regions(i, j, k) = LH;
            } else {
              regions(i, j, k) = RH;
            }
          }
        }
      }
    }

    if (debug) {
      regions.Write("debug_output+rl_closed.nii.gz");
    }
  }

  // ---------------------------------------------------------------------------
  // Change CB(+BS) labels to GM(/BG) when bordering with WM away from the BS
  // cutting plane to ensure that there is some separating cortex between
  // white surface and BS+CB region which may be misclassified as cerebellar GM
  if (gmmask_name || !gmmask_labels.empty()) {
    const double min_dist_bs = 10. * ds;
    const double min_dist_cb =  2. * ds;
    bool is_empty = true;
    BinaryImage boundary(attr, 1);
    NeighborhoodOffsets offsets(&regions, CONNECTIVITY_18);
    for (int k = 1; k < regions.Z()-1; ++k)
    for (int j = 1; j < regions.Y()-1; ++j)
    for (int i = 1; i < regions.X()-1; ++i) {
      const auto &region = regions(i, j, k);
      if (region == CB || region == BS) {
        p = Point(i, j, k);
        regions.ImageToWorld(p);
        if (abs(bs_plane.SignedDistance(p)) > (region == CB ? min_dist_cb : min_dist_bs)) {
          const auto data = regions.Data(i, j, k);
          for (int n = 0; n < offsets.Size(); ++n) {
            const auto &label = *(data + offsets(n));
            if (label == LH || label == RH) {
              boundary(i, j, k) = 1;
              is_empty = false;
              break;
            }
          }
        }
      } else if (region == GM) {
        bool next_to_wm = false;
        bool next_to_cb = false;
        const auto data = regions.Data(i, j, k);
        for (int n = 0; n < offsets.Size(); ++n) {
          const auto &label = *(data + offsets(n));
          if (label == LH || label == RH) {
            next_to_wm = true;
          } else if (label == CB || label == BS) {
            next_to_cb = true;
          }
        }
        if (next_to_wm && next_to_cb) {
          boundary(i, j, k) = 1;
        }
      }
    }
    if (!is_empty) {
      Dilate<BinaryPixel>(&boundary, 1, CONNECTIVITY_18);
      for (int vox = 0; vox < nvox; ++vox) {
        auto &region = regions(vox);
        if (boundary(vox) != 0 && (region == CB || region == BS || region == BG)) {
          region = GM;
        }
      }
      if (debug) {
        regions.Write("debug_output+wm-cb_gap.nii.gz");
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Change WM labels to BG when bordering with BS below the cutting plane to
  // ensure at least one voxel gap between white surface and brainstem
  // (+cerebellum) surface mesh for merge-surfaces to be able to insert a
  // single cutting plane as divider.
  if (gmmask_name || !gmmask_labels.empty()) {
    const double ds = max(max(regions.XSize(), regions.YSize()), regions.ZSize());
    const double min_dist_bs = 2. * ds;
    bool is_empty = true;
    BinaryImage boundary(attr, 1);
    NeighborhoodOffsets offsets(&regions, CONNECTIVITY_18);
    for (int k = 1; k < regions.Z()-1; ++k)
    for (int j = 1; j < regions.Y()-1; ++j)
    for (int i = 1; i < regions.X()-1; ++i) {
      const auto &region = regions(i, j, k);
      if (region == BS) {
        p = Point(i, j, k);
        regions.ImageToWorld(p);
        if (abs(bs_plane.SignedDistance(p)) > min_dist_bs) {
          const auto data = regions.Data(i, j, k);
          for (int n = 0; n < offsets.Size(); ++n) {
            const auto label = data + offsets(n);
            if (*label == LH || *label == RH) {
              const auto vox = static_cast<int>(label - regions.Data());
              boundary(vox) = 1;
              is_empty = false;
            }
          }
        }
      }
    }
    if (!is_empty) {
      Dilate<BinaryPixel>(&boundary, 1, CONNECTIVITY_18);
      for (int vox = 0; vox < nvox; ++vox) {
        auto &region = regions(vox);
        if (boundary(vox) != 0 && (region == LH || region == RH)) {
          region = BG;
        }
      }
      if (debug) {
        regions.Write("debug_output+wm-bs_gap.nii.gz");
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Merge BS and CB labels
  if (merge_bs_cb) {
    const auto old_label = max(BS, CB);
    const auto new_label = min(BS, CB);
    for (int vox = 0; vox < nvox; ++vox) {
      auto &region = regions(vox);
      if (region == old_label) region = new_label;
    }
  }

  // ---------------------------------------------------------------------------
  // Write output labels
  regions.Write(output_name);

  // ---------------------------------------------------------------------------
  // Create distance map for interior of cortical surface
  if (depth_name) {

    // Determine yz plane next to mid-plane cut
    int bi;
    for (bi = 0; bi < attr._x; ++bi) {
      p = Point(bi, 0, 0);
      regions.ImageToWorld(p);
      if (rl_plane.SignedDistance(p) > 0.) break;
    }
    if (bi > 0) --bi;

    // Initialize mask of cortical surface interior
    BinaryImage mask(attr, 1);
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      if (regions(i, j, k) == RH || regions(i, j, k) == LH || regions(i, j, k) == GM) {
        mask(i, j, k) = 1;
      } else {
        mask(i, j, k) = 0;
      }
    }

    // Close small holes of BG labeled voxels
    Close<BinaryPixel>(&mask, 5, CONNECTIVITY_18);

    // Ensure a clear separation of the two hemispheres outside the subcortical
    // structures such as Corpus Callosum even if it cuts through cortex
    for (int k =  0; k <  attr._z; ++k)
    for (int j =  0; j <  attr._y; ++j)
    for (int i = bi; i <= bi + 1;  ++i) {
      if (regions(i, j, k) != RH && regions(i, j, k) != LH) {
        mask(i, j, k) = 0;
      }
    }
    if (debug) {
      mask.Write("debug_interior_mask.nii.gz");
    }

    // Compute distance transform
    RealImage input(attr, 1), output(attr, 1);
    for (int vox = 0; vox < nvox; ++vox) {
      input(vox) = (mask(vox) != 0 ? 0. : 1.);
    }
    EuclideanDistanceTransform<RealPixel> filter;
    filter.Input (&input);
    filter.Output(&output);
    filter.Run();
    for (int vox = 0; vox < nvox; ++vox) {
      output(vox) = sqrt(output(vox));
    }

    output.Write(depth_name);
  }

  return 0;
}
