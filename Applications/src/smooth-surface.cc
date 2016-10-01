/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Area weighted Laplacian smoothing proposed by Tosun in MedIA 2004
 * implemented by Paul Aljabar.
 *
 * Copyright Paul Aljabar
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

#include "mirtk/EdgeTable.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/SurfaceCurvature.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkIdList.h"
#include "vtkPolyDataNormals.h"
#include "vtkCurvatures.h"
#include "vtkDataSetAttributes.h"
#include "vtkWindowedSincPolyDataFilter.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << "       " << name << " <input> <output> <iterations> [lambda] [mu] [options] (deprecated)" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Smooths the node positions and/or scalar data of a surface mesh." << endl;
  cout << endl;
  cout << "  When node positions are smoothed, iteratively move mesh nodes towards the" << endl;
  cout << "  centroid of the adjacent nodes. The relaxation factor (RF) determines how" << endl;
  cout << "  far each node moves towards the local centroid (0 <= lambda <= 1). The new position" << endl;
  cout << "  of a node is a weighted combination of its previous position and the neighbours' centroid." << endl;
  cout << "  If lambda = 1, a node is moved all the way to the centroid, while lambda = 0 keeps a node" << endl;
  cout << "  at its previous position." << endl;
  cout << endl;
  cout << "  For a low-pass filtering of the surface mesh, set :option:`-lambda` to a positive value in" << endl;
  cout << "  the range (0, 1), and :option:`-mu` to a negative value greater in magnitude than" << endl;
  cout << "  lambda, i.e., 0 < lambda < mu. For a detailed explanation of the parameters," << endl;
  cout << "  cf. Taubin, Curve and Surface Smoothing without Shrinkage, ICCV 1995." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input        File name of input surface mesh." << endl;
  cout << "  output       File name of output surface mesh." << endl;
  cout << "  iterations   No. of smoothing iterations. (deprecated: see :option:`-iterations` option)" << endl;
  cout << "  lambda       Relaxation/scale factor for odd  iterations. (deprecated: see -lambda option)" << endl;
  cout << "  mu           Magnitude of scale factor for even iterations. (deprecated: see :option:`-mu` option)" << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -iterations <n>          Number of smoothing iterations. (default: 1)" << endl;
  cout << "  -lambda <float>          Relaxation/scale factor for odd  iterations. (default: 1)" << endl;
  cout << "  -mu <float>              Relaxation/scale factor for even iterations. (default: lambda)" << endl;
  cout << "  -points                  Smooth node positions. (default if no :option:`-scalars` specified)" << endl;
  cout << "  -scalars [<name>]        Name of scalar point data array to smooth. If no <name> given," << endl;
  cout << "                           all non-attribute scalar arrays and the SCALARS if any are smoothed." << endl;
  cout << "  -exclnode                Exclude node itself from smoothing, only average adjacent values." << endl;
  cout << "  -inclnode                Include node itself in the smoothing." << endl;
  cout << "  -areaweighted            Area weighted Laplacian relaxation of surface node positions. (default)" << endl;
  cout << "  -windowedsinc [<band>]   Adjust point positions using a windowed sinc interpolation function." << endl;
  cout << "  -combinatorial           Combinatorial weighting of adjacent nodes." << endl;
  cout << "  -distance [<sigma>]      Inverse distance weighted smoothing kernel." << endl;
  cout << "  -gaussian [<sigma>]      Gaussian smoothing kernel." << endl;
  cout << endl;
  cout << "  -anisotropic ([<sigma>] [<tensor_array>] | <sigma1> <sigma2> [(<e2_array> | <e1_array> <e2_array>)])" << endl;
  cout << "      Anisotropic Gaussian smoothing kernel given input tensor field point data array" << endl;
  cout << "      named <tensor_array>. The default tensor field used with no arguments or only <sigma> specified" << endl;
  cout << "      is the curvature tensor output of the calculate-surface-attributes -tensor command. Alternatively," << endl;
  cout << "      specify <sigma1> and <sigma2> together with the direction of maximum curvature in <e2_array> or both," << endl;
  cout << "      the direction of minimum and maximum curvature direction in point data arrays named <e1_array> and <e2_array>." << endl;
  cout << endl;
  cout << "  -track [<name>]          Track the signed distances traversed by points (positive is outwards)." << endl;
  cout << "                           The optional arguments specifies the name of output output point data array." << endl;
  cout << "                           (default: \"smoothingDists\" if option is given without <name>)." << endl;
  cout << "  -threshold <value>       A value indicating the smoothness as the L_2 norm of H^2 over the surface." << endl;
  cout << "                           See Tosun, MedIA, 2004. Iterations stop if the norm for the surface drops" << endl;
  cout << "                           below the given value. Only used for :option:`-areaweighted` smoothing." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Tosun, MedIA, 2004
// =============================================================================

// -----------------------------------------------------------------------------
double hSquareRobustMean(vtkPolyData* surface)
{
  const vtkIdType noOfPoints = surface->GetNumberOfPoints();

  const double lo =  5.0;
  const double hi = 95.0;

  // Compute mean curvature
  vtkCurvatures *curve = vtkCurvatures::New();
  SetVTKInput(curve, surface);
  curve->SetCurvatureTypeToMean();
  curve->Update();

  vtkPolyData  *output  = curve->GetOutput();
  vtkDataArray *scalars = output->GetPointData()->GetScalars("Mean_Curvature");

  // Sort mean curvature values
  Array<double> data(noOfPoints);
  for (vtkIdType i = 0; i < noOfPoints; ++i) data[i] = scalars->GetTuple1(i);
  sort(data.begin(), data.end());

  // Get 5th and 95th percentile
  double minVal = data[int(round(lo * (noOfPoints - 1) / 100.0))];
  double maxVal = data[int(round(hi * (noOfPoints - 1) / 100.0))];

  // Compute robust squared mean within percentile range
  double sum = .0, val;
  int    num = 0;

  for (vtkIdType i = 0; i < noOfPoints; ++i) {
    val = scalars->GetTuple1(i);
    if (minVal <= val && val <= maxVal) {
      sum += val * val;
      ++num;
    }
  }

  return sum / num;
}

// -----------------------------------------------------------------------------
double surfaceArea(vtkPolyData* surface)
{
  vtkCellArray* facets = surface->GetPolys();
  vtkTriangle* facet = vtkTriangle::New();

  double A = 0.0;
  double v0[3], v1[3], v2[3];

  vtkIdType f, *vert=0;
  facets->InitTraversal();
  while (facets->GetNextCell(f,vert)){

    surface->GetPoint(vert[0],v0);
    surface->GetPoint(vert[1],v1);
    surface->GetPoint(vert[2],v2);

    A += double(facet->TriangleArea(v0,v1,v2));
  }

  return A;
}

// -----------------------------------------------------------------------------
void getCoG(vtkPolyData *input, double*cog)
{
  int i, noOfPoints;
  double cofgx, cofgy, cofgz;
  double point[3];

  cofgx = cofgy = cofgz = 0.0;

  noOfPoints = input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    input->GetPoint (i, point);
    cofgx += point[0];
    cofgy += point[1];
    cofgz += point[2];
  }

  if (noOfPoints > 0){
    cofgx /= noOfPoints;
    cofgy /= noOfPoints;
    cofgz /= noOfPoints;
  }

  cog[0] = cofgx;
  cog[1] = cofgy;
  cog[2] = cofgz;
}

// -----------------------------------------------------------------------------
void shiftAndScalePolyData(vtkPolyData* input, double *shift, double factor)
{
  int i, j, noOfPoints;
  noOfPoints = input->GetNumberOfPoints();
  double vOld[3], vNew[3];

  for (i = 0; i < noOfPoints; ++i){
    input->GetPoint(i, vOld);

    for (j = 0; j < 3; ++j){
      vNew[j] = factor * (vOld[j] + shift[j]);
    }

    input->GetPoints()->SetPoint(i, vNew);
  }
}

// -----------------------------------------------------------------------------
double meanRadius(vtkPolyData* input, double*cog)
{
  int i, noOfPoints;
  double point[3];
  double r, rSum = 0.0;

  noOfPoints = input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    input->GetPoint (i, point);
    r = sqrt((point[0]-cog[0])*(point[0]-cog[0]) +
             (point[1]-cog[1])*(point[1]-cog[1]) +
             (point[2]-cog[2])*(point[2]-cog[2]));
    rSum += r;
  }

  return rSum / noOfPoints;
}

// -----------------------------------------------------------------------------
void AreaWeightedLaplacianSmoothing(vtkPolyData *input, vtkDataArray *mask,
                                    int    noOfIterations,
                                    double lambda, double mu,
                                    double smoothnessThreshold = -1.0,
                                    bool   trackingOn = false)
{
  int i, j, k;
  double E_H2, area;

  double currPos[3];
  unsigned short noOfCells = 0;
  vtkIdType* cells = NULL;
  vtkTriangle* triangle = NULL;
  double totalArea = 0;
  double update[3];
  vtkIdList* ptIds = NULL;
  double v1[3], v2[3], v3[3], centre[3];
  double triangleArea = 0;
  double dx, dy, dz;
  double dist, val;
  double *normal;
  double h2norm;

  double cogOld[3];
  double cogNew[3];
  double radiusOld, radiusNew;
  double shift[3];
  double scaleFactor;

  const int noOfPoints = static_cast<int>(input->GetNumberOfPoints());
  vtkDataArray * const normals = input->GetPointData()->GetNormals();

  // Output points
  Array<double> pts(3 * noOfPoints);
  for (j = 0; j < noOfPoints; ++j) {
    input->GetPoint(j, pts.data() + 3 * j);
  }

  // Sum of signed distances traveled by a point
  vtkSmartPointer<vtkDataArray> dists;
  if (trackingOn) {
    dists = vtkSmartPointer<vtkFloatArray>::New();
    dists->SetName("smoothingDists");
    dists->SetNumberOfComponents(1);
    dists->SetNumberOfTuples(noOfPoints);
    dists->FillComponent(0, 0.);
    input->GetPointData()->AddArray(dists);
  }

  // Smoothing iterations
  for (i = 0; i <= noOfIterations; ++i) {
    if (verbose) cout << "iteration  " << i << " ";

    double relaxationFactor = lambda;
    if (i % 2 == 1) relaxationFactor = mu;

    // Estimate \int H^2 dA by multiplying E(H^2) with Area.
    E_H2 = hSquareRobustMean(input);
    area = surfaceArea(input);

    // The L_2 norm using the Tosun formulation (MedIA 2004)
    h2norm = sqrt(E_H2 * area / 4.0 / pi);
    if (h2norm < smoothnessThreshold) break;
    if (verbose > 1) cout << h2norm << endl;

    getCoG(input, cogOld);
    radiusOld = meanRadius(input, cogOld);

    // Loop over surface.
    for (j = 0; j < noOfPoints; ++j) {

      // Skip excluded points
      if (mask && mask->GetComponent(j, 0) == 0.) continue;

      // What cells does this node adjoin?
      input->GetPointCells(j, noOfCells, cells);
      if (noOfCells == 0) continue;

      // Store the current position of the node.
      input->GetPoint(j, currPos);

      // Initialisation for current point.
      totalArea = 0;
      update[0] = 0;
      update[1] = 0;
      update[2] = 0;

      for (k = 0; k < noOfCells; ++k) {
        triangle = vtkTriangle::SafeDownCast(input->GetCell(cells[k]));
        if (triangle == nullptr) continue;
        ptIds = triangle->GetPointIds();

        input->GetPoint(ptIds->GetId(0), v1);
        input->GetPoint(ptIds->GetId(1), v2);
        input->GetPoint(ptIds->GetId(2), v3);

        triangleArea = vtkTriangle::TriangleArea(v1, v2, v3);
        vtkTriangle::TriangleCenter(v1, v2, v3, centre);

        totalArea += triangleArea;

        update[0] += triangleArea * centre[0];
        update[1] += triangleArea * centre[1];
        update[2] += triangleArea * centre[2];
      }

      if (totalArea <= 0.) {
        update[0] = currPos[0];
        update[1] = currPos[1];
        update[2] = currPos[2];
      } else {
      	update[0] /= totalArea;
      	update[1] /= totalArea;
      	update[2] /= totalArea;
      }

      dx = relaxationFactor * (update[0] - currPos[0]);
      dy = relaxationFactor * (update[1] - currPos[1]);
      dz = relaxationFactor * (update[2] - currPos[2]);

      pts[j*3  ] = currPos[0] + dx;
      pts[j*3+1] = currPos[1] + dy;
      pts[j*3+2] = currPos[2] + dz;

      if (dists) {
        dist = sqrt(dx*dx + dy*dy + dz*dz);
        normal = normals->GetTuple3(j);
        val = normal[0]*dx + normal[1]*dy + normal[2]*dz;
        if (val < 0) dist = -dist;
        dists->SetComponent(j, 0, dists->GetComponent(j, 0) + dist);
      }
    }

    for (j = 0; j < noOfPoints; ++j) {
      input->GetPoints()->SetPoint(j, pts.data() + 3*j);
    }

    // update radius and centre of gravity
    getCoG(input, cogNew);
    radiusNew = meanRadius(input, cogNew);

    shift[0] = cogOld[0] - cogNew[0];
    shift[1] = cogOld[1] - cogNew[1];
    shift[2] = cogOld[2] - cogNew[2];

    scaleFactor = radiusOld / radiusNew;

    shiftAndScalePolyData(input, shift, scaleFactor);
    if (verbose) cout << endl;
  }

  if (verbose) {
    cout << "Final iterations : " << i << endl;
    cout << "Final L_2 norm of H^2 (threshold) : " << h2norm << " (" << smoothnessThreshold << ")" << endl;
  }
}

// =============================================================================
// vtkWindowedSincPolyDataFilter
// =============================================================================

// -----------------------------------------------------------------------------
/// Compute centroid of point set
void GetCentroid(vtkPolyData *mesh, double centroid[3])
{
  double p[3];
  centroid[0] = centroid[1] = centroid[2] = .0;
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    mesh->GetPoint(ptId, p);
    centroid[0] += p[0];
    centroid[1] += p[1];
    centroid[2] += p[2];
  }
  centroid[0] /= npoints;
  centroid[1] /= npoints;
  centroid[2] /= npoints;
}

// -----------------------------------------------------------------------------
/// Get approximate scale of point set
void GetScale(vtkPolyData *mesh, const double centroid[3], double scale[3])
{
  double p[3];
  scale[0] = scale[1] = scale[2] = .0;
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    mesh->GetPoint(ptId, p);
    scale[0] += abs(p[0] - centroid[0]);
    scale[1] += abs(p[1] - centroid[1]);
    scale[2] += abs(p[2] - centroid[2]);
  }
  scale[0] /= npoints;
  scale[1] /= npoints;
  scale[2] /= npoints;
}

// -----------------------------------------------------------------------------
// Perform low-pass filtering using windowed sinc interpolation function
//
// The translation and scale "fix" is due to a bug in vtkWindowedSincPolyDataFilter:
// http://vtk.1045678.n5.nabble.com/Bug-in-vtkWindowedSincPolyDataFilter-td1234055.html
void WindowedSincSmoothing(vtkSmartPointer<vtkPolyData> surface, vtkDataArray *mask, int niter, double band)
{
  double c1[3], s1[3], c2[3], s2[3], p1[3], p2[3];

  GetCentroid(surface, c1);
  GetScale(surface, c1, s1);

  vtkNew<vtkWindowedSincPolyDataFilter> filter;
  filter->SetPassBand(band);
  filter->SetNumberOfIterations(niter);
  filter->NormalizeCoordinatesOn();
  SetVTKInput(filter, surface);
  filter->Update();

  GetCentroid(filter->GetOutput(), c2);
  GetScale(filter->GetOutput(), c2, s2);

  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    if (!mask || mask->GetComponent(ptId, 0) != 0.) {
      surface->GetPoint(ptId, p1);
      filter->GetOutput()->GetPoint(ptId, p2);
      p2[0] = c1[0] + s1[0] * (p2[0] - c2[0]) / s2[0];
      p2[1] = c1[1] + s1[1] * (p2[1] - c2[1]) / s2[1];
      p2[2] = c1[2] + s1[2] * (p2[2] - c2[2]) / s2[2];
      surface->GetPoints()->SetPoint(ptId, p2);
    }
  }
}

// =============================================================================
// Cortical surface smoothing
// =============================================================================

// -----------------------------------------------------------------------------
/// Smooth gyral points along maximum curvature direction
void SmoothInMaximumCurvatureDirection(vtkSmartPointer<vtkPolyData> surface,
                                       vtkDataArray *mask = nullptr,
                                       double lambda = 1.0,
                                       bool adjacent_only = true)
{
  vtkDataArray *k2_array = surface->GetPointData()->GetArray("Maximum_Curvature");
  vtkDataArray *e2_array = surface->GetPointData()->GetArray("Maximum_Curvature_Direction");

  if (!k2_array) {
    FatalError("Input surface mesh has no point data array named Maximum_Curvature");
  }
  if (!e2_array) {
    FatalError("Input surface mesh has no point data array named Maximum_Curvature_Direction");
  }

  double k2_range[2];
  k2_array->GetRange(k2_range);

  EdgeTable  edgeTable(surface);
  int        numAdjPts;
  const int *adjPtIds;
  double     alpha, beta;

  vtkSmartPointer<vtkPoints> points;
  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(surface->GetNumberOfPoints());

  double k2, e2[3], p1[3], p2[3], p[3], e[3], d, w, wsum;
  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    surface->GetPoint(ptId, p1);
    if (!mask || mask->GetComponent(ptId, 0) != 0.) {
      k2 = k2_array->GetComponent(ptId, 0);
      if (k2 > .0) {
        e2_array->GetTuple(ptId, e2);
        vtkMath::Normalize(e2);
        if (adjacent_only) {
          p[0] = p[1] = p[2] = wsum = .0;
        } else {
          p[0] = p1[0], p[1] = p1[1], p[2] = p1[2];
          wsum = 1.0;
        }
        edgeTable.GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        for (int i = 0; i < numAdjPts; ++i) {
          surface->GetPoint(adjPtIds[i], p2);
          vtkMath::Subtract(p2, p1, e);
          d = vtkMath::Norm(e);
          vtkMath::MultiplyScalar(e, 1.0 / d);
          w = abs(vtkMath::Dot(e, e2));
          p[0] += w * p2[0];
          p[1] += w * p2[1];
          p[2] += w * p2[2];
          wsum += w;
        }
        if (wsum > .0) {
          alpha = lambda * (k2 / k2_range[1]);
          beta  = alpha / wsum;
          alpha = 1.0 - alpha;
          p1[0] = alpha * p1[0] + beta * p[0];
          p1[1] = alpha * p1[1] + beta * p[1];
          p1[2] = alpha * p1[2] + beta * p[2];
        }
      }
    }
    points->SetPoint(ptId, p1);
  }

  surface->SetPoints(points);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
enum WeightFunction {
  Default,
  Combinatorial       = MeshSmoothing::Combinatorial,
  InverseDistance     = MeshSmoothing::InverseDistance,
  Gaussian            = MeshSmoothing::Gaussian,
  AnisotropicGaussian = MeshSmoothing::AnisotropicGaussian,
  AreaWeighted,
  WindowedSinc,
  GyralWeights
};

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  int noOfIterations = 1;
  if (NUM_POSARGS > 2) {
    if (!FromString(POSARG(3), noOfIterations)) {
      FatalError("Invalid no. of iterations argument: " << POSARG(3));
    }
  }

  double lambda = 1.0;
  if (NUM_POSARGS > 3) {
    if (!FromString(POSARG(4), lambda)) {
      FatalError("Invalid lambda value: " << POSARG(4));
    }
  }

  double mu = numeric_limits<double>::quiet_NaN();
  if (NUM_POSARGS > 4) {
    if (!FromString(POSARG(5), mu)) {
      FatalError("Invalid mu value: " << POSARG(5));
    }
    mu = -abs(mu);
  }

  if (NUM_POSARGS > 5) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Read the input mesh
  vtkSmartPointer<vtkPolyData> polydata = ReadPolyData(input_name);
  vtkSmartPointer<vtkDataArray> mask;

  // Parse remaining arguments
  Array<string> scalar_names;
  WeightFunction weighting   = Default;
  const char *tensor_name    = nullptr;
  const char *e1_name        = nullptr;
  const char *e2_name        = nullptr;
  bool   smooth_points       = false;
  double smoothnessThreshold = -1.0;
  double sigma1              = .0;
  double sigma2              = .0;
  bool   trackingOn          = false;
  bool   adjacent_only       = false;

  for (ALL_OPTIONS) {
    if (OPTION("-mask")) {
      const char *mask_name = ARGUMENT;
      mask = GetArrayByCaseInsensitiveName(polydata->GetPointData(), mask_name);
      if (!mask) {
        FatalError("Surface has no point data array named: " << mask_name);
      }
    }
    else if (OPTION("-nomask")) mask = nullptr;
    else if (OPTION("-points")) smooth_points = true;
    else if (OPTION("-scalars")) {
      if (HAS_ARGUMENT) {
        do {
          scalar_names.push_back(ARGUMENT);
        } while (HAS_ARGUMENT);
      } else {
        for (int i = 0; i < polydata->GetPointData()->GetNumberOfArrays(); ++i) {
          int type = polydata->GetPointData()->IsArrayAnAttribute(i);
          if (polydata->GetPointData()->GetArrayName(i) &&
              type != vtkDataSetAttributes::NORMALS &&
              type != vtkDataSetAttributes::TCOORDS &&
              type != vtkDataSetAttributes::PEDIGREEIDS &&
              type != vtkDataSetAttributes::GLOBALIDS &&
              type != vtkDataSetAttributes::EDGEFLAG) {
            scalar_names.push_back(polydata->GetPointData()->GetArrayName(i));
          }
        }
      }
    }
    else if (OPTION("-iterations")) PARSE_ARGUMENT(noOfIterations);
    else if (OPTION("-lambda"))     PARSE_ARGUMENT(lambda);
    else if (OPTION("-mu"))         PARSE_ARGUMENT(mu);
    else if (OPTION("-combinatorial") || OPTION("-umbrella")) {
      weighting = Combinatorial;
    }
    else if (OPTION("-distance") || OPTION("inversedistance")) {
      weighting = InverseDistance;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(sigma1);
      else sigma1 = .0;
    }
    else if (OPTION("-gaussian")) {
      weighting = Gaussian;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(sigma1);
      else sigma1 = .0;
    }
    else if (OPTION("-gyri")) weighting = GyralWeights;
    // -anisotropic
    // -anisotropic <sigma>
    // -anisotropic <sigma> <tensor_name>
    // -anisotropic <sigma1> <sigma2>
    // -anisotropic <sigma1> <sigma2> <e2_name>
    // -anisotropic <sigma1> <sigma2> <e1_name> <e2_name>
    else if (OPTION("-anisotropic")) {
      weighting = AnisotropicGaussian;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(sigma1);
      else sigma1 = .0;
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        if (FromString(arg, sigma2)) {
          if (HAS_ARGUMENT) {
            e2_name = ARGUMENT;
            if (HAS_ARGUMENT) {
              e1_name = e2_name;
              e2_name = ARGUMENT;
            } else {
              if (polydata->GetPointData()->HasArray(SurfaceCurvature::MINIMUM_DIRECTION)) {
                e1_name = SurfaceCurvature::MINIMUM_DIRECTION;
              }
            }
          } else {
            if (polydata->GetPointData()->HasArray(SurfaceCurvature::MINIMUM_DIRECTION)) {
              e1_name = SurfaceCurvature::MINIMUM_DIRECTION;
            }
            if (polydata->GetPointData()->HasArray(SurfaceCurvature::MAXIMUM_DIRECTION)) {
              e2_name = SurfaceCurvature::MAXIMUM_DIRECTION;
            }
          }
        } else {
          sigma2      = sigma1;
          tensor_name = ARGUMENT;
        }
      } else {
        if (polydata->GetPointData()->HasArray(SurfaceCurvature::TENSOR)) {
          tensor_name = SurfaceCurvature::TENSOR;
        }
      }
    }
    else if (OPTION("-exclnode") || OPTION("-adjacent")) adjacent_only = true;
    else if (OPTION("-inclnode"))                        adjacent_only = false;
    else if (OPTION("-areaweighted")) weighting = AreaWeighted;
    else if (OPTION("-windowedsinc")) {
      weighting = WindowedSinc;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(lambda);
    }
    else if (OPTION("-track")) trackingOn = true;
    else if (OPTION("-threshold")) PARSE_ARGUMENT(smoothnessThreshold);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (IsNaN(mu)) mu = lambda;

  if (weighting == Default) {
    if (scalar_names.empty()) {
      weighting = AreaWeighted;
    } else {
      weighting = Gaussian;
    }
  }

  if (!smooth_points && scalar_names.empty()) smooth_points = true;

  if (!scalar_names.empty() && (weighting == AreaWeighted ||
                                weighting == WindowedSinc ||
                                weighting == GyralWeights)) {
    FatalError("Cannot smooth scalar point data using this smoothing function");
  }

  const bool output_normals = (polydata->GetPointData()->GetNormals() != NULL);

  if (verbose) {
    cout << "Input:          " << input_name << endl;
    cout << "Mask:           " << (mask ? mask->GetName() : "None") << endl;
    cout << "Output:         " << output_name << endl;
    cout << "Iterations:     " << noOfIterations << endl;
    cout << "Lambda / mu:    " << lambda << " / " << mu << endl;
    cout << "Smooth points:  " << ToString(smooth_points) << endl;
    cout << "Smooth scalars:";
    if (scalar_names.empty()) cout << " No";
    for (size_t i = 0; i < scalar_names.size(); ++i) cout << " " << scalar_names[i];
    cout << endl;
  }

  // Normals are required by Tosun's method
  if (weighting == AreaWeighted) {
    vtkNew<vtkPolyDataNormals> filter;
    filter->SplittingOff();
    SetVTKInput(filter, polydata);
    filter->Update();
    polydata = filter->GetOutput();
  }

  // Smooth gyral points
  if (weighting == GyralWeights) {

    for (int iter = 0; iter < noOfIterations; ++iter) {
      if (iter % 2 == 0 || IsNaN(mu)) {
        SmoothInMaximumCurvatureDirection(polydata, mask, lambda, adjacent_only);
      } else {
        SmoothInMaximumCurvatureDirection(polydata, mask, mu, adjacent_only);
      }
    }

  }
  // Smooth node positions and/or scalar data
  else if (weighting == AreaWeighted) {

    AreaWeightedLaplacianSmoothing(polydata, mask, noOfIterations, lambda, mu, smoothnessThreshold, trackingOn);

  }
  // Smooth node positions using windowed sinc filter
  else if (weighting == WindowedSinc) {

    WindowedSincSmoothing(polydata, mask, noOfIterations, lambda);

  }
  // Perform Laplacian smoothing
  else {

    if (verbose) cout << endl;
    if (smooth_points && trackingOn) {
      FatalError("-track only implemented for -areaweighted Laplacian smoothing");
    }
    MeshSmoothing smoother;
    smoother.Input(polydata);
    smoother.Mask(mask);
    smoother.NumberOfIterations(noOfIterations);
    smoother.Lambda(lambda);
    smoother.Mu(mu);
    smoother.Sigma(-sigma1); // negative: multiple of avg. edge length
    smoother.MaximumDirectionSigma(-sigma2);
    smoother.Weighting(static_cast<MeshSmoothing::WeightFunction>(weighting));
    if (tensor_name) smoother.GeometryTensorName(tensor_name);
    if (e1_name    ) smoother.MinimumDirectionName(e1_name);
    if (e1_name    ) smoother.MaximumDirectionName(e2_name);
    smoother.SmoothPoints(smooth_points);
    smoother.SmoothArrays(scalar_names);
    smoother.AdjacentValuesOnly(adjacent_only);
    smoother.Verbose(verbose);
    smoother.Run();
    polydata = smoother.Output();

  }

  // Recompute normals if output requested and node positions changed
  if (output_normals) {
    if (smooth_points || polydata->GetPointData()->GetNormals() == NULL) {
      if (verbose) cout << "\nRecalculating normals...", cout.flush();
      vtkNew<vtkPolyDataNormals> filter;
      filter->SplittingOff();
      filter->NonManifoldTraversalOff();
      filter->ConsistencyOff();
      filter->AutoOrientNormalsOff();
      SetVTKInput(filter, polydata);
      filter->Update();
      polydata = filter->GetOutput();
      if (verbose) cout << " done" << endl;
    }
  } else {
    polydata->GetPointData()->SetNormals(NULL);
  }

  // Save output mesh
  return WritePolyData(output_name, polydata) ? 0 : 1;
}
