/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
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

#include "mirtk/Vector.h"
#include "mirtk/Matrix.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/GenericImage.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/NearestNeighborInterpolateImageFunction.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/EuclideanDistanceTransform.h"
#include "mirtk/PointCorrespondence.h"
#include "mirtk/RegisteredSurface.h"

#include "mirtk/Vtk.h"
#include "vtkGenericCell.h"
#include "vtkPolyData.h"
#include "vtkShortArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkFeatureEdges.h"
#include "vtkTriangleFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCellLocator.h"
#include "vtkOBBTree.h"
#include "vtkModifiedBSPTree.h"
#include "vtkPointLocator.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage:  " << name << " <input> <output> -image <file> [options]\n";
  cout << "        " << name << " <input> <output> -labels <file> [options]\n";
  cout << "        " << name << " <input> <output> -surface <file> [-scalars <scalars> [<scalars_new_name>] ]\n";
  cout << "        " << name << " <input> <output> -boundary [-nolabel-points] [-nolabel-cells] [-name <scalars>]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Assign scalars or labels to either the vertices or the cells of a surface mesh.\n";
  cout << "  When the input is a real-valued image, the values are linearly interpolated.\n";
  cout << "  When the input is a segmentation image, the value assigned to a vertex/cell is the\n";
  cout << "  label of the nearest voxel in the given segmentation image. For the projection of\n";
  cout << "  cortical labels onto the WM/GM or GM/CSF boundary, use :option:`-white` or :option:`-pial`.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input surface mesh.\n";
  cout << "  output   Output surface mesh.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -white\n";
  cout << "      Input surface is cortical WM/GM boundary. (default: off)\n";
  cout << "  -pial\n";
  cout << "      Input surface is cortical GM/CSF boundary. (default: off)\n";
  cout << "  -image <file>\n";
  cout << "      Input real-valued scalar/vector image.\n";
  cout << "  -labels <file>\n";
  cout << "      Input segmentation image with positive integer labels.\n";
  cout << "  -surface <file>\n";
  cout << "      Input surface from which to project scalars.\n";
  cout << "  -name <name>\n";
  cout << "      Name of output scalar array. (default: Scalars or Labels)\n";
  cout << "  -[no]celldata\n";
  cout << "      Assign values to cells of input surface. (default: off)\n";
  cout << "  -[no]pointdata\n";
  cout << "      Assign values to points of input surface. (default: on)\n";
  cout << "  -[no]fill\n";
  cout << "      Fill holes in projected surface parcellation.\n";
  cout << "  -smooth <n>\n";
  cout << "      Number of iterations to smooth. (default: 0)\n";
  cout << "  -min-size <n>\n";
  cout << "      Surface patches with less than n points are removed. (default: 0)\n";
  cout << "  -min-ratio <ratio>\n";
  cout << "      Keep only components that are larger than ratio times the size\n";
  cout << "      of the largest connected component per label, with 0 < ratio <= 1. (default: 0)\n";
  cout << "  -scalars <input_name> [<output_name>]\n";
  cout << "      Scalars to be projected from the input surface (can be defined multiple times). \n";
  cout << "      <name2> can be specified to set the name of the output scalar array\n";
  cout << "  -boundary\n";
  cout << "      Output boundary lines between surface parcels. (default: off)\n";
  cout << "      When no :option:`-image` or :option:`-labels` input file is\n";
  cout << "      specified, the boundaries of the input parcellation given by\n";
  cout << "      the labels array with the specified :option:`-name` are extracted.\n";
  cout << "  -dilation-radius <float>\n";
  cout << "      Maximum distance of voxel from input label image boundary to be assigned\n";
  cout << "      this label during dilation before projecting these labels onto the surface. (default: inf)\n";
  cout << "  -write-dilated-labels <file>\n";
  cout << "      Write image of dilated labels.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Types
// =============================================================================

typedef float                        DistanceType;
typedef float                        ScalarsType;
typedef short                        LabelType;
typedef GenericImage<DistanceType>   DistanceImage;
typedef GenericImage<ScalarsType>    ScalarsImage;
typedef GenericImage<LabelType>      LabelImage;
typedef vtkShortArray                LabelArray;
typedef OrderedMap<LabelType, long>  CountMap;
typedef CountMap::const_iterator     CountIter;
typedef OrderedSet<LabelType>        LabelSet;
typedef LabelSet::const_iterator     LabelIter;

// =============================================================================
// Project scalars (i.e., image intensities) onto surface
// =============================================================================

// -----------------------------------------------------------------------------
/// Assign linearly interpolated real-valued image values to surface points
void AssignValuesToPoints(vtkPolyData *surface, const ScalarsImage &scalars, const char *name = "Scalars")
{
  const vtkIdType noOfPoints = surface->GetNumberOfPoints();

  if (noOfPoints == 0) {
    Warning("Cannot project scalars onto surface mesh without any points!");
    return;
  }

  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(scalars.T());
  array->SetNumberOfTuples(noOfPoints);
  surface->GetPointData()->AddArray(array);

  GenericLinearInterpolateImageFunction<ScalarsImage> interp;
  interp.Input(&scalars);
  interp.Initialize();

  double p[3], *values = new double[scalars.T()];
  for (vtkIdType i = 0; i < noOfPoints; ++i) {
    surface->GetPoint(i, p);
    scalars.WorldToImage(p[0], p[1], p[2]);
    interp.Evaluate(values, p[0], p[1], p[2]);
    array->SetTuple(i, values);
  }
  delete[] values;
}

// -----------------------------------------------------------------------------
/// Assign linearly interpolated real-valued image values to surface cells
void AssignValuesToCells(vtkPolyData *surface, const ScalarsImage &scalars, const char *name = "Scalars")
{
  const vtkIdType noOfCells = surface->GetNumberOfCells();

  if (noOfCells == 0) {
    Warning("Cannot project scalars onto surface mesh without any cells!");
    return;
  }

  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(scalars.T());
  array->SetNumberOfTuples(noOfCells);
  surface->GetCellData()->AddArray(array);

  GenericLinearInterpolateImageFunction<ScalarsImage> interp;
  interp.Input(&scalars);
  interp.Initialize();

  int subId;
  vtkCell *cell;
  double p[3], pcoords[3];
  double *weights = new double[surface->GetMaxCellSize()];
  double *values  = new double[scalars.T()];

  for (vtkIdType i = 0; i < noOfCells; ++i) {
    cell  = surface->GetCell(i);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p, weights);
    scalars.WorldToImage(p[0], p[1], p[2]);
    interp.Evaluate(values, p[0], p[1], p[2]);
    array->SetTuple(i, values);
  }

  delete[] values;
  delete[] weights;
}

// =============================================================================
// Project segmentation labels onto surface
// =============================================================================

// -----------------------------------------------------------------------------
/// Returns a floating point 0/1 image to show a label.  The image needs to
/// be floating point so that it can later be used in a distance map filter.
void InitializeLabelMask(const LabelImage &labels, DistanceImage &mask, LabelType label)
{
  const int noOfVoxels = labels.NumberOfVoxels();

  DistanceType    *ptr2mask   = mask  .Data();
  const LabelType *ptr2labels = labels.Data();

  for (int i = 0; i < noOfVoxels; ++i, ++ptr2mask, ++ptr2labels) {
    *ptr2mask = static_cast<float>(*ptr2labels == label);
  }
}

// -----------------------------------------------------------------------------
LabelImage DilateLabels(const LabelImage &labels, double max_distance = inf)
{
  // Count up different labels so we can identify the number of distinct labels.
  CountMap labelCount;
  const int noOfVoxels = labels.NumberOfVoxels();
  const LabelType *ptr2label  = labels.Data();
  for (int i = 0; i < noOfVoxels; ++i, ++ptr2label) {
    if (*ptr2label > 0) ++labelCount[*ptr2label];
  }
  const int noOfLabels = static_cast<int>(labelCount.size());

  if (verbose) {
    cout << "No. of voxels        = " << noOfVoxels << endl;
    if (verbose > 1) {
      cout << "Label Counts " << endl;
      for (CountIter iter = labelCount.begin(); iter != labelCount.end(); ++iter) {
        cout << iter->first << "\t" << iter->second << endl;
      }
    }
    cout << "No. of labels        = " << noOfLabels << endl;
  }

  // Using the distance maps.
  DistanceImage minDmap(labels.Attributes());
  DistanceImage curDmap(labels.Attributes());
  DistanceImage curMask(labels.Attributes());

  // Initialise the minimum distance map.
  minDmap = numeric_limits<DistanceType>::max();

  // Note that the dilated labels are initialised to the given label image.
  // I.e. the original labels are left alone and we seek to assign labels to
  // the zero voxels based on closest labeled voxels.
  LabelImage dilatedLabels = labels;

  // Single distance transform filter for all labels.
  typedef EuclideanDistanceTransform<DistanceType> DistanceTransform;
  DistanceTransform edt(DistanceTransform::DT_3D);

  if (verbose) {
    cout << "Finding distance maps ";
    if (verbose > 1) cout << "...\nCurrent label =";
  }
  int niter = 0;
  for (CountIter iter = labelCount.begin(); iter != labelCount.end(); ++iter, ++niter) {
    const LabelType &curLabel = iter->first;
    if (verbose == 1) {
      if (niter > 0 && niter % 65 == 0) cout << "\n               ";
      cout << '.';
      cout.flush();
    } else if (verbose > 1) {
      if (niter > 0 && niter % 20 == 0) cout << "\n               ";
      cout << " " << curLabel;
      cout.flush();
    }

    // There is a new operator used for currLabelMask in the following function call.
    InitializeLabelMask(labels, curMask, curLabel);

    edt.Input (&curMask);
    edt.Output(&curDmap);
    edt.Run();

    DistanceType    *ptr2minDmap      = minDmap.Data();
    DistanceType    *ptr2dmap         = curDmap.Data();
    const LabelType *ptr2label        = labels.Data();
    LabelType       *ptr2dilatedLabel = dilatedLabels.Data();

    for (int i = 0; i < noOfVoxels; ++i, ++ptr2minDmap, ++ptr2dmap, ++ptr2label, ++ptr2dilatedLabel) {
      if (*ptr2label == 0 && *ptr2dmap < *ptr2minDmap && *ptr2dmap <= max_distance) {
        *ptr2minDmap      = *ptr2dmap;
        *ptr2dilatedLabel = curLabel;
      }
    }
  }
  if (verbose) {
    if (verbose > 1) cout << "\nFinding distance maps ...";
    cout << " done" << endl;
  }
  return dilatedLabels;
}

// -----------------------------------------------------------------------------
void LabelPoints(vtkPolyData *surface, const LabelImage &labels, const char *name = "Labels")
{
  const vtkIdType noOfPoints = surface->GetNumberOfPoints();

  if (noOfPoints == 0) {
    Warning("Cannot label points of surface mesh without any points!");
    return;
  }

  GenericNearestNeighborInterpolateImageFunction<LabelImage> nn;
  nn.Input(&labels);
  nn.Initialize();

  vtkSmartPointer<LabelArray> point_labels;
  point_labels = vtkSmartPointer<LabelArray>::New();
  point_labels->SetName(name);
  point_labels->SetNumberOfComponents(1);
  point_labels->SetNumberOfTuples(noOfPoints);

  double p[3], label;
  for (vtkIdType i = 0; i < noOfPoints; ++i) {
    surface->GetPoint(i, p);
    labels.WorldToImage(p[0], p[1], p[2]);
    label = round(nn.Evaluate(p[0], p[1], p[2]));
    point_labels->SetComponent(i, 0, label);
  }

  surface->GetPointData()->AddArray(point_labels);
}

// -----------------------------------------------------------------------------
void LabelCells(vtkPolyData *surface, const LabelImage &labels, const char *name = "Labels")
{
  const vtkIdType noOfCells = surface->GetNumberOfCells();

  if (noOfCells == 0) {
    Warning("Warning: Cannot label cells of surface mesh without any cells!");
    return;
  }

  GenericNearestNeighborInterpolateImageFunction<LabelImage> nn;
  nn.Input(&labels);
  nn.Initialize();

  vtkSmartPointer<LabelArray> cell_labels;
  cell_labels = vtkSmartPointer<LabelArray>::New();
  cell_labels->SetName(name);
  cell_labels->SetNumberOfComponents(1);
  cell_labels->SetNumberOfTuples(noOfCells);

  int subId;
  vtkCell *cell;
  double pcoords[3], p[3], label;
  Array<double> weights(surface->GetMaxCellSize());

  for (vtkIdType i = 0; i < noOfCells; ++i) {
    cell  = surface->GetCell(i);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p, weights.data());
    labels.WorldToImage(p[0], p[1], p[2]);
    label = round(nn.Evaluate(p[0], p[1], p[2]));
    cell_labels->SetComponent(i, 0, round(label));
  }

  surface->GetCellData()->AddArray(cell_labels);
}

// -----------------------------------------------------------------------------
/// Calculate surface normals
vtkSmartPointer<vtkPolyData> ComputeCellNormals(vtkPolyData *surface)
{
  vtkSmartPointer<vtkSmoothPolyDataFilter> smoother;
  smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  smoother->SetNumberOfIterations(20);
  SetVTKInput(smoother, surface);

  vtkSmartPointer<vtkPolyDataNormals> calc_normals;
  calc_normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  calc_normals->ConsistencyOff();
  calc_normals->AutoOrientNormalsOn();
  calc_normals->SplittingOff();
  calc_normals->ComputePointNormalsOff();
  calc_normals->ComputeCellNormalsOn();
  calc_normals->FlipNormalsOff();
  SetVTKConnection(calc_normals, smoother);
  calc_normals->Update();

  calc_normals->GetOutput()->GetPoints()->DeepCopy(surface->GetPoints());
  return calc_normals->GetOutput();
}

// -----------------------------------------------------------------------------
void LabelCortex(vtkPolyData *surface, const LabelImage &labels, int nsteps = 10, double h = .25, const char *name = "Labels")
{
  const vtkIdType noOfCells = surface->GetNumberOfCells();

  if (noOfCells == 0) {
    Warning("Warning: Cannot label cells of surface mesh without any cells!");
    return;
  }

  GenericNearestNeighborInterpolateImageFunction<LabelImage> nn;
  nn.Input(&labels);
  nn.Initialize();

  vtkSmartPointer<LabelArray> cell_labels;
  cell_labels = vtkSmartPointer<LabelArray>::New();
  cell_labels->SetName(name);
  cell_labels->SetNumberOfComponents(1);
  cell_labels->SetNumberOfTuples(noOfCells);

  vtkSmartPointer<vtkPolyData> mesh = ComputeCellNormals(surface);
  vtkDataArray * normals = mesh->GetCellData()->GetNormals();

  Matrix R = labels.Attributes().GetWorldToImageOrientation();
  Vector d(3);

  vtkCell *cell;
  int      subId;
  double   pcoords[3], p[3], n[3];
  double  *pweights = new double[mesh->GetMaxCellSize()];
  CountMap hist;
  LabelType label;

  for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
    // Get cell center and normal
    cell  = mesh->GetCell(cellId);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p, pweights);
    normals->GetTuple(cellId, n);
    // Convert to voxel units and scale direction vector
    labels.WorldToImage(p[0], p[1], p[2]);
    d(0) = h * n[0], d(1) = h * n[1], d(2) = h * n[2];
    d = R * d;
    // Create histogram of labels along ray in normal direction
    hist.clear();
    for (int i = 0; i < nsteps; ++i) {
      label = static_cast<LabelType>(round(nn.Evaluate(p[0], p[1], p[2])));
      if (label > 0) ++hist[label];
      p[0] += d(0), p[1] += d(1), p[2] += d(2);
    }
    // Assign label with highest frequency
    label = 0;
    long max_count = 0;
    for (CountIter i = hist.begin(); i != hist.end(); ++i) {
      if (i->second > max_count) {
        label     = i->first;
        max_count = i->second;
      }
    }
    cell_labels->SetComponent(cellId, 0, static_cast<double>(label));
  }

  surface->GetCellData()->AddArray(cell_labels);
  delete[] pweights;
}

// -----------------------------------------------------------------------------
/// Consistently label WM/cGM and cGM/CSF surfaces simultaneously
void LabelCortex(vtkPolyData *white_surface, vtkPolyData *pial_surface,
                 const LabelImage &labels, int nsamples = 10,
                 const char *name = "Labels")
{
  if (pial_surface->GetNumberOfCells() == 0 || white_surface->GetNumberOfCells() == 0) {
    Warning("Warning: Cannot label cells of surface mesh without any cells!");
    return;
  }

  GenericNearestNeighborInterpolateImageFunction<LabelImage> nn;
  nn.Input(&labels);
  nn.Initialize();

  vtkSmartPointer<LabelArray> white_labels;
  white_labels = vtkSmartPointer<LabelArray>::New();
  white_labels->SetName(name);
  white_labels->SetNumberOfComponents(1);
  white_labels->SetNumberOfTuples(white_surface->GetNumberOfCells());

  vtkSmartPointer<LabelArray> pial_labels;
  pial_labels = vtkSmartPointer<LabelArray>::New();
  pial_labels->SetName(name);
  pial_labels->SetNumberOfComponents(1);
  pial_labels->SetNumberOfTuples(pial_surface->GetNumberOfCells());

  Matrix R = labels.Attributes().GetWorldToImageOrientation();
  Vector d(3);

  const int max_cell_size = max(white_surface->GetMaxCellSize(),
                                pial_surface ->GetMaxCellSize());

  vtkCell *cell;
  int      subId;
  double   pcoords[3], p1[3], p2[3], p[3], n[3], t;
  double  *pweights = new double[max_cell_size];
  CountMap hist;
  LabelType label;
  vtkIdType whiteCellId, pialCellId;

  vtkSmartPointer<vtkPolyData> mesh = ComputeCellNormals(white_surface);
  vtkDataArray *normals = mesh->GetCellData()->GetNormals();

  vtkSmartPointer<vtkModifiedBSPTree> locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
  locator->SetDataSet(pial_surface);
  locator->BuildLocator();

  for (whiteCellId = 0; whiteCellId < mesh->GetNumberOfCells(); ++whiteCellId) {
    // Get cell center and normal
    cell  = mesh->GetCell(whiteCellId);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p1, pweights);
    normals->GetTuple(whiteCellId, n);
    // Find intersection with pial surface
    p2[0] = p1[0] + 5.0 * n[0];
    p2[1] = p1[1] + 5.0 * n[1];
    p2[2] = p1[2] + 5.0 * n[2];
    if (locator->IntersectWithLine(p1, p2, .05, t, p, pcoords, subId, pialCellId) == 0) {
      t = 2.5 / nsamples;
    } else {
      n[0] = p[0] - p1[0];
      n[1] = p[1] - p1[1];
      n[2] = p[2] - p1[2];
      t /= (nsamples - 1);
    }
    // Convert to voxel units and scale direction vector
    labels.WorldToImage(p1[0], p1[1], p1[2]);
    d(0) = t * n[0], d(1) = t * n[1], d(2) = t * n[2];
    d = R * d;
    // Create histogram of labels along ray in normal direction
    hist.clear();
    for (int i = 0; i < nsamples; ++i) {
      label = static_cast<LabelType>(round(nn.Evaluate(p1[0], p1[1], p1[2])));
      if (label > 0) ++hist[label];
      p1[0] += d(0), p1[1] += d(1), p1[2] += d(2);
    }
    // Assign label with highest frequency
    label = 0;
    long max_count = 0;
    for (CountIter i = hist.begin(); i != hist.end(); ++i) {
      if (i->second > max_count) {
        label     = i->first;
        max_count = i->second;
      }
    }
    white_labels->SetComponent(whiteCellId, 0, static_cast<double>(label));
    if (pialCellId != -1) {
      pial_labels->SetComponent(pialCellId, 0, static_cast<double>(label));
    }
  }

  mesh    = ComputeCellNormals(pial_surface);
  normals = mesh->GetCellData()->GetNormals();

  locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
  locator->SetDataSet(white_surface);
  locator->BuildLocator();

  for (vtkIdType pialCellId = 0; pialCellId < mesh->GetNumberOfCells(); ++pialCellId) {
    // Get cell center and normal
    cell  = mesh->GetCell(pialCellId);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p1, pweights);
    normals->GetTuple(pialCellId, n);
    // Find intersection with pial surface
    p2[0] = p1[0] - 10.0 * n[0];
    p2[1] = p1[1] - 10.0 * n[1];
    p2[2] = p1[2] - 10.0 * n[2];
    if (locator->IntersectWithLine(p1, p2, .05, t, p, pcoords, subId, whiteCellId) == 0) {
      t = 5.0 / nsamples;
    } else {
      n[0] = p[0] - p1[0];
      n[1] = p[1] - p1[1];
      n[2] = p[2] - p1[2];
      t /= (nsamples - 1);
    }
    // Convert to voxel units and scale direction vector
    labels.WorldToImage(p1[0], p1[1], p1[2]);
    d(0) = t * n[0], d(1) = t * n[1], d(2) = t * n[2];
    d = R * d;
    // Create histogram of labels along ray in normal direction
    hist.clear();
    label = static_cast<LabelType>(pial_labels->GetComponent(pialCellId, 0));
    if (label > 0) ++hist[label];
    for (int i = 0; i < nsamples; ++i) {
      label = static_cast<LabelType>(round(nn.Evaluate(p1[0], p1[1], p1[2])));
      if (label > 0) ++hist[label];
      p1[0] += d(0), p1[1] += d(1), p1[2] += d(2);
    }
    // Assign label with highest frequency
    label = 0;
    long max_count = 0;
    for (CountIter i = hist.begin(); i != hist.end(); ++i) {
      if (i->second > max_count) {
        label     = i->first;
        max_count = i->second;
      }
    }
    if (whiteCellId != -1) {
      if (label != static_cast<LabelType>(white_labels->GetComponent(whiteCellId, 0))) {
        label = 0;
      }
    }
    pial_labels->SetComponent(pialCellId, 0, static_cast<double>(label));
    if (whiteCellId != -1) {
      white_labels->SetComponent(whiteCellId, 0, static_cast<double>(label));
    }
  }

  white_surface->GetCellData()->AddArray(white_labels);
  pial_surface ->GetCellData()->AddArray(pial_labels);
  delete[] pweights;
}

// -----------------------------------------------------------------------------
/// Set label of unlabeled cells which are not part of the cut between
/// cortical hemispheres to -1 so they are filled in by FillHoles or SmoothLabels
void MarkHoles(vtkPolyData *surface, const char *scalars_name = "Labels", bool using_cells = true)
{
  vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
  double pcoords[3], p[3];
  int    subId;
  vtkCell *cell;

  // Get labels
  vtkSmartPointer<vtkDataArray> labels;
  if (using_cells) {
    labels = surface->GetCellData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " cell data, re-run with -celldata.");
    }
  } else {
    labels = surface->GetPointData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " point data, re-run with -pointdata.");
    }
  }

  // Keep pointer to original point scalars
  vtkSmartPointer<vtkDataArray> original_point_scalars;
  original_point_scalars = surface->GetPointData()->GetScalars();

  // Replace scalars by temporary point label mask used by connectivity filter
  vtkSmartPointer<LabelArray> point_mask;
  point_mask = vtkSmartPointer<LabelArray>::New();
  point_mask->SetNumberOfComponents(1);
  point_mask->SetNumberOfTuples(surface->GetNumberOfPoints());
  surface->GetPointData()->SetScalars(point_mask);

  // Initialize connectivity filter
  vtkSmartPointer<vtkPolyDataConnectivityFilter> filter;
  filter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  filter->SetScalarRange(0.9, 1.1);
  filter->ScalarConnectivityOn();
  filter->FullScalarConnectivityOn();
  SetVTKInput(filter, surface);

  // Label points as either unlabeled or not
  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    point_mask->SetComponent(ptId, 0, 1.);
    if (using_cells) {
      surface->GetPointCells(ptId, ids);
    } else {
      ids->SetNumberOfIds(1);
      ids->InsertId(0, ptId);
    }
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
      if (static_cast<LabelType>(labels->GetComponent(ids->GetId(i), 0)) != 0) {
        point_mask->SetComponent(ptId, 0, 0.);
        break;
      }
    }
  }

  // Find largest connected unlabeled region
  filter->SetExtractionModeToLargestRegion();
  filter->Update();

  // Replace label of not extracted unlabeled cells by -1 (these can be filled in)
  if (using_cells) {
    vtkNew<vtkCellLocator> locator;
    locator->SetDataSet(filter->GetOutput());
    locator->BuildLocator();
    Array<double> pweights(surface->GetMaxCellSize());
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      if (static_cast<LabelType>(labels->GetComponent(cellId, 0)) == 0) {
        cell  = surface->GetCell(cellId);
        subId = cell->GetParametricCenter(pcoords);
        cell->EvaluateLocation(subId, pcoords, p, pweights.data());
        if (locator->FindCell(p) == -1) {
          labels->SetComponent(cellId, 0, -1.);
        }
      }
    }
  } else {
    vtkNew<vtkPointLocator> locator;
    locator->SetDataSet(filter->GetOutput());
    locator->BuildLocator();
    for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
      if (static_cast<LabelType>(labels->GetComponent(ptId, 0)) == 0) {
        surface->GetPoint(ptId, pcoords); 
        if (locator->IsInsertedPoint(p) == -1) {
          labels->SetComponent(ptId, 0, -1.);
        }
      }
    }
  }

  // Reset original point data scalars
  surface->GetPointData()->SetScalars(original_point_scalars);
}

// -----------------------------------------------------------------------------
void MarkUnvisited(OrderedSet<vtkIdType> &cellIds, vtkDataArray *regions, vtkDataArray *labels, LabelType label)
{
  for (vtkIdType cellId = 0; cellId < regions->GetNumberOfTuples(); ++cellId) {
    if (static_cast<LabelType>(labels->GetComponent(cellId, 0)) == label) {
      regions->SetComponent(cellId, 0, -1.);
      cellIds.insert(cellId);
    } else {
      regions->SetComponent(cellId, 0, -2.);
    }
  }
}

// -----------------------------------------------------------------------------
vtkIdType NextSeed(OrderedSet<vtkIdType> &cellIds, vtkDataArray *regions)
{
  for (auto cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId) {
    if (static_cast<LabelType>(regions->GetComponent(*cellId, 0)) == -1) return *cellId;
  }
  return -1;
}

// -----------------------------------------------------------------------------
vtkIdType GrowRegion(vtkPolyData *surface, vtkDataArray *regions, LabelType regionId, vtkIdType id, bool using_cells = true)
{
  vtkSmartPointer<vtkIdList> cellPointIds    = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> idList          = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);

  EdgeTable edgeTable;
  if (!using_cells) edgeTable.Initialize(surface);

  Queue<vtkIdType> active;
  active.push(id);
  vtkIdType region_size = 0;
  while (!active.empty()) {
    id = active.front();
    active.pop();
    if (static_cast<LabelType>(regions->GetComponent(id, 0)) == -1) {
      ++region_size;
      regions->SetComponent(id, 0, regionId);
      if (using_cells) {
        surface->GetCellPoints(id, cellPointIds);
        for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); ++i) {
          idList->SetId(0, cellPointIds->GetId(i));
          if (i == cellPointIds->GetNumberOfIds() - 1) {
            idList->SetId(1, cellPointIds->GetId(0));
          } else {
            idList->SetId(1, cellPointIds->GetId(i + 1));
          }
          surface->GetCellNeighbors(id, idList, neighborCellIds);
          for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); ++j) {
            const vtkIdType neighborCellId = neighborCellIds->GetId(j);
            if (static_cast<LabelType>(regions->GetComponent(neighborCellId, 0)) == -1) {
              active.push(neighborCellId);
            }
          }
        }
      } else {
        int        adjPts;
        const int *adjIds;
        edgeTable.GetAdjacentPoints(id, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          if (static_cast<LabelType>(regions->GetComponent(adjIds[i], 0)) == -1) {
            active.push(adjIds[i]);
          }
        }
      }
    }
  }

  return region_size;
}

// -----------------------------------------------------------------------------
void MarkSmallRegions(vtkPolyData *surface, int min_region_size, const char *scalars_name = "Labels", bool using_cells = true)
{
  vtkSmartPointer<vtkDataArray> labels;
  vtkSmartPointer<LabelArray>   regions;
  regions = vtkSmartPointer<LabelArray>::New();
  regions->SetNumberOfComponents(1);

  if (using_cells) {
    labels = surface->GetCellData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " cell data, re-run with -celldata.");
    }
    regions->SetNumberOfTuples(surface->GetNumberOfCells());
  } else {
    labels = surface->GetPointData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " point data, re-run with -pointdata.");
    }
    regions->SetNumberOfTuples(surface->GetNumberOfPoints());
  }

  LabelSet label_set;
  for (vtkIdType i = 0; i < labels->GetNumberOfTuples(); ++i) {
    label_set.insert(static_cast<LabelType>(labels->GetComponent(i, 0)));
  }
  label_set.erase(0);
  label_set.erase(-1);

  vtkIdType             id;
  OrderedSet<vtkIdType> ids;
  Array<vtkIdType>      regionSz;

  surface->BuildLinks();
  for (LabelIter label = label_set.begin(); label != label_set.end(); ++label) {
    MarkUnvisited(ids, regions, labels, *label);
    while ((id = NextSeed(ids, regions)) != -1) {
      LabelType regionId = static_cast<LabelType>(regionSz.size());
      regionSz.push_back(GrowRegion(surface, regions, regionId, id, using_cells));
    }
    for (size_t i = 0; i < regionSz.size(); ++i) {
      if (regionSz[i] < static_cast<vtkIdType>(min_region_size)) {
        for (id = 0; id < surface->GetNumberOfCells(); ++id) {
          if (static_cast<size_t>(regions->GetComponent(id, 0)) == i) {
            labels->SetComponent(id, 0, -1.);
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void FillHoles(vtkPolyData *surface, const char *scalars_name = "Labels", bool using_cells = true)
{
  // Structured needed for "label front propagation"
  vtkSmartPointer<vtkIdList> idList          = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellPointIds    = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkIdType neighborCellId;
  vtkIdType id;
  LabelType label;

  idList->SetNumberOfIds(2);

  // Get cell labels
  vtkSmartPointer<vtkDataArray> labels;

  // Ensure links are build for better efficiency
  surface->BuildLinks();

  EdgeTable edgeTable;
  if (!using_cells) edgeTable.Initialize(surface);

  // Determine unlabeled boundary cells
  Queue<vtkIdType> active;
  if (using_cells) {
    labels = surface->GetCellData()->GetArray(scalars_name);
    if (labels == nullptr) {
      FatalError("Surface has no " <<  scalars_name << " cell data, re-run with -celldata.");
    }
    for (id = 0; id < surface->GetNumberOfCells(); ++id) {
      label = static_cast<LabelType>(round(labels->GetComponent(id, 0)));
      if (label != -1) continue;
      // Iterate over cell edges
      surface->GetCellPoints(id, cellPointIds);
      bool is_boundary_cell = false;
      for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); ++i) {
        // Get cells sharing these edge points
        idList->SetId(0, cellPointIds->GetId(i));
        if (i == cellPointIds->GetNumberOfIds() - 1) {
          idList->SetId(1, cellPointIds->GetId(0));
        } else {
          idList->SetId(1, cellPointIds->GetId(i + 1));
        }
        surface->GetCellNeighbors(id, idList, neighborCellIds);
        // Add if any neighboring cell is labeled
        for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); ++j) {
          neighborCellId = neighborCellIds->GetId(j);
          if (static_cast<LabelType>(labels->GetComponent(neighborCellId, 0)) != -1) {
            is_boundary_cell = true;
            break;
          }
        }
        if (is_boundary_cell) {
          active.push(id);
          break;
        }
      }
    }
  } else {
    labels = surface->GetPointData()->GetArray(scalars_name);
    if (labels == nullptr) {
      FatalError("Surface has no " <<  scalars_name << " point data, re-run with -pointdata.");
    }
    for (id = 0; id < surface->GetNumberOfPoints(); ++id) {
      label = static_cast<LabelType>(round(labels->GetComponent(id, 0)));
      if (label != -1) continue;
      int        adjPts;
      const int *adjIds;
      edgeTable.GetAdjacentPoints(id, adjPts, adjIds);
      for (int i = 0; i < adjPts; ++i) {
        if (static_cast<LabelType>(labels->GetComponent(adjIds[i], 0)) != -1) {
          active.push(id);
          break;
        }
      }
    }
  }

  // Propagate labels of neighboring cells
  CountMap hist;
  while (!active.empty()) {
    id = active.front();
    active.pop();
    // Get label of current cell
    label = static_cast<LabelType>(labels->GetComponent(id, 0));
    if (label != -1) continue;
    hist.clear();
    if (using_cells) {
      // Iterate over cell edges
      surface->GetCellPoints(id, cellPointIds);
      for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); ++i) {
        // Get cells sharing these edge points
        idList->SetId(0, cellPointIds->GetId(i));
        if (i == cellPointIds->GetNumberOfIds() - 1) {
          idList->SetId(1, cellPointIds->GetId(0));
        } else {
          idList->SetId(1, cellPointIds->GetId(i + 1));
        }
        surface->GetCellNeighbors(id, idList, neighborCellIds);
        // Count labels of neighboring cells and add unlabeled ones to queue
        for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); ++j) {
          neighborCellId = neighborCellIds->GetId(j);
          label = static_cast<LabelType>(labels->GetComponent(neighborCellId, 0));
          if (label == -1) active.push(neighborCellId);
          else ++hist[label];
        }
      }
    } else {
      int        adjPts;
      const int *adjIds;
      edgeTable.GetAdjacentPoints(id, adjPts, adjIds);
      for (int i = 0; i < adjPts; ++i) {
        label = static_cast<LabelType>(labels->GetComponent(adjIds[i], 0));
        if (label == -1) active.push(adjIds[i]);
        else ++hist[label];
      }
    }
    // Assign label with highest frequency
    label = -1;
    long max_count = 0;
    for (CountIter i = hist.begin(); i != hist.end(); ++i) {
      if (i->second > max_count) {
        label     = i->first;
        max_count = i->second;
      }
    }
    if (label == -1) active.push(id);
    else labels->SetComponent(id, 0, label);
  }
}

// -----------------------------------------------------------------------------
void SmoothLabels(vtkPolyData *surface, int niter, const char *scalars_name = "Labels", bool using_cells = true)
{
  if (niter < 1) return;

  // Structured needed for iterating over neighboring cells
  vtkSmartPointer<vtkIdList> idList          = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellPointIds    = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkIdType id, neighborCellId;
  LabelType label,  neighborLabel;

  idList->SetNumberOfIds(2);

  // Get labels
  vtkSmartPointer<vtkDataArray> labels;
  vtkSmartPointer<LabelArray> new_labels;
  new_labels = vtkSmartPointer<LabelArray>::New();
  new_labels->SetName(scalars_name);
  new_labels->SetNumberOfComponents(1);
  if (using_cells) {
    labels = surface->GetCellData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " cell data, re-run with -celldata.");
    }
    new_labels->SetNumberOfTuples(surface->GetNumberOfCells());
  } else {
    labels = surface->GetPointData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " point data, re-run with -pointdata.");
    }
    new_labels->SetNumberOfTuples(surface->GetNumberOfPoints());
  }

  // Ensure links are build for better efficiency
  surface->BuildLinks();

  EdgeTable edgeTable;
  if (!using_cells) edgeTable.Initialize(surface);

  // Replace labels by majority of neighboring cell labels
  CountMap hist;
  for (int iter = 0; iter < niter; ++iter) {
    for (id = 0; id < labels->GetNumberOfTuples(); ++id) {
      hist.clear();
      label = static_cast<LabelType>(labels->GetComponent(id, 0));
      if (label != -1) hist[label] = 1;

      if (using_cells) {
        // Iterate over cell edges
        surface->GetCellPoints(id, cellPointIds);
        for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); ++i) {
          // Get cells sharing these edge points
          idList->SetId(0, cellPointIds->GetId(i));
          if (i == cellPointIds->GetNumberOfIds() - 1) {
            idList->SetId(1, cellPointIds->GetId(0));
          } else {
            idList->SetId(1, cellPointIds->GetId(i + 1));
          }
          surface->GetCellNeighbors(id, idList, neighborCellIds);
          // Count labels of neighboring cells and add unlabeled ones to queue
          for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); ++j) {
            neighborCellId = neighborCellIds->GetId(j);
            neighborLabel  = static_cast<LabelType>(labels->GetComponent(neighborCellId, 0));
            if (neighborLabel != -1) ++hist[neighborLabel];
          }
        }
      } else {
        int        adjPts;
        const int *adjIds;
        edgeTable.GetAdjacentPoints(id, adjPts, adjIds);
        for (int i = 0; i < adjPts; ++i) {
          neighborLabel = static_cast<LabelType>(labels->GetComponent(adjIds[i], 0));
          if (neighborLabel != -1) ++hist[neighborLabel];
        }
      }
      // Assign label with highest frequency
      long max_count = 0;
      for (CountIter i = hist.begin(); i != hist.end(); ++i) {
        if (i->second > max_count) {
          label     = i->first;
          max_count = i->second;
        }
      }
      new_labels->SetComponent(id, 0, static_cast<double>(label));
    }
    labels->DeepCopy(new_labels);
  }
}


// -----------------------------------------------------------------------------
void KeepLargestRegionRatio(vtkPolyData *surface, double min_region_ratio, const char *scalars_name = "Labels", bool using_cells = true)
{
  vtkSmartPointer<vtkDataArray> labels;
  vtkSmartPointer<LabelArray>   regions;
  regions = vtkSmartPointer<LabelArray>::New();
  regions->SetNumberOfComponents(1);

  if (using_cells) {
    labels = surface->GetCellData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " cell data, re-run with -celldata.");
    }
    regions->SetNumberOfTuples(surface->GetNumberOfCells());
  } else {
    labels = surface->GetPointData()->GetArray(scalars_name);
    if (labels == NULL) {
      FatalError("Surface has no " <<  scalars_name << " point data, re-run with -pointdata.");
    }
    regions->SetNumberOfTuples(surface->GetNumberOfPoints());
  }

  LabelSet label_set;
  for (vtkIdType i = 0; i < labels->GetNumberOfTuples(); ++i) {
    label_set.insert(static_cast<LabelType>(labels->GetComponent(i, 0)));
  }
  label_set.erase(0);
  label_set.erase(-1);

  vtkIdType             id;
  OrderedSet<vtkIdType> ids;

  surface->BuildLinks();

  for (LabelIter label = label_set.begin(); label != label_set.end(); ++label) {
    Array<vtkIdType> regionSz;
    MarkUnvisited(ids, regions, labels, *label);
    while ((id = NextSeed(ids, regions)) != -1) {
      LabelType regionId = static_cast<LabelType>(regionSz.size());
      regionSz.push_back(GrowRegion(surface, regions, regionId, id, using_cells));
    }
    vtkIdType min_region_size = -1;
    for (size_t i = 0; i < regionSz.size(); ++i) {
      if (min_region_size < regionSz[i]) {
        min_region_size = regionSz[i];
      }
    }
    min_region_size = iround(min_region_size * min_region_ratio);
    for (size_t i = 0; i < regionSz.size(); ++i) {
      if (regionSz[i] < min_region_size) {
        for (id = 0; id < regions->GetNumberOfTuples(); ++id) {
          if (static_cast<size_t>(regions->GetComponent(id, 0)) == i) {
            labels->SetComponent(id, 0, -1.);
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ExtractLabelBoundaries(vtkPolyData *surface, const char *scalars_name = "Labels")
{
  vtkSmartPointer<vtkPolyData>  mesh;
  vtkSmartPointer<vtkCellArray> edges;
  edges = vtkSmartPointer<vtkCellArray>::New();

  // Clean input surface
  if (verbose) cout << "Cleaning surface mesh ...", cout.flush();
  vtkSmartPointer<vtkCleanPolyData> cleaner;
  cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  SetVTKInput(cleaner, surface);
  cleaner->PointMergingOff();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetTolerance(.0);
  cleaner->Update();
  mesh = cleaner->GetOutput();
  if (verbose) cout << " done" << endl;

  // Triangulate input surface
  if (verbose) cout << "Triangulating surface mesh ...", cout.flush();
  vtkSmartPointer<vtkTriangleFilter> triangulate;
  triangulate = vtkSmartPointer<vtkTriangleFilter>::New();
  SetVTKInput(triangulate, mesh);
  triangulate->PassVertsOff();
  triangulate->PassLinesOff();
  triangulate->Update();
  mesh = triangulate->GetOutput();
  if (verbose) cout << " done" << endl;

  // Ensure links are build for better efficiency
  mesh->BuildLinks();

  // Get cell labels
  vtkSmartPointer<vtkDataArray> cell_labels;
  cell_labels = mesh->GetCellData()->GetArray(scalars_name);
  if (cell_labels == NULL) {
    FatalError("Surface has no " <<  scalars_name << " cell data, re-run with -celldata.");
  }

  // Iterate over cells
  if (verbose) cout << "Determining boundary edges ...", cout.flush();

  LabelType label;
  bool is_boundary_edge;
  vtkSmartPointer<vtkIdList> idList          = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellPointIds    = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkIdType neighborCellId;

  idList->SetNumberOfIds(2);

  OrderedSet<Pair<vtkIdType, vtkIdType> > boundaryEdgeIds;
  for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
    // Get cell label and point IDs
    label = static_cast<LabelType>(round(cell_labels->GetComponent(cellId, 0)));
    cellPointIds->Reset();
    mesh->GetCellPoints(cellId, cellPointIds);
    // Iterate over cell edges
    for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); ++i) {
      is_boundary_edge = false;
      // Get cells sharing these edge points
      idList->SetId(0, cellPointIds->GetId(i));
      if (i == cellPointIds->GetNumberOfIds() - 1) {
        idList->SetId(1, cellPointIds->GetId(0));
      } else {
        idList->SetId(1, cellPointIds->GetId(i + 1));
      }
      neighborCellIds->Reset();
      mesh->GetCellNeighbors(cellId, idList, neighborCellIds);
      // Compare labels of neighboring cells to this cell's label
      for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); ++j) {
        neighborCellId = neighborCellIds->GetId(j);
        if (label != static_cast<LabelType>(round(cell_labels->GetComponent(neighborCellId, 0)))) {
          is_boundary_edge = true;
        }
      }
      // Add edge if it separates differently labeled cells
      if (is_boundary_edge) {
        boundaryEdgeIds.insert(MakePair(idList->GetId(0), idList->GetId(1)));
      }
    }
  }
  OrderedSet<Pair<vtkIdType, vtkIdType> >::iterator i = boundaryEdgeIds.begin();
  while (i != boundaryEdgeIds.end()) {
    idList->SetId(0, i->first);
    idList->SetId(1, i->second);
    edges->InsertNextCell(idList);
    ++i;
  }
  if (verbose) cout << " done: #edges = " << edges->GetNumberOfCells() << endl;

  mesh->SetVerts(NULL);
  mesh->SetLines(edges);
  mesh->SetPolys(NULL);
  mesh->SetStrips(NULL);
  mesh->GetCellData()->Initialize();

  cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  SetVTKInput(cleaner, mesh);
  cleaner->PointMergingOff();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetTolerance(.0);
  cleaner->Update();
  return cleaner->GetOutput();
}


// -----------------------------------------------------------------------------
void CheckBounds(vtkPolyData *surface, const LabelImage &labels, bool verbose, string surface_name="")
{
  double bounds[6];
  for (int b=0; b<6; b+=2) bounds[b] = numeric_limits<double>::max();
  for (int b=1; b<6; b+=2) bounds[b] = numeric_limits<double>::min();

  Point p;
  const vtkIdType npoints = surface->GetNumberOfPoints();
  for (vtkIdType i = 0; i < npoints; ++i) {
    surface->GetPoint(i, p);
    labels.WorldToImage(p);
    for (int d = 0; d < 3; ++d) {
      if (p(d) < bounds[d*2  ]) bounds[d*2  ] = p(d);
      if (p(d) > bounds[d*2+1]) bounds[d*2+1] = p(d);
    }
  }

  const double &xmin = bounds[0];
  const double &xmax = bounds[1];
  const double &ymin = bounds[2];
  const double &ymax = bounds[3];
  const double &zmin = bounds[4];
  const double &zmax = bounds[5];

  if (verbose) { 
    cout << "Bounds of " << surface_name << "surface in image coordinates = ";
    cout << "(" << xmin << ", " << ymin << ", " << zmin << ") and ";
    cout << "(" << xmax << ", " << ymax << ", " << zmax << ")" << endl;
  }

  if (xmin < -0.5 || xmax > labels.X() - .5 ||
      ymin < -0.5 || ymax > labels.Y() - .5 ||
      zmin < -0.5 || zmax > labels.Z() - .5) {
    FatalError("Surface outside bounds of input image!");
  }
} 

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Positional arguments
  REQUIRES_POSARGS(2);

  const char *input_image_name   = NULL;
  const char *input_labels_name  = NULL;
  const char *input_surface_name = NULL;
  const char *output_name        = NULL;
  const char *input_white_name   = NULL;
  const char *input_pial_name    = NULL;
  const char *output_white_name  = NULL;
  const char *output_pial_name   = NULL;

  bool csf_gm_boundary = false;
  bool gm_wm_boundary  = false;

  if (NUM_POSARGS == 2) {
    input_surface_name = POSARG(1);
    output_name        = POSARG(2);
  } else if (NUM_POSARGS == 4) {
    input_white_name   = POSARG(1);
    input_pial_name    = POSARG(2);
    output_white_name  = POSARG(3);
    output_pial_name   = POSARG(4);
    csf_gm_boundary = gm_wm_boundary = true;
  } else {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Optional arguments
  const char *output_label_image_name = NULL;
  const char *output_scalars_name     = NULL;
  bool        output_boundary_edges   = false;
  bool        label_cells             = false;
  bool        label_points            = false;
  int         min_region_size         = 0;
  double      min_region_ratio        = 0;
  bool        fill_holes              = true;
  int         smoothing_iterations    = 0;
  double      max_dilation_distance   = inf;

  // arguments for scalars projection from surface
  const char *input_surface_proj_name   = NULL;
  Array<const char *> surface_proj_input_scalars_names;
  Array<const char *> surface_proj_output_scalars_names;

  // Parse remaining arguments
  for (ALL_OPTIONS) {
    if      (OPTION("-write-dilated-labels")) output_label_image_name = ARGUMENT;
    else if (OPTION("-dilation-radius")) PARSE_ARGUMENT(max_dilation_distance);
    else if (OPTION("-name"))   output_scalars_name = ARGUMENT;
    else if (OPTION("-image"))  input_image_name    = ARGUMENT;
    else if (OPTION("-labels")) input_labels_name   = ARGUMENT;
    else if (OPTION("-min-size")) PARSE_ARGUMENT(min_region_size);
    else if (OPTION("-min-ratio")) PARSE_ARGUMENT(min_region_ratio);
    else if (OPTION("-smooth")) PARSE_ARGUMENT(smoothing_iterations);
    else if (OPTION("-pial")) {
      csf_gm_boundary = true;
      gm_wm_boundary  = false;
    }
    else if (OPTION("-white")) {
      csf_gm_boundary = false;
      gm_wm_boundary  = true;
    }
    else HANDLE_BOOLEAN_OPTION("fill", fill_holes);
    else HANDLE_BOOLEAN_OPTION("celldata", label_cells);
    else HANDLE_BOOLEAN_OPTION("pointdata", label_points);
    else HANDLE_BOOLEAN_OPTION("boundary", output_boundary_edges);
    else if (OPTION("-surface"))  input_surface_proj_name = ARGUMENT;
    else if (OPTION("-scalars")){
      char* scalars_to_copy = ARGUMENT;
      surface_proj_input_scalars_names.push_back(scalars_to_copy);
      if (HAS_ARGUMENT) scalars_to_copy = ARGUMENT;
      surface_proj_output_scalars_names.push_back(scalars_to_copy);
    } 
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (input_image_name) {
    if (!output_scalars_name) output_scalars_name = "Scalars";
  } else if (input_labels_name) {
    if (!output_scalars_name) output_scalars_name = "Labels";
    if (!label_cells && !label_points) {
      // Labeling the cells makes more sense, but for backward compatibility reasons
      // the default is to label the points...
      if (csf_gm_boundary || gm_wm_boundary) {
        label_cells = true;
      } else {
        label_points = true;
      }
    }
  } else if (!input_surface_proj_name && !output_boundary_edges) {
    FatalError("Input -image or -labels or -surface or -boundary required!");
  }

  // Read input surface
  vtkSmartPointer<vtkPolyData> surface, white_surface, pial_surface;
  if (input_surface_name) {
    if (verbose) cout << "Reading surface...", cout.flush();
    surface = ReadPolyData(input_surface_name);
    if (verbose) cout << " done" << endl;
  }
  if (input_white_name) {
    if (verbose) cout << "Reading  WM/cGM surface...", cout.flush();
    white_surface = ReadPolyData(input_white_name);
    if (verbose) cout << " done" << endl;
  }
  if (input_pial_name) {
    if (verbose) cout << "Reading cGM/CSF surface...", cout.flush();
    pial_surface = ReadPolyData(input_pial_name);
    if (verbose) cout << " done" << endl;
  }

  // Read values from input image
  ScalarsImage image;
  LabelImage   labels;
  if (input_image_name) {
    if (verbose) cout << "Reading image...", cout.flush();
    image.Read(input_image_name);
    if (verbose) cout << " done" << endl;
  }
  if (input_labels_name) {
    if (verbose) cout << "Reading labels...", cout.flush();
    labels.Read(input_labels_name);
    if (verbose) cout << " done" << endl;
  }

  // Check that surface does not go outside fov of label image.
  if (surface && (input_image_name || input_labels_name) && (label_points || label_cells)) {
    CheckBounds(surface, labels, verbose > 0, "");
  }
  if (white_surface && (input_image_name || input_labels_name) && (label_points || label_cells)) {
    CheckBounds(white_surface, labels, verbose > 0, "WM/cGM");
  }
  if (pial_surface && (input_image_name || input_labels_name) && (label_points || label_cells)) {
    CheckBounds(white_surface, labels, verbose > 0, "cGM/CSF");
  }

  // Compute dilated labels
  LabelImage dilatedLabels;
  if (input_labels_name && (label_points || (label_cells && !csf_gm_boundary && !gm_wm_boundary))) {
    dilatedLabels = DilateLabels(labels, max_dilation_distance);
    if (output_label_image_name) dilatedLabels.Write(output_label_image_name);
  }

  // Assign labels to points (vertices)
  if (label_points) {
    if (input_labels_name) {
      if (surface) {
        if (verbose) cout << "Labeling the input vertices...", cout.flush();
        LabelPoints(surface, dilatedLabels, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
      if (white_surface) {
        if (verbose) cout << "Labeling vertices of  WM/cGM surface...", cout.flush();
        LabelPoints(white_surface, dilatedLabels, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
      if (pial_surface) {
        if (verbose) cout << "Labeling vertices of cGM/CSF surface...", cout.flush();
        LabelPoints(pial_surface, dilatedLabels, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
    } else if (input_image_name) {
      if (surface) {
        if (verbose) cout << "Assigning values to the vertices of the input surface...", cout.flush();
        AssignValuesToPoints(surface, image, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
      if (white_surface) {
        if (verbose) cout << "Assigning values to the vertices of the WM/cGM surface...", cout.flush();
        AssignValuesToPoints(white_surface, image, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
      if (pial_surface) {
        if (verbose) cout << "Assigning values to the vertices of the cGM/CSF surface...", cout.flush();
        AssignValuesToPoints(pial_surface, image, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
    }
  }

  // Assign labels to cells (faces, triangle, ...)
  if (label_cells) {
    if (input_labels_name) {
      if (csf_gm_boundary && gm_wm_boundary) {
        if (verbose) cout << "Labeling cortical surfaces...", cout.flush();
        LabelCortex(white_surface, pial_surface, labels, 10, output_scalars_name);
      } else if (gm_wm_boundary) {
        if (verbose) cout << "Labeling WM/cGM surface...", cout.flush();
        LabelCortex(surface, labels,  5, +.2, output_scalars_name);
      } else if (csf_gm_boundary) {
        if (verbose) cout << "Labeling cGM/CSF surface...", cout.flush();
        LabelCortex(surface, labels, 10, -.2, output_scalars_name);
      } else {
        if (verbose) cout << "Labeling the input cells...", cout.flush();
        LabelCells(surface, dilatedLabels, output_scalars_name);
      }
      if (verbose) cout << " done" << endl;
    } else if (input_image_name) {
      if (surface) {
        if (verbose) cout << "Assigning values to cells of input surface...", cout.flush();
        AssignValuesToCells(surface, image, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
      if (white_surface) {
        if (verbose) cout << "Assigning values to cells of WM/cGM surface...", cout.flush();
        AssignValuesToCells(white_surface, image, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
      if (pial_surface) {
        if (verbose) cout << "Assigning values to cells of cGM/CSF surface...", cout.flush();
        AssignValuesToCells(pial_surface, image, output_scalars_name);
        if (verbose) cout << " done" << endl;
      }
    }
  }


  if(input_labels_name && (label_cells || label_points)){
    bool using_cells;
    for (int i=0; i<2; i++){
      if(i==0 && !label_points) continue;
      if(i==1 && !label_cells) continue;
      using_cells = i;

      if (fill_holes || smoothing_iterations > 0) {
        if (verbose) cout << "Marking holes in parcellation...", cout.flush();
        if (surface)       MarkHoles(surface,       output_scalars_name, using_cells);
        if (white_surface) MarkHoles(white_surface, output_scalars_name, using_cells);
        if (pial_surface)  MarkHoles(pial_surface,  output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }
      if (fill_holes) {
        if (verbose) cout << "Filling holes...", cout.flush();
        if (surface)       FillHoles(surface,       output_scalars_name, using_cells);
        if (white_surface) FillHoles(white_surface, output_scalars_name, using_cells);
        if (pial_surface)  FillHoles(pial_surface,  output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }

      if (min_region_size > 1) {
        if (verbose) cout << "Marking parcels with less than " << min_region_size << " cells...", cout.flush();
        if (surface)       MarkSmallRegions(surface,       min_region_size, output_scalars_name, using_cells);
        if (white_surface) MarkSmallRegions(white_surface, min_region_size, output_scalars_name, using_cells);
        if (pial_surface)  MarkSmallRegions(pial_surface,  min_region_size, output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }else if (min_region_ratio > 0){
        if (verbose) cout << "Keep largest parcel per label ...", cout.flush();
        if (surface)       KeepLargestRegionRatio(surface,       min_region_ratio, output_scalars_name, using_cells);
        if (white_surface) KeepLargestRegionRatio(white_surface, min_region_ratio, output_scalars_name, using_cells);
        if (pial_surface)  KeepLargestRegionRatio(pial_surface,  min_region_ratio, output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }

      if (smoothing_iterations > 0) {
        if (verbose) cout << "Smoothing/filling parcels (niter=" << smoothing_iterations << ")...", cout.flush();
        if (surface)       SmoothLabels(surface, smoothing_iterations,       output_scalars_name, using_cells);
        if (white_surface) SmoothLabels(white_surface, smoothing_iterations, output_scalars_name, using_cells);
        if (pial_surface)  SmoothLabels(pial_surface,  smoothing_iterations, output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }

      if (min_region_size > 1) {
        if (verbose) cout << "Marking parcels with less than " << min_region_size << " cells...", cout.flush();
        if (surface)       MarkSmallRegions(surface,       min_region_size, output_scalars_name, using_cells);
        if (white_surface) MarkSmallRegions(white_surface, min_region_size, output_scalars_name, using_cells);
        if (pial_surface)  MarkSmallRegions(pial_surface,  min_region_size, output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }else if (min_region_ratio > 0){
        if (verbose) cout << "Keep largest parcel per label ...", cout.flush();
        if (surface)       KeepLargestRegionRatio(surface,       min_region_ratio, output_scalars_name, using_cells);
        if (white_surface) KeepLargestRegionRatio(white_surface, min_region_ratio, output_scalars_name, using_cells);
        if (pial_surface)  KeepLargestRegionRatio(pial_surface,  min_region_ratio, output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }

      if (fill_holes) {
        if (verbose) cout << "Filling remaining holes...", cout.flush();
        if (surface)       FillHoles(surface,       output_scalars_name, using_cells);
        if (white_surface) FillHoles(white_surface, output_scalars_name, using_cells);
        if (pial_surface)  FillHoles(pial_surface,  output_scalars_name, using_cells);
        if (verbose) cout << " done" << endl;
      }
    }
  }

  // scalars projection from surface
  if (surface && input_surface_proj_name) {
    vtkSmartPointer<vtkPolyData> surface_proj = ReadPolyData(input_surface_proj_name);

    if (surface_proj_input_scalars_names.empty()) {
      vtkDataArray *scalars_array = surface_proj->GetPointData()->GetScalars();
      if (scalars_array == nullptr) {
        FatalError("Scalars need to be specified with the -scalars option");
      }
      surface_proj_input_scalars_names .push_back(scalars_array->GetName());
      surface_proj_output_scalars_names.push_back(scalars_array->GetName());
    }

    RegisteredSurface target, source;
    target.InputSurface(surface);
    source.InputSurface(surface_proj);
    target.Initialize();
    source.Initialize();
    target.Update();
    source.Update();

    UniquePtr<PointCorrespondence> cmap;
    cmap.reset(PointCorrespondence::New(PointCorrespondence::ClosestPoint));
    cmap->FromTargetToSource(true);
    cmap->FromSourceToTarget(false);
    cmap->Target(&target);
    cmap->Source(&source);
    cmap->Initialize();
    cmap->Update();

    const vtkIdType noOfPoints = surface->GetNumberOfPoints();
    for (size_t a = 0; a < surface_proj_input_scalars_names.size(); ++a) {
      vtkDataArray *source_array = surface_proj->GetPointData()->GetArray(surface_proj_input_scalars_names[a]);
      if (source_array == nullptr) {
        FatalError("Surface has no scalars: " <<  surface_proj_input_scalars_names[a]);
      }
      vtkSmartPointer<vtkDataArray> array;
      array.TakeReference(source_array->NewInstance());
      array->SetName(surface_proj_output_scalars_names[a]);
      array->SetNumberOfComponents(source_array->GetNumberOfComponents());
      array->SetNumberOfTuples(noOfPoints);
      surface->GetPointData()->AddArray(array);
      for (vtkIdType i = 0; i < noOfPoints; ++i) {
        array->SetTuple(i, source_array->GetTuple(cmap->GetIndex(i)));
      }
    }
  }

  // Write dataset with line segments at boundary between differently labeled regions
  const char *what = "output";
  if (input_labels_name) what = "labeled";
  if (output_boundary_edges) {
    what = "boundaries of labeled regions of";
    if (surface) {
      surface = ExtractLabelBoundaries(surface, output_scalars_name);
    }
    if (white_surface) {
      white_surface = ExtractLabelBoundaries(white_surface, output_scalars_name);
    }
    if (pial_surface) {
      pial_surface = ExtractLabelBoundaries(pial_surface, output_scalars_name);
    }
  }
  if (white_surface && output_white_name) {
    if (verbose) cout << "Writing " << what << "  WM/cGM surface ...", cout.flush();
    WritePolyData(output_white_name, white_surface);
    if (verbose) cout << " done" << endl;
  }
  if (pial_surface && output_pial_name) {
    if (verbose) cout << "Writing " << what << " cGM/CSF surface ...", cout.flush();
    WritePolyData(output_pial_name, pial_surface);
    if (verbose) cout << " done" << endl;
  }
  if (surface && output_name) {
    if (verbose) cout << "Writing " << what << " surface ...", cout.flush();
    WritePolyData(output_name, surface);
    if (verbose) cout << " done" << endl;
  }

  return 0;
}
