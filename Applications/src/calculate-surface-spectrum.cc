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

#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/RegisteredSurface.h"
#include "mirtk/PointCorrespondence.h"
#include "mirtk/FuzzyCorrespondence.h"
#include "mirtk/SpectralDecomposition.h"

#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"

using namespace mirtk;
using namespace mirtk::SpectralDecomposition;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << "       " << name << " <input1> <input2> <output1> <output2> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Performs a spectral analysis of the general graph Laplacian matrix" << endl;
  cout << "  computed from the given input surface mesh(es)." << endl;
  cout << endl;
  cout << "  If :option:`-target` surface is specified, the eigenvectors are reordered" << endl;
  cout << "  and their sign flipped to match the eigenvectors of the target surface mesh." << endl;
  cout << endl;
  cout << "  If two input and output surface meshes are specified, a spectral analysis" << endl;
  cout << "  of the joint graph Laplacian is performed after an initial surface match." << endl;
  cout << "  The type of initial point correspondences is defined by :option:`-corr`." << endl;
  cout << endl;
  cout << "  The implementation is based on the MATLAB code by Herve Lombaert of [1]_." << endl;
  cout << endl;
  cout << "  .. [1] H. Lombaert, J. Sporring, and K. Siddiqi." << endl;
  cout << "         Diffeomorphic Spectral Matching of Cortical Surfaces." << endl;
  cout << "         In the 23rd Image Processing in Medical Imaging (IPMI), 2013." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    File name of input  surface (vtkPolyData)." << endl;
  cout << "  output   File name of output surface (vtkPolyData)." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -target <file>            File name of dataset (vtkPolyData) whose spectral components" << endl;
  cout << "                            should be used to adjust the sign and order of the eigenmodes." << endl;
  cout << "  -corr <string>            Type of point correspondences to use for initial correspondences:" << endl;
  cout << "                            \"closest point\" (\"CP\"), \"spectral match\" (\"SM\"). (default: SM)" << endl;
  cout << "  -corrpar <name>=<value>   Name/value pair of correspondence parameter." << endl;
  cout << "  -corrin <file>            Text file naming for each point in input1 the index of the corresponding" << endl;
  cout << "                            point in the second dataset input2, one index per line. (default: :option:`-corr`)" << endl;
  cout << "  -points                   Replace output vertex coordinates by first three spectral coordinates." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Write point set with surface eigenmodes
void Write(const char *fname, vtkPointSet *pointset, vtkPolyData *surface, FileOption fopt, bool as_points = false)
{
  vtkSmartPointer<vtkPointSet> output = pointset;

  vtkDataArray *modes = surface->GetPointData()->GetArray("joint_eigenmodes");
  if (!modes)   modes = surface->GetPointData()->GetArray("eigenmodes");
  if (!modes) {
    FatalError("Output surface mesh has no eigenmodes point data!");
  }

  if (as_points) {

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(surface->GetNumberOfPoints());

    double *p = new double[max(3, int(modes->GetNumberOfComponents()))];
    for (int c = 0; c < 3; ++c) p[c] = .0;
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
      modes ->GetTuple(i, p);
      points->SetPoint(i, p);
    }
    delete[] p;

    output = vtkSmartPointer<vtkPolyData>::New();
    output->ShallowCopy(surface);
    output->SetPoints(points);

  } else if (pointset != surface) {

    vtkSmartPointer<vtkDataArray> output_modes;
    output_modes.TakeReference(modes->NewInstance());
    output_modes->SetName(modes->GetName());
    output_modes->SetNumberOfComponents(modes->GetNumberOfComponents());
    output_modes->SetNumberOfTuples(pointset->GetNumberOfPoints());

    for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
      for (int j = 0; j < modes->GetNumberOfComponents(); ++j) {
        output_modes->SetComponent(ptId, j, .0);
      }
    }

    vtkIdType origPtId;
    vtkDataArray *origPtIds = surface->GetPointData()->GetArray("vtkOriginalPointIds");
    for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
      origPtId = static_cast<vtkIdType>(origPtIds->GetComponent(ptId, 0));
      for (int j = 0; j < modes->GetNumberOfComponents(); ++j) {
        output_modes->SetComponent(origPtId, j, modes->GetComponent(ptId, j));
      }
    }

    output->GetPointData()->AddArray(output_modes);
  }

  WritePointSet(fname, output, fopt);
}

// -----------------------------------------------------------------------------
/// Compute eigenmodes of joint intra- and inter-mesh graph Laplacian matrix
int ComputeJointEigenmodes(vtkPolyData *target, const Array<int> &target_sample,
                           vtkPolyData *source, const Array<int> &source_sample,
                           const PointCorrespondence *corr, int k,
                           Matrix &modes, Vector &freq)
{
  using SpectralDecomposition::SparseMatrix;
  // Set intra-mesh affinity weights
  const int m = target->GetNumberOfPoints();
  const int n = source->GetNumberOfPoints();
  SparseMatrix::Entries *cols = new SparseMatrix::Entries[m + n];
  AdjacencyMatrix(cols, SparseMatrix::CCS, 0, 0, target);
  AdjacencyMatrix(cols, SparseMatrix::CCS, m, m, source);
  // Set inter-mesh affinity weights
  const FuzzyCorrespondence *fuzzy;
  if ((fuzzy = dynamic_cast<const FuzzyCorrespondence *>(corr))) {
    ConnectivityMatrix(cols, SparseMatrix::CCS, 0, m, target, &target_sample,
                       source, &source_sample, &fuzzy->Weight(), false);
    ConnectivityMatrix(cols, SparseMatrix::CCS, m, 0, source, &source_sample,
                       target, &target_sample, &fuzzy->Weight(), true );
  } else {
    int i, j;
    double w, p1[3], p2[3];
    for (size_t s = 0; s < target_sample.size(); ++s) {
      i = target_sample[s];
      j = corr->GetSourceIndex(i);
      if (j < 0) {
        FatalError("Point correspondence must be fuzzy or map point indices!");
      }
      target->GetPoint(i, p1);
      source->GetPoint(j, p2);
      w = 1.0 / (sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + EPSILON);
      cols[i  ].push_back(MakePair(j+m, w));
      cols[j+m].push_back(MakePair(i,   w));
    }
    for (size_t s = 0; s < source_sample.size(); ++s) {
      j = source_sample[s];
      i = corr->GetTargetIndex(j);
      if (i < 0) {
        FatalError("Point correspondence must be fuzzy or map point indices!");
      }
      target->GetPoint(i, p1);
      source->GetPoint(j, p2);
      w = 1.0 / (sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + EPSILON);
      cols[j+m].push_back(MakePair(i,   w));
      cols[i  ].push_back(MakePair(j+m, w));
    }
  }
  // Compute graph Laplacian of joint connectivity graph
  SparseMatrix L(SparseMatrix::CCS);
  L.Initialize(m + n, m + n, cols);
  delete[] cols;
  NormalizedLaplacian(L, L);
  // Compute eigenmodes of joint graph Laplacian
  return ComputeEigenmodes(L, k+1, modes, freq);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  const char *input_name [2] = {nullptr};
  const char *output_name[2] = {nullptr};
  int ndatasets = 0;

  switch (NUM_POSARGS) {
    case 2:
      input_name [0] = POSARG(1);
      output_name[0] = POSARG(2);
      ndatasets = 1;
      break;
    case 4:
      input_name [0] = POSARG(1);
      input_name [1] = POSARG(2);
      output_name[0] = POSARG(3);
      output_name[1] = POSARG(4);
      ndatasets = 2;
      break;
    default:
      PrintHelp(EXECNAME);
      FatalError("Invalid number of positional arguments (" << NUM_POSARGS << ")");
  };

  // Optional arguments
  const char *target_name = nullptr;
  const char *dofin_name  = nullptr;
  int         k           = 5;
  FileOption  fopt        = FO_Default;
  bool        as_points   = false;

  PointCorrespondence::TypeId ctype = PointCorrespondence::ClosestPoint;
  ParameterList               cparam;
  Array<string>               cfeature_name;
  Array<double>               cfeature_weight;

  for (ALL_OPTIONS) {
    if      (OPTION("-target")) target_name  = ARGUMENT;
    else if (OPTION("-dofin"))  dofin_name   = ARGUMENT;
    else if (OPTION("-k") || OPTION("-dim")) k = atoi(ARGUMENT);
    else if (OPTION("-corr")) PARSE_ARGUMENT(ctype);
    else if (OPTION("-corrin")) {
      ctype = PointCorrespondence::FiducialMatch;
      Insert(cparam, "Correspondence map", ARGUMENT);
    }
    else if (OPTION("-corrpar")) {
      Insert(cparam, ARGUMENT, ARGUMENT);
    }
    else if (OPTION("-feature")) {
      cfeature_name.push_back(ARGUMENT);
      if (HAS_ARGUMENT) {
        double w;
        PARSE_ARGUMENT(w);
        cfeature_weight.push_back(w);
      }
      else cfeature_weight.push_back(1.0);
    }
    else if (OPTION("-points")) as_points = true;
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Default correspondence features and normalization of names to lowercase
  if (cfeature_name.empty()) {
    cfeature_name  .push_back("eigenmodes");
    cfeature_weight.push_back(1.0);
  } else {
    for (size_t i = 0; i < cfeature_name.size(); ++i) {
      transform(cfeature_name[i].begin(), cfeature_name[i].end(),
                cfeature_name[i].begin(), ::tolower);
    }
  }

  // Ensure number of spectral components is positive
  if (k <= 0) {
    FatalError("Argument for option -k/-dim must be positive!");
  }
  if (as_points && k < 3) {
    FatalError("Argument for option -k/-dim must be >= 3 when -points should be written!");
  }

  // Read input datasets
  vtkSmartPointer<vtkPointSet> dataset[2];
  vtkSmartPointer<vtkPolyData> surface[2];
  for (int i = 0; i < ndatasets; ++i) {
    if (verbose > 1) cout << "Reading point set " << (i+1) << " from " << input_name[i] << endl;
    FileOption opt;
    dataset[i] = ReadPointSet(input_name[i], opt);
    if (fopt == FO_Default) fopt = opt;
    if (dataset[i] == NULL || dataset[i]->GetNumberOfPoints() == 0) {
      FatalError("Failed to read dataset " << (i+1) << " or dataset is empty!");
    }
    surface[i] = DataSetSurface(dataset[i], true);
  }

  // Transform first surface by specified transformation
  if (dofin_name) {
    UniquePtr<Transformation> dof(Transformation::New(dofin_name));
    RegisteredSurface transformed;
    transformed.InputSurface(surface[0]);
    transformed.Transformation(dof.get());
    transformed.Initialize();
    transformed.Update();
    surface[0] = transformed.Surface();
  }

  // Skip individual spectral analysis if joint analysis requested based on
  // pre-defined correspondences, non-spectral matches, or known eigenmodes
  bool individual_analysis[2] = {false};
  bool spectral_match         =  false;
  if (ndatasets == 1) {
    spectral_match = !as_points; // -points option can be applied to single input dataset
  } else if (ndatasets > 1) {
    switch (ctype) {
      case PointCorrespondence::FiducialMatch:
        spectral_match = false;
        break;
      case PointCorrespondence::SpectralMatch:
        spectral_match = true;
        break;
      default:
        spectral_match = false;
        for (size_t i = 0; i < cfeature_name.size(); ++i) {
          if (cfeature_name[i].find("spectral" ) != string::npos ||
              cfeature_name[i].find("eigenmode") != string::npos) {
            spectral_match = true;
            break;
          }
        }
    }
  }
  if (spectral_match) {
    for (int i = 0; i < ndatasets; ++i) {
      individual_analysis[i] = (!surface[i]->GetPointData()->HasArray("eigenmodes") &&
                                !surface[i]->GetPointData()->HasArray("joint_eigenmodes"));
    }
  }

  // Individual spectral analysis of input datasets
  Matrix modes[2];
  Vector freq [2];
  for (int i = 0; i < ndatasets; ++i) {
    if (individual_analysis[i]) {
      if (ComputeEigenmodes(surface[i], k, modes[i], freq[i]) < k) {
        FatalError("Failed to find " << k << " eigenmodes of dataset " << (i+1));
      }
      if (verbose) {
        cout << "Frequencies";
        if (ndatasets > 1) cout << " of dataset " << (i+1);
        cout << ": ";
        for (int c = 0; c < k; ++c) {
          if (c > 0) cout << "  ";
          cout << scientific << freq[i](c);
        }
        cout << endl;
      }
    }
  }

  // Joint spectral analysis
  if (ndatasets == 2) {

    // Adjust sign and order of individually computed eigenmodes
    if (spectral_match) {
      Vector cost = MatchEigenmodes(surface[0]->GetPoints(), modes[0], freq[0],
                                        surface[1]->GetPoints(), modes[1], freq[1]);

      // TODO: Transform modes[1] using coherent point drift algorithm

      // Add individual eigenmodes to output point data
      for (int i = 0; i < ndatasets; ++i) {
        SetEigenmodes(surface[i], modes[i], "eigenmodes");
        modes[i].Clear();
        freq [i].Clear();
      }
    }

    // Obtain correspondence map for initial links
    RegisteredSurface target, source;
    target.InputSurface(surface[0]);
    source.InputSurface(surface[1]);
    target.Initialize();
    source.Initialize();
    target.Update();
    source.Update();

    UniquePtr<PointCorrespondence> cmap(PointCorrespondence::New(ctype));
    cmap->FromTargetToSource(true);
    cmap->FromSourceToTarget(true);
    cmap->Parameter(cparam);
    cmap->Target(&target);
    cmap->Source(&source);
    for (size_t i = 0; i < cfeature_name.size(); ++i) {
      cmap->AddFeature(cfeature_name[i].c_str(), cfeature_weight[i]);
    }
    cmap->Initialize();
    cmap->Update();

    // Draw samples for which to add connectivity links
    Array<int> target_sample;
    Array<int> source_sample;
    using PointCorrespondenceUtils::SamplePoints;
    SamplePoints(&target, target_sample, target.NumberOfPoints() / 10);
    SamplePoints(&source, source_sample, source.NumberOfPoints() / 10);

    // Perform spectral analysis of joint graph Laplacian
    Matrix modes;
    Vector freq;
    if (ComputeJointEigenmodes(target.Surface(), target_sample,
                               source.Surface(), source_sample,
                               cmap.get(), k, modes, freq) < k) {
      FatalError("Error: Failed to find " << k << " eigenmodes of joint graph Laplacian!");
    }

    // Weight eigenmodes
    Vector w(k);
    double wsum = .0;
    for (int i = 0; i < k; ++i) {
      wsum += 1.0 / sqrt(freq(i));
    }
    for (int i = 0; i < k; ++i) {
      modes.ScaleCol(i, (1.0 / sqrt(freq(i))) / wsum);
    }

    // Add eigenmodes to output point data
    SetEigenmodes(surface[0], modes, 0,                       k, "joint_eigenmodes");
    SetEigenmodes(surface[1], modes, target.NumberOfPoints(), k, "joint_eigenmodes");

  } else {

    // Add individual eigenmodes to output point data
    for (int i = 0; i < ndatasets; ++i) {
      if (individual_analysis[i]) {
        SetEigenmodes(surface[i], modes[i], "eigenmodes");
        modes[i].Clear();
        freq [i].Clear();
      }
    }

  }

  // Match spectral components of input datasets with those of target dataset
  if (target_name) {
    if (verbose > 1) cout << "Reading target from " << target_name << endl;
    const bool exit_on_failure = false;
    vtkSmartPointer<vtkPointSet> target = ReadPointSet(target_name, exit_on_failure);
    if (!target || target->GetNumberOfPoints() == 0) {
      Warning("Failed to read target dataset or dataset has no points!"
              " Skipping matching of spectral components therefore.");
    } else {
      for (int i = 0; i < ndatasets; ++i) {
        // TODO: Match eigenmodes
      }
    }
  }

  // Write datasets together with their spectral components
  for (int i = 0; i < ndatasets; ++i) {
    Write(output_name[i], dataset[i], surface[i], fopt, as_points);
    if (verbose > 1) cout << "Wrote result dataset to " << output_name[i] << endl;
  }

  return 0;
}
