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
#include "mirtk/SurfaceCurvature.h"
#include "mirtk/MeshSmoothing.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkPolyDataNormals.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "Calculate attributes of input surface such as normals and curvature." << endl;
  cout << "If required, as in case of the curvature calculations, the input mesh" << endl;
  cout << "is triangulated beforehand if it contains non-triangular faces." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input surface mesh." << endl;
  cout << "  output   Output surface mesh." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -normals, -point-normals        Surface point normals." << endl;
  cout << "  -cell-normals                   Surface cell normals." << endl;
  cout << "  -[no]auto-orient                Enable/disable auto-orientation of normals. (default: on)" << endl;
  cout << "  -[no]splitting                  Enable/disable splitting of sharp edges. (default: off)" << endl;
  cout << "  -[no]consistency                Enable/disable enforcement of vertex order consistency. (default: on)" << endl;
  cout << endl;
  cout << "Curvature output options:" << endl;
  cout << "  -H [<name>]                     Mean curvature." << endl;
  cout << "  -K [<name>]                     Gauss curvature." << endl;
  cout << "  -C [<name>]                     Curvedness." << endl;
  cout << "  -k1 [<name>]                    Minimum curvature." << endl;
  cout << "  -k2 [<name>]                    Maximum curvature." << endl;
  cout << "  -k1k2 [<name>] [<name>]         Principal curvatures." << endl;
  cout << "  -e1 [<name>]                    Direction of minimum curvature." << endl;
  cout << "  -e2 [<name>]                    Direction of maximum curvature." << endl;
  cout << "  -normalize                      Normalize curvature using volume of convex hull." << endl;
  cout << "  -vtk-curvatures                 Use vtkCurvatures when possible." << endl;
  cout << "  -robust-curvatures              Do not use vtkCurvatures. Instead, estimate the curvature" << endl;
  cout << "                                  tensor field and decompose it to obtain principle curvatures. (default)" << endl;
  cout << endl;
  cout << "  -smooth-iterations [<niter>]    Number of smoothing iterations" << endl;
  cout << "  -smooth-weighting <name> [options]" << endl;
  cout << "      Smooth calculated scalar curvature measures according to the weighting specified in <name>:" <<endl;
  cout << "      - 'Gaussian': using a Gaussian smoothing kernel (default)." << endl;
  cout << "                  Options: [<sigma>]" << endl;
  cout << "                  If sigma is not specified, it is automatically determined from the edges." << endl;
  cout << "      - 'AnisotropicGaussian': using an anisotropic Gaussian smoothing kernel." << endl;
  cout << "                  Options: [<sigma>] [<sigma2>]" << endl;
  cout << "                  If sigma is not specified, it is automatically determined from the edges." << endl;
  cout << "                  If sigma2 is specified, an anisotropic kernel with standard deviation" << endl;
  cout << "                  sigma along the direction of minimum curvature, and sigma2 in the" << endl;
  cout << "                  direction of maximum curvature is used." << endl;
  cout << "                  If sigma2 is not specified, an isotropic Gaussian kernel used that is oriented" << endl;
  cout << "                  and scaled along each local geometry axis using the curvature tensor." << endl;
  cout << "      - 'InverseDistance': using the inverse node distance." << endl;
  cout << "                  Options: [<bias>]" << endl;
  cout << "                  If the bias is specified the distance is estimated as 1/(dist+bias)" << endl;
  cout << "      - 'Combinatorial': using uniform node weights." << endl;


  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  const char *kmin_name       = SurfaceCurvature::MINIMUM;
  const char *kmax_name       = SurfaceCurvature::MAXIMUM;
  const char *gauss_name      = SurfaceCurvature::GAUSS;
  const char *mean_name       = SurfaceCurvature::MEAN;
  const char *curvedness_name = SurfaceCurvature::CURVEDNESS;
  const char *e1_name         = SurfaceCurvature::MINIMUM_DIRECTION;
  const char *e2_name         = SurfaceCurvature::MAXIMUM_DIRECTION;
  const char *tensor_name     = SurfaceCurvature::TENSOR;
  const char *inverse_name    = SurfaceCurvature::INVERSE_TENSOR;

  bool   point_normals        = false;
  bool   cell_normals         = false;
  bool   auto_orient_normals  = true;
  bool   splitting            = false;
  bool   consistency          = true;
  bool   use_vtkCurvatures    = false;
  int    curvatures           = 0;
  int    tensor_averaging     = 3;
  bool   normalize            = false;
  int    smooth_iterations    = 0;
  double smooth_sigma         = .0;
  double smooth_sigma2        = .0;
  bool   smooth_along_tensor  = false;
  MeshSmoothing::WeightFunction weighting = MeshSmoothing::Default;

  for (ALL_OPTIONS) {
    if (OPTION("-point-normals") || OPTION("-normals")) point_normals = true;
    else if (OPTION("-cell-normals")) cell_normals = true;
    else if (OPTION("-auto-orient"))   auto_orient_normals = true;
    else if (OPTION("-noauto-orient")) auto_orient_normals = false;
    else if (OPTION("-consistency"))   consistency = true;
    else if (OPTION("-noconsistency")) consistency = false;
    else if (OPTION("-k1")) {
      curvatures |= SurfaceCurvature::Minimum;
      if (HAS_ARGUMENT) kmin_name = ARGUMENT;
    }
    else if (OPTION("-k2")) {
      curvatures |= SurfaceCurvature::Maximum;
      if (HAS_ARGUMENT) kmax_name = ARGUMENT;
    }
    else if (OPTION("-k1k2")) {
      curvatures |= (SurfaceCurvature::Minimum | SurfaceCurvature::Maximum);
      if (HAS_ARGUMENT) {
        kmin_name = ARGUMENT;
        kmax_name = ARGUMENT;
      }
    }
    else if (OPTION("-H")) {
      curvatures |= SurfaceCurvature::Mean;
      if (HAS_ARGUMENT) mean_name = ARGUMENT;
    }
    else if (OPTION("-K")) {
      curvatures |= SurfaceCurvature::Gauss;
      if (HAS_ARGUMENT) gauss_name = ARGUMENT;
    }
    else if (OPTION("-C")) {
      curvatures |= SurfaceCurvature::Curvedness;
      if (HAS_ARGUMENT) curvedness_name = ARGUMENT;
    }
    else if (OPTION("-n") || OPTION("-normal")) {
      curvatures |= SurfaceCurvature::Normal;
    }
    else if (OPTION("-e1")) {
      curvatures |= SurfaceCurvature::MinimumDirection;
      if (HAS_ARGUMENT) e1_name = ARGUMENT;
    }
    else if (OPTION("-e2")) {
      curvatures |= SurfaceCurvature::MaximumDirection;
      if (HAS_ARGUMENT) e2_name = ARGUMENT;
    }
    else if (OPTION("-tensor")) {
      curvatures |= SurfaceCurvature::Tensor;
      if (HAS_ARGUMENT) tensor_name = ARGUMENT;
    }
    else if (OPTION("-inverse-tensor")) {
      curvatures |= SurfaceCurvature::InverseTensor;
      if (HAS_ARGUMENT) inverse_name = ARGUMENT;
    }
    else if (OPTION("-tensor-averaging")) {
      PARSE_ARGUMENT(tensor_averaging);
    }
    else if (OPTION("-normalize")) normalize = true;
    else if (OPTION("-vtk-curvatures"))    use_vtkCurvatures = true;
    else if (OPTION("-robust-curvatures")) use_vtkCurvatures = false;
    else if (OPTION("-smooth-iterations")) PARSE_ARGUMENT(smooth_iterations);
    else if (OPTION("-smooth-weighting")){
      if(smooth_iterations == 0) smooth_iterations = 1;

      if(!PARSE_ARGUMENT(weighting)){
        FatalError("Invalid -smooth-weighting <name> argument specified");
      }

      if( weighting ==  MeshSmoothing::InverseDistance || weighting ==  MeshSmoothing::Gaussian || weighting ==  MeshSmoothing::AnisotropicGaussian){
        if (HAS_ARGUMENT) PARSE_ARGUMENT(smooth_sigma);

        if (weighting ==  MeshSmoothing::AnisotropicGaussian){
          smooth_along_tensor = true;
          if(HAS_ARGUMENT) {
            smooth_along_tensor = false;
            PARSE_ARGUMENT(smooth_sigma2);
          }
        }

      }
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (curvatures == 0 && !point_normals && !cell_normals) {
    point_normals = true;
    curvatures = SurfaceCurvature::Scalars;
  }

  int curvature_type = curvatures;
  if (smooth_along_tensor) {
    curvature_type |= SurfaceCurvature::Tensor;
  } else if (smooth_sigma2 != .0) {
    curvature_type |= SurfaceCurvature::MinimumDirection;
    curvature_type |= SurfaceCurvature::MaximumDirection;
  }

  // Read input surface
  vtkSmartPointer<vtkPolyData> input = ReadPolyData(input_name);

  // Triangulate input surface
  vtkSmartPointer<vtkPolyData> surface;
  if (curvature_type != 0) surface = Triangulate(input);
  else                     surface = input;

  // Calculate normals
  if (point_normals || cell_normals) {
    if (verbose) cout << "Calculating surface normals...", cout.flush();
    vtkNew<vtkPolyDataNormals> filter;
    SetVTKInput(filter, surface);
    filter->SetSplitting(splitting);
    filter->SetComputePointNormals(point_normals);
    filter->SetComputeCellNormals(cell_normals);
    filter->SetConsistency(consistency);
    filter->SetAutoOrientNormals(auto_orient_normals);
    filter->FlipNormalsOff();
    filter->Update();
    surface = filter->GetOutput();
    if (verbose) cout << " done" << endl;
  }
 
  // Calculate curvatures
  if (curvature_type != 0) {
    if (verbose) cout << "Calculating surface curvature measure(s)...", cout.flush();

    // Compute curvature
    SurfaceCurvature curvature;
    curvature.Input(surface);
    curvature.CurvatureType(curvature_type);
    curvature.VtkCurvatures(use_vtkCurvatures);
    curvature.TensorAveraging(tensor_averaging);
    curvature.Normalize(normalize);
    curvature.Run();
    surface = curvature.Output();

    // Get output arrays
    vtkPointData *pd         = surface->GetPointData();
    vtkDataArray *kmin       = pd->GetArray(SurfaceCurvature::MINIMUM);
    vtkDataArray *kmax       = pd->GetArray(SurfaceCurvature::MAXIMUM);
    vtkDataArray *mean       = pd->GetArray(SurfaceCurvature::MEAN);
    vtkDataArray *gauss      = pd->GetArray(SurfaceCurvature::GAUSS);
    vtkDataArray *curvedness = pd->GetArray(SurfaceCurvature::CURVEDNESS);
    vtkDataArray *e1         = pd->GetArray(SurfaceCurvature::MINIMUM_DIRECTION);
    vtkDataArray *e2         = pd->GetArray(SurfaceCurvature::MAXIMUM_DIRECTION);
    vtkDataArray *tensor     = pd->GetArray(SurfaceCurvature::TENSOR);
    vtkDataArray *inverse    = pd->GetArray(SurfaceCurvature::INVERSE_TENSOR);

    // Rename output arrays
    if (kmin)       kmin      ->SetName(kmin_name);
    if (kmax)       kmax      ->SetName(kmax_name);
    if (mean)       mean      ->SetName(mean_name);
    if (gauss)      gauss     ->SetName(gauss_name);
    if (curvedness) curvedness->SetName(curvedness_name);
    if (e1)         e1        ->SetName(e1_name);
    if (e2)         e2        ->SetName(e2_name);
    if (tensor)     tensor    ->SetName(tensor_name);
    if (inverse)    inverse   ->SetName(inverse_name);

    if (verbose) cout << " done" << endl;

    // Smooth calculated attributes
    if (smooth_iterations) {
      if (verbose) cout << "Smoothing scalar curvature measures...", cout.flush();

      MeshSmoothing smoother;
      smoother.Input(surface);
      smoother.SmoothPointsOff();
      if (kmin)       smoother.SmoothArray(kmin_name);
      if (kmax)       smoother.SmoothArray(kmax_name);
      if (mean)       smoother.SmoothArray(mean_name);
      if (gauss)      smoother.SmoothArray(gauss_name);
      if (curvedness) smoother.SmoothArray(curvedness_name);
      if (e1)         smoother.SmoothArray(e1_name);
      if (e2)         smoother.SmoothArray(e2_name);
      smoother.NumberOfIterations(smooth_iterations);
      smoother.Sigma(-smooth_sigma); // negative: multiple of avg. edge length
      smoother.Weighting(weighting);
      if (weighting ==  MeshSmoothing::AnisotropicGaussian){
        if (smooth_along_tensor) {
          smoother.GeometryTensorName(tensor_name);
        } else {
          smoother.MinimumDirectionName(e2_name);
          smoother.MaximumDirectionName(e1_name);
        }
        smoother.MaximumDirectionSigma(-smooth_sigma2);
      }
      smoother.Run();
      vtkPointData *pd = smoother.Output()->GetPointData();
      if (kmin) kmin->DeepCopy(pd->GetArray(kmin_name));
      if (kmax) kmax->DeepCopy(pd->GetArray(kmax_name));
      if (mean) mean->DeepCopy(pd->GetArray(mean_name));
      if (gauss) gauss->DeepCopy(pd->GetArray(gauss_name));
      if (curvedness) curvedness->DeepCopy(pd->GetArray(curvedness_name));
      if (e1) e1->DeepCopy(pd->GetArray(e1_name));
      if (e2) e2->DeepCopy(pd->GetArray(e2_name));

      if (verbose) cout << " done" << endl;
    }
  }

  // Remove not requested output arrays which were used for anisotropic smoothing
  if ((curvatures & SurfaceCurvature::Tensor) == 0) {
    surface->GetPointData()->RemoveArray(tensor_name);
  }
  if ((curvatures & SurfaceCurvature::MinimumDirection) == 0) {
    surface->GetPointData()->RemoveArray(e1_name);
  }
  if ((curvatures & SurfaceCurvature::MaximumDirection) == 0) {
    surface->GetPointData()->RemoveArray(e2_name);
  }

  // Write output surface
  if (verbose) cout << "Writing output surface to file " << output_name << "...", cout.flush();
  if (!WritePolyData(output_name, surface)) {
    if (verbose) cout << " failed" << endl;
    FatalError("Failed to write output surface to file " << output_name);
  }
  if (verbose) cout << " done" << endl;

  return 0;
}
