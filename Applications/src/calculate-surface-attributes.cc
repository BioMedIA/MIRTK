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
#include "mirtk/ImageSurfaceStatistics.h"
#include "mirtk/GradientImageFilter.h"
#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/FastCubicBSplineInterpolateImageFunction.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkPolyDataNormals.h"

using namespace mirtk;
using namespace mirtk::data::statistic;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Calculate attributes of input surface such as normals and curvature.\n";
  cout << "  If required, as in case of the curvature calculations, the input mesh\n";
  cout << "  is triangulated beforehand if it contains non-triangular faces.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input  surface mesh.\n";
  cout << "  output   Output surface mesh.\n";
  cout << "\n";
  cout << "Normals options:\n";
  cout << "  -normals, -point-normals   Surface point normals.\n";
  cout << "  -cell-normals              Surface cell normals.\n";
  cout << "  -[no]auto-orient           Enable/disable auto-orientation of normals. (default: on)\n";
  cout << "  -[no]splitting             Enable/disable splitting of sharp edges. (default: off)\n";
  cout << "  -[no]consistency           Enable/disable enforcement of vertex order consistency. (default: on)\n";
  cout << "\n";
  cout << "Curvature options:\n";
  cout << "  -k1k2 [<name>] [<name>]    Principal curvatures.\n";
  cout << "  -k1 [<name>]               Minimum curvature.\n";
  cout << "  -k2 [<name>]               Maximum curvature.\n";
  cout << "  -e1 [<name>]               Direction of minimum curvature.\n";
  cout << "  -e2 [<name>]               Direction of maximum curvature.\n";
  cout << "  -H [<name>]                Mean curvature:  H = .5 * (k1 + k2).\n";
  cout << "  -K [<name>]                Gauss curvature: K = k1 * k2.\n";
  cout << "  -C [<name>]                Curvedness:      C = sqrt(.5 * (k1^2 + k2^2)).\n";
  cout << "  -normalize                 Normalize curvature using volume of convex hull.\n";
  cout << "  -vtk-curvatures            Use vtkCurvatures when possible.\n";
  cout << "  -robust-curvatures         Do not use vtkCurvatures. Instead, estimate the curvature\n";
  cout << "                             tensor field and decompose it to obtain principle curvatures. (default)\n";
  cout << "\n";
  cout << "Local image options:\n";
  cout << "  -image <file>\n";
  cout << "      Input image file required by -gradient* and -patch* options.\n";
  cout << "  -gradient-normal [<name>]\n";
  cout << "      Compute image derivative in normal direction using cubic B-spline interpolation.\n";
  cout << "      The <name> of the output point data array is by default 'ImageGradientNormal'. (default: off)\n";
  cout << "  -gradient-angle [<name>]\n";
  cout << "      Compute cosine of angle made up by image gradient and normal vector using cubic\n";
  cout << "      B-spline interpolation for computing the image derivatives. The <name> of the output\n";
  cout << "      point data array is by default 'ImageGradientAngle' (default: off)\n";
  cout << "  -patch-name <name>\n";
  cout << "      Name of output point data array storing patch image statistics.\n";
  cout << "      (default: LocalImageStatistics)\n";
  cout << "  -patch-size <nx> [<ny> [<nz>]]\n";
  cout << "      Size of image patches. When only <nx> is given, an image patch of size\n";
  cout << "      nx = ny = nz is used. When only <nz> is omitted, a 2D patch is used.\n";
  cout << "  -patch-spacing <dx> [<dy> [<dz>]]\n";
  cout << "      Spacing between patch sample points. When only <dx> is given, an isotropic\n";
  cout << "      sampling in all three dimensions of <dx> is used. When only <dz> is omitted,\n";
  cout << "      a 2D patch spacing is used with dz=0.\n";
  cout << "  -patch-space image|world|tangent\n";
  cout << "      Coordinate system of patch. (default: tangent)\n";
  cout << "      - world:   Patch is aligned with world coordinate system.\n";
  cout << "      - image:   Patch is algined with image coordinate system.\n";
  cout << "      - tangent: Each patch is aligned with the coordinate system made up\n";
  cout << "                 by the normal vector and two orthonormal tangent vectors.\n";
  cout << "  -[no]patch-samples\n";
  cout << "      Whether to store individual intensities interpolated at patch sample points.\n";
  cout << "  -demean-patch\n";
  cout << "      Substract mean intensity from individual :option:`-patch-samples`.\n";
  cout << "  -whiten-patch\n";
  cout << "      Dividide individual :option:`-patch-samples` by standard deviation.\n";
  cout << "  -patch-min\n";
  cout << "      Append minimum patch intensity to output point data array.\n";
  cout << "  -patch-max\n";
  cout << "      Append maximum patch intensity to output point data array.\n";
  cout << "  -patch-min-abs\n";
  cout << "      Append minimum absolute patch intensity to output point data array.\n";
  cout << "  -patch-max-abs\n";
  cout << "      Append maximum absolute patch intensity to output point data array.\n";
  cout << "  -patch-mean\n";
  cout << "      Append mean patch intensity to output point data array.\n";
  cout << "  -patch-sigma\n";
  cout << "      Append standard deviation of patch intensities to output point data array.\n";
  cout << "\n";
  cout << "Smoothing options:\n";
  cout << "  -smooth-iterations [<niter>]\n";
  cout << "      Number of smoothing iterations.\n";
  cout << "  -smooth-weighting <name> [options]\n";
  cout << "      Smooth scalar attributes using the named weighting function:\n";
  cout << "      - 'Gaussian': Isotropic Gaussian smoothing kernel. (default)\n";
  cout << "        - Options: [<sigma>]\n";
  cout << "        - If sigma is not specified, it is automatically determined from the edges.\n";
  cout << "      - 'AnisotropicGaussian': Anisotropic Gaussian smoothing kernel.\n";
  cout << "        - Options: [<sigma>] [<sigma2>]\n";
  cout << "        - If sigma is not specified, it is automatically determined from the edges.\n";
  cout << "        - If sigma2 is specified, an anisotropic kernel with standard deviation\n";
  cout << "          sigma along the direction of minimum curvature, and sigma2 in the\n";
  cout << "          direction of maximum curvature is used.\n";
  cout << "        - If sigma2 is not specified, an isotropic Gaussian kernel used that is oriented\n";
  cout << "          and scaled along each local geometry axis using the curvature tensor.\n";
  cout << "      - 'InverseDistance': Inverse node distance.\n";
  cout << "        - Options: [<bias>]\n";
  cout << "        - If the bias is specified, the distance is estimated as 1/(dist+bias).\n";
  cout << "      - 'Combinatorial': Uniform node weights.\n";
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Compute surface normals when not available
vtkSmartPointer<vtkDataArray> Normals(vtkPolyData *surface)
{
  vtkSmartPointer<vtkDataArray> normals = surface->GetPointData()->GetNormals();
  if (normals == nullptr) {
    vtkNew<vtkPolyDataNormals> filter;
    SetVTKInput(filter, surface);
    filter->SplittingOff();
    filter->ComputePointNormalsOn();
    filter->ComputeCellNormalsOff();
    filter->ConsistencyOn();
    filter->AutoOrientNormalsOff();
    filter->FlipNormalsOff();
    filter->Update();
    normals = filter->GetOutput()->GetPointData()->GetNormals();
  }
  return normals;
}

// -----------------------------------------------------------------------------
/// Directional image gradient in normal direction
void EvaluateImageDerivative(vtkPolyData *surface, const RealImage &image, const char *name, bool normalize)
{
  const int npoints = static_cast<int>(surface->GetNumberOfPoints());

  Point   p;
  Vector3 n, g;
  Matrix  jac(1, 3);
  double  value;

  GenericFastCubicBSplineInterpolateImageFunction<RealImage> f;
  f.Input(&image);
  f.Initialize();

  vtkSmartPointer<vtkDataArray> output;
  output = NewVtkDataArray(VTK_FLOAT, npoints, 1, name);
  surface->GetPointData()->AddArray(output);

  vtkSmartPointer<vtkDataArray> normals = Normals(surface);
  for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
    surface->GetPoint(ptId, p);
    normals->GetTuple(ptId, n);
    f.WorldToImage(p);
    f.Jacobian3D(jac, p.x, p.y, p.z);
    g.x = jac(0, 0);
    g.y = jac(0, 1);
    g.z = jac(0, 2);
    f.ImageToWorld(g);
    if (normalize) g.Normalize();
    value = n.Dot(g);
    output->SetComponent(ptId, 0, value);
  }
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

  ImageSurfaceStatistics image_stats;
  SharedPtr<InterpolateImageFunction> image_func;
  image_func.reset(InterpolateImageFunction::New(Interpolation_Linear));
  const char *image_name             = nullptr;
  const char *cosine_image_gradient  = nullptr;
  const char *normal_image_gradient  = nullptr;
  bool        calc_image_gradient    = false;

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
      if (smooth_iterations == 0) smooth_iterations = 1;
      PARSE_ARGUMENT(weighting);
      if (weighting ==  MeshSmoothing::InverseDistance ||
          weighting ==  MeshSmoothing::Gaussian ||
          weighting ==  MeshSmoothing::AnisotropicGaussian) {
        if (HAS_ARGUMENT) PARSE_ARGUMENT(smooth_sigma);
        if (weighting ==  MeshSmoothing::AnisotropicGaussian) {
          smooth_along_tensor = true;
          if (HAS_ARGUMENT) {
            smooth_along_tensor = false;
            PARSE_ARGUMENT(smooth_sigma2);
          }
        }
      }
    }
    else if (OPTION("-image")) {
      image_name = ARGUMENT;
    }
    else if (OPTION("-gradient-normal")) {
      calc_image_gradient = true;
      if (HAS_ARGUMENT) normal_image_gradient = ARGUMENT;
      else              normal_image_gradient = "ImageGradientNormal";
    }
    else if (OPTION("-gradient-angle")) {
      calc_image_gradient = true;
      if (HAS_ARGUMENT) cosine_image_gradient = ARGUMENT;
      else              cosine_image_gradient = "ImageGradientAngle";
    }
    else if (OPTION("-patch-name")) {
      image_stats.ArrayName(ARGUMENT);
    }
    else if (OPTION("-patch-size")) {
      int nx, ny, nz;
      PARSE_ARGUMENT(nx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(ny);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(nz);
        else nz = 1;
      } else {
        ny = nz = nx;
      }
      image_stats.PatchSize(make_int3(nx, ny, nz));
    }
    else if (OPTION("-patch-spacing")) {
      double dx, dy, dz;
      PARSE_ARGUMENT(dx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(dy);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(dz);
        else dz = 0.;
      } else {
        dy = dz = dx;
      }
      image_stats.PatchSpacing(make_double3(dx, dy, dz));
    }
    else if (OPTION("-patch-space")) {
      const char *arg = ARGUMENT;
      if (strcmp(arg, "image") == 0) {
        image_stats.PatchSpace(ImageSurfaceStatistics::ImageSpace);
      } else if (strcmp(arg, "world") == 0) {
        image_stats.PatchSpace(ImageSurfaceStatistics::WorldSpace);
      } else if (strcmp(arg, "tangent") == 0) {
        image_stats.PatchSpace(ImageSurfaceStatistics::TangentSpace);
      } else {
        FatalError("Invalid -patch-space argument: " << arg);
      }
    }
    else if (OPTION("-patch-samples")) {
      image_stats.PatchSamples(true);
    }
    else if (OPTION("-nopatch-samples")) {
      image_stats.PatchSamples(false);
    }
    else if (OPTION("-demean-patch")) {
      image_stats.DemeanSamples(true);
    }
    else if (OPTION("-whiten-patch")) {
      image_stats.WhitenSamples(true);
    }
    else if (OPTION("-patch-mean")) {
      image_stats.Statistics().push_back(NewShared<Mean>());
    }
    else if (OPTION("-patch-sigma")) {
      image_stats.Statistics().push_back(NewShared<StDev>());
    }
    else if (OPTION("-patch-min")) {
      image_stats.Statistics().push_back(NewShared<Min>());
    }
    else if (OPTION("-patch-min-abs")) {
      image_stats.Statistics().push_back(NewShared<MinAbs>());
    }
    else if (OPTION("-patch-max")) {
      image_stats.Statistics().push_back(NewShared<Max>());
    }
    else if (OPTION("-patch-max-abs")) {
      image_stats.Statistics().push_back(NewShared<MaxAbs>());
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (calc_image_gradient && !image_name) {
    FatalError("The -gradient* options require an input -image!");
  }

  bool calc_normals     = (point_normals || cell_normals);
  bool calc_image_stats = (image_name && (image_stats.PatchSamples() || !image_stats.Statistics().empty()));
  if (curvatures == 0 && !calc_normals && !calc_image_stats && !calc_image_gradient) {
    point_normals = true;
    curvatures    = SurfaceCurvature::Scalars;
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

  // Read input image
  RealImage image;
  if (calc_image_stats || calc_image_gradient) {
    InitializeIOLibrary();
    image.Read(image_name);
    if (calc_image_stats) {
      image_func->Input(&image);
      image_func->Initialize();
    }
  }

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

    // Smooth calculated curvatures
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
      if (e1)         smoother.SmoothArray(e1_name, vtkDataSetAttributes::VECTORS);
      if (e2)         smoother.SmoothArray(e2_name, vtkDataSetAttributes::VECTORS);
      smoother.NumberOfIterations(smooth_iterations);
      smoother.Sigma(-smooth_sigma); // negative: multiple of avg. edge length
      smoother.Weighting(weighting);
      if (weighting == MeshSmoothing::AnisotropicGaussian) {
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

  // Calculate local image patch statistics
  if (calc_image_stats) {
    if (verbose) cout << "Calculating local image statistics...", cout.flush();
    image_stats.Image(image_func);
    image_stats.Input(surface);
    image_stats.Run();
    surface = image_stats.Output();
    if (verbose) cout << " done" << endl;

    // Smooth local image statistics
    if (smooth_iterations) {
      if (verbose) cout << "Smoothing local image statistics...", cout.flush();
      const char *name = image_stats.ArrayName().c_str();

      MeshSmoothing smoother;
      smoother.Input(surface);
      smoother.SmoothPointsOff();
      smoother.SmoothArray(name);
      smoother.NumberOfIterations(smooth_iterations);
      smoother.Sigma(-smooth_sigma); // negative: multiple of avg. edge length
      smoother.Weighting(weighting);
      if (weighting == MeshSmoothing::AnisotropicGaussian){
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
      surface->GetPointData()->GetArray(name)->DeepCopy(pd->GetArray(name));

      if (verbose) cout << " done" << endl;
    }
  }

  // Calculate image gradient
  if (calc_image_gradient) {
    if (verbose) cout << "Evaluating image gradient...", cout.flush();
    if (normal_image_gradient) {
      EvaluateImageDerivative(surface, image, normal_image_gradient, false);
    }
    if (cosine_image_gradient) {
      EvaluateImageDerivative(surface, image, cosine_image_gradient, true);
    }
    if (verbose) cout << " done" << endl;
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
