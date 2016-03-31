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

#include "mirtk/ImplicitSurfaceUtils.h"

#include "mirtk/Vtk.h"
#include "mirtk/Math.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/Resampling.h"
#include "mirtk/LinearInterpolateImageFunction.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkMarchingCubes.h"
#include "vtkStructuredPoints.h"

#include <algorithm>


namespace mirtk { namespace ImplicitSurfaceUtils {


// =============================================================================
// Contouring
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> Isosurface(const DistanceImage &dmap, double offset,
                                        double blurring, bool isotropic, bool close,
                                        bool normals, bool gradients)
{
  typedef DistanceImage::VoxelType VoxelType;
  const DistanceImage *distance_image = &dmap;

  // Blur distance image
  if (blurring > .0) {
    DistanceImage *blurred_image = new DistanceImage(*distance_image);
    GaussianBlurring<VoxelType> blur(blurring);
    blur.Input (distance_image);
    blur.Output(blurred_image);
    blur.Run();
    if (distance_image != &dmap) delete distance_image;
    distance_image = blurred_image;
  }

  // Resample to isotropic voxel size
  bool is_isotropic = fequal(dmap.GetXSize(), dmap.GetYSize()) &&
                      fequal(dmap.GetXSize(), dmap.GetZSize());
  if (isotropic && !is_isotropic) {
    double ds = min(min(dmap.GetXSize(), dmap.GetYSize()), dmap.GetZSize());
    DistanceImage *isotropic_image = new DistanceImage(dmap);
    Resampling<VoxelType> resampling(ds, ds, ds);
    LinearInterpolateImageFunction interpolator;
    resampling.Input (distance_image);
    resampling.Output(isotropic_image);
    resampling.Interpolator(&interpolator);
    resampling.Run();
    if (distance_image != &dmap) delete distance_image;
    distance_image = isotropic_image;
  }

  int mx = 0, my = 0, mz = 0; // width of margin
  int nx = distance_image->X();
  int ny = distance_image->Y();
  int nz = distance_image->Z();
  if (close) {
    mx = 1, my = 1, mz = 1;
    nx += mx + mx, ny += my + my, nz += mz + mz;
  }

  // Convert image to VTK (in voxel units due to missing orientation)
  vtkSmartPointer<vtkStructuredPoints> vtkimage;
  vtkimage = vtkSmartPointer<vtkStructuredPoints>::New();
  vtkimage->SetOrigin(-mx, -my, -mz);
  vtkimage->SetDimensions(nx, ny, nz);
  vtkimage->SetSpacing(1.0, 1.0, 1.0);
#if VTK_MAJOR_VERSION >= 6
  vtkimage->AllocateScalars(distance_image->ImageToVTKScalarType(), 1);
#else
  vtkimage->SetScalarType(distance_image->ImageToVTKScalarType());
  vtkimage->AllocateScalars();
#endif

  const double boundary_value = offset + 10.0;
  const VoxelType *d = distance_image->Data();
  VoxelType       *o = reinterpret_cast<VoxelType *>(vtkimage->GetScalarPointer());

  const int i1 = mx, i2 = nx - mx - 1;
  const int j1 = my, j2 = ny - my - 1;
  const int k1 = mz, k2 = nz - mz - 1;

  for (int k = 0; k < nz; ++k)
  for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i, ++o) {
    if (i1 <= i && i <= i2 &&
        j1 <= j && j <= j2 &&
        k1 <= k && k <= k2) {
	    *o = *d, ++d;
    } else {
      *o = boundary_value;
    }
  }

  // Extract isosurface
  vtkNew<vtkMarchingCubes> mcubes;
  mcubes->SetComputeNormals(normals);
  mcubes->SetComputeGradients(gradients);
  mcubes->SetNumberOfContours(1);
  mcubes->SetValue(0, offset);
  SetVTKInput(mcubes, vtkimage);
  mcubes->Update();

  // Map points from voxel indices to world coordinates
  double p[3];
  vtkPoints *points = mcubes->GetOutput()->GetPoints();
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    distance_image->ImageToWorld(p[0], p[1], p[2]);
    points->SetPoint(ptId, p);
  }

  mcubes->GetOutput()->GetPointData()->SetScalars(nullptr);

  return mcubes->GetOutput();
}


} } // namespace mirtk::ImplicitSurfaceUtils
