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

#include "mirtk/ImageSurfaceStatistics.h"

#include "mirtk/PointSetUtils.h"
#include "mirtk/Matrix3x3.h"

#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkPolyDataNormals.h"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace ImageSurfaceStatisticsUtils {


// -----------------------------------------------------------------------------
/// Base class of ImageSurfaceStatistics::Execute function body
struct CalculatePatchStatistics
{
  vtkPoints                      *_Points;
  const InterpolateImageFunction *_Image;
  int3                            _Size;
  int                             _N;
  Array<const data::Statistic *>  _Statistics;
  vtkDataArray                   *_Output;
  bool                            _Samples;
  bool                            _Demean;
  bool                            _Whiten;

protected:

  /// Calculate patch values and statistics at given surface point
  void Execute(vtkIdType ptId, Array<double> &patch_mem, Array<double> &stats_mem,
               const Vector3 &dx, const Vector3 &dy, const Vector3 &dz) const
  {
    Point o, p;
    patch_mem.resize(_N);
    double * const patch = patch_mem.data();

    // Sample patch values
    _Points->GetPoint(ptId, o);
    _Image->WorldToImage(o._x, o._y, o._z);
    o -= (_Size.x / 2.) * dx;
    o -= (_Size.y / 2.) * dy;
    o -= (_Size.z / 2.) * dz;
    double *v = patch;
    for (int k = 0; k < _Size.z; ++k)
    for (int j = 0; j < _Size.y; ++j) {
      p = o + j * dy + k * dz;
      for (int i = 0; i < _Size.x; ++i, ++v, p += dx) {
        (*v) = _Image->Evaluate(p._x, p._y, p._z);
      }
    }

    // Evaluate statistics (**before** subtracting mean or dividing by standard deviation)
    if (!_Statistics.empty()) {
      int f = (_Samples ? _N : 0);
      Array<double> &stats = stats_mem;
      for (size_t i = 0; i < _Statistics.size(); ++i) {
        _Statistics[i]->Evaluate(stats, _N, patch);
        for (size_t j = 0; j < stats.size(); ++j, ++f) {
          _Output->SetComponent(ptId, f, stats[j]);
        }
      }
    }

    // Store patch samples (**after** evaluating statistics of unmodified samples)
    if (_Samples) {
      if (_Demean || _Whiten) {
        double mean, sigma;
        data::statistic::NormalDistribution::Calculate(mean, sigma, _N, patch);
        if (_Demean && _Whiten) {
          for (int i = 0; i < _N; ++i) {
            patch[i] = (patch[i] - mean) / sigma;
          }
        } else if (_Demean) {
          for (int i = 0; i < _N; ++i) {
            patch[i] -= mean;
          }
        } else {
          for (int i = 0; i < _N; ++i) {
            patch[i] = mean + (patch[i] - mean) / sigma;
          }
        }
      }
      for (int i = 0; i < _N; ++i) {
        _Output->SetComponent(ptId, i, patch[i]);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate image statistics for patch aligned with global coordinate axes
struct CalculateImagePatchStatistics : public CalculatePatchStatistics
{
  Vector3 _DirX, _DirY, _DirZ;
  void operator ()(const blocked_range<vtkIdType> ptIds) const
  {
    Array<double> patch_mem(_N), stats_mem;
    stats_mem.reserve(10);
    for (vtkIdType ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      Execute(ptId, patch_mem, stats_mem, _DirX, _DirY, _DirZ);
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate image statistics for patch aligned with local tangent vectors
struct CalculateTangentPatchStatistics : public CalculatePatchStatistics
{
  vtkDataArray *_Normals;
  Matrix3x3     _Rotation;
  double3       _Scaling;

  void operator ()(const blocked_range<vtkIdType> ptIds) const
  {
    Vector3 dx, dy, dz;
    Array<double> patch_mem(_N), stats_mem;
    stats_mem.reserve(10);
    for (vtkIdType ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Normals->GetTuple(ptId, dx);
      dx = _Rotation * dx;
      ComputeTangents(dx, dy, dz);
      dx *= _Scaling.x;
      dy *= _Scaling.y;
      dz *= _Scaling.z;
      Execute(ptId, patch_mem, stats_mem, dx, dy, dz);
    }
  }
};


} // ImageSurfaceStatisticsUtils
using namespace ImageSurfaceStatisticsUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ImageSurfaceStatistics::CopyAttributes(const ImageSurfaceStatistics &other)
{
  _Image         = other._Image;
  _ArrayName     = other._ArrayName;
  _PatchSpace    = other._PatchSpace;
  _PatchSize     = other._PatchSize;
  _PatchSpacing  = other._PatchSpacing;
  _PatchSamples  = other._PatchSamples;
  _DemeanSamples = other._DemeanSamples;
  _WhitenSamples = other._WhitenSamples;
}

// -----------------------------------------------------------------------------
ImageSurfaceStatistics::ImageSurfaceStatistics()
:
  _PatchSpace(TangentSpace),
  _PatchSize(make_int3(5)),
  _PatchSpacing(make_double3(1.)),
  _PatchSamples(false),
  _DemeanSamples(false),
  _WhitenSamples(false)
{
}

// -----------------------------------------------------------------------------
ImageSurfaceStatistics::ImageSurfaceStatistics(const ImageSurfaceStatistics &other)
:
  SurfaceFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
ImageSurfaceStatistics &ImageSurfaceStatistics::operator =(const ImageSurfaceStatistics &other)
{
  if (this != &other) {
    SurfaceFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
ImageSurfaceStatistics::~ImageSurfaceStatistics()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class Body>
void ImageSurfaceStatistics::ExecuteInParallel(Body &body)
{
  body._Points  = _Output->GetPoints();
  body._Image   = _Image.get();
  body._Size    = _PatchSize;
  body._N       = _PatchSize.x * _PatchSize.y * _PatchSize.z;
  body._Samples = _PatchSamples;
  body._Demean  = _DemeanSamples;
  body._Whiten  = _WhitenSamples;

  body._Statistics.resize(_Statistics.size());
  for (size_t i = 0; i < _Statistics.size(); ++i) {
    body._Statistics[i] = _Statistics[i].get();
  }

  int nfeatures = static_cast<int>(_Statistics.size());
  if (_PatchSamples) nfeatures += body._N;
  if (nfeatures <= 0) {
    Throw(ERR_LogicError, "Execute", "No surface patch samples and statistics to evaluate");
  }

  vtkSmartPointer<vtkDataArray> output;
  output = NewArray("LocalImageStatistics", nfeatures);
  if (!_ArrayName.empty()) output->SetName(_ArrayName.c_str());
  body._Output = output;

  if (body._Samples) {
    char name[64];
    for (int i = 0; i < body._N; ++i) {
      snprintf(name, 64, "Sample %d", i+1);
      output->SetComponentName(i, name);
    }
  }
  int f = (_PatchSamples ? body._N : 0);
  for (size_t i = 0; i < _Statistics.size(); ++i) {
    const Array<string> &names = _Statistics[i]->Names();
    for (size_t j = 0; j < names.size(); ++j, ++f) {
      output->SetComponentName(f, names[j].c_str());
    }
  }

  parallel_for(blocked_range<vtkIdType>(0, _Output->GetNumberOfPoints()), body);
  _Output->GetPointData()->AddArray(output);
}

// -----------------------------------------------------------------------------
void ImageSurfaceStatistics::Execute()
{
  switch (_PatchSpace) {

    // Patches aligned with image coordinate axes
    case ImageSpace: {
      CalculateImagePatchStatistics body;
      body._DirX = Vector3(_PatchSpacing.x / _Image->XSize(), 0., 0.);
      body._DirY = Vector3(0., _PatchSpacing.y / _Image->YSize(), 0.);
      body._DirZ = Vector3(0., 0., _PatchSpacing.z / _Image->ZSize());
      ExecuteInParallel(body);
    } break;

    // Patches aligned with world coordinate axes
    case WorldSpace: {
      Matrix A(3, 3);
      A(0, 0) = _PatchSpacing.x / _Image->XSize();
      A(1, 1) = _PatchSpacing.y / _Image->YSize();
      A(2, 2) = _PatchSpacing.z / _Image->ZSize();
      A = _Image->Attributes().GetWorldToImageOrientation() * A;
      CalculateImagePatchStatistics body;
      body._DirX = Vector3(A.Col(0));
      body._DirY = Vector3(A.Col(1));
      body._DirZ = Vector3(A.Col(2));
      ExecuteInParallel(body);
    } break;

    // Patches aligned with normal and tangent vectors
    case TangentSpace: {
      vtkSmartPointer<vtkDataArray> normals = _Input->GetPointData()->GetNormals();
      if (!normals) {
        vtkNew<vtkPolyDataNormals> calc_normals;
        SetVTKInput(calc_normals, _Input);
        calc_normals->SplittingOff();
        calc_normals->Update();
        normals = calc_normals->GetOutput()->GetPointData()->GetNormals();
      }
      const Matrix R = _Image->Attributes().GetWorldToImageOrientation();
      CalculateTangentPatchStatistics body;
      body._Normals   = normals;
      body._Rotation  = Matrix3x3(R(0, 0), R(0, 1), R(0, 2),
                                  R(1, 0), R(1, 1), R(1, 2),
                                  R(2, 0), R(2, 1), R(2, 2));
      body._Scaling.x = _PatchSpacing.x / _Image->XSize();
      body._Scaling.y = _PatchSpacing.y / _Image->YSize();
      body._Scaling.z = _PatchSpacing.z / _Image->ZSize();
      ExecuteInParallel(body);
    } break;

    default: {
      Throw(ERR_LogicError, __func__, "Invalid patch space enumeration value: ", _PatchSpace);
    } break;
  }
}


} // namespace mirtk
