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

#include "mirtk/SurfaceDistance.h"

#include "mirtk/Assert.h"
#include "mirtk/Array.h"
#include "mirtk/Vector3D.h"
#include "mirtk/PointSetIO.h"

#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkVertex.h"

#include <cstdio>


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
SurfaceDistance::SurfaceDistance(const char *name, double weight)
:
  PointSetDistance(name, weight)
{
  _ParameterPrefix.push_back("Surface distance ");
}

// -----------------------------------------------------------------------------
SurfaceDistance::SurfaceDistance(const SurfaceDistance &other, int m, int n)
:
  PointSetDistance(other,
    other._GradientWrtTarget && other._Target ? other._Target->NumberOfSurfacePoints() : 0,
    other._GradientWrtSource && other._Source ? other._Source->NumberOfSurfacePoints() : 0
  )
{
}

// -----------------------------------------------------------------------------
void SurfaceDistance::CopyAttributes(const SurfaceDistance &other, int m, int n)
{
  PointSetDistance::CopyAttributes(other,
    m == -1 ? (other._GradientWrtTarget && other._Target ? other._Target->NumberOfSurfacePoints() : 0) : m,
    n == -1 ? (other._GradientWrtSource && other._Source ? other._Source->NumberOfSurfacePoints() : 0) : n
  );
}

// -----------------------------------------------------------------------------
SurfaceDistance &SurfaceDistance::operator =(const SurfaceDistance &other)
{
  CopyAttributes(other);
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceDistance::~SurfaceDistance()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceDistance::Initialize()
{
  // Check inputs
  if (_Target == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Target dataset is NULL" << endl;
    exit(1);
  }
  if (_Source == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Source dataset is NULL" << endl;
    exit(1);
  }

  // Initialize base class
  PointSetDistance::Initialize(_Target->NumberOfSurfacePoints(),
                               _Source->NumberOfSurfacePoints());
}

// -----------------------------------------------------------------------------
void SurfaceDistance::Reinitialize()
{
  // Reinitialize base class
  PointSetDistance::Reinitialize(_Target->NumberOfSurfacePoints(),
                                 _Source->NumberOfSurfacePoints());
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceDistance
::ParametricGradient(const RegisteredPointSet *wrt_pset,
                     const Vector3D<double>   *np_gradient,
                     double                   *gradient,
                     double                    weight)
{
  const class Transformation * const T = wrt_pset->Transformation();
  mirtkAssert(T != NULL, "point set is being transformed");
  const double t0 = wrt_pset->InputTime();
  const double t  = wrt_pset->Time();
  T->ParametricGradient(wrt_pset->InputSurfacePoints(), np_gradient, gradient, t, t0, weight);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkFloatArray> ToFloatArray(const SurfaceDistance::GradientType *v, int n)
{
  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(n);
  for (int i = 0; i < n; ++i) array->SetTuple3(i, v[i]._x, v[i]._y, v[i]._z);
  return array;
}

// -----------------------------------------------------------------------------
void SurfaceDistance::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  if (_Target->Transformation() || all) {
    snprintf(fname, sz, "%starget%s.vtp", prefix, suffix);
    WritePolyData(fname, _Target->Surface());
  }
  if (_Source->Transformation() || all) {
    snprintf(fname, sz, "%ssource%s.vtp", prefix, suffix);
    WritePolyData(fname, _Source->Surface());
  }
}

// -----------------------------------------------------------------------------
void SurfaceDistance::WriteGradient(const char               *fname,
                                    const RegisteredPointSet *data,
                                    const GradientType       *g,
                                    const Array<int>         *sample) const
{
  bool samples_only = sample && !sample->empty();
  const vtkIdType n = (samples_only ? sample->size() : data->NumberOfSurfacePoints());

  vtkSmartPointer<vtkPoints>    points;
  vtkSmartPointer<vtkCellArray> vertices;
  vtkSmartPointer<vtkPolyData>  output;

  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(n);

  vertices = vtkSmartPointer<vtkCellArray>::New();
  vertices->Allocate(n);

  double p[3];
  for (vtkIdType i = 0; i < n; ++i) {
    data->GetInputSurfacePoint((samples_only ? (*sample)[i] : i), p);
    points  ->SetPoint(i, p);
    vertices->InsertNextCell(1, &i);
  }

  vtkSmartPointer<vtkDataArray> gradient = ToFloatArray(g, n);
  gradient->SetName("gradient");

  output = vtkSmartPointer<vtkPolyData>::New();
  output->SetPoints(points);
  output->SetVerts(vertices);
  output->GetPointData()->AddArray(gradient);

  WritePolyData(fname, output);
}


} // namespace mirtk
