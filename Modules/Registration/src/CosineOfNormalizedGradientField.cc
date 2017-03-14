/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Stefan Pszczolkowski Parraguez, Andreas Schuh
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

#include "mirtk/CosineOfNormalizedGradientField.h"

#include "mirtk/Math.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(CosineOfNormalizedGradientField);


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace CosineOfNormalizedGradientFieldUtils {


// -----------------------------------------------------------------------------
struct EvaluateCosineOfNormalizedGradientFieldSimilarity : public VoxelReduction
{
private:

  const CosineOfNormalizedGradientField *_Similarity;

  int _dx; ///< Offset of 1st order derivative w.r.t x
  int _dy; ///< Offset of 1st order derivative w.r.t y
  int _dz; ///< Offset of 1st order derivative w.r.t z

  double _Sum;
  int    _Cnt;

public:

  EvaluateCosineOfNormalizedGradientFieldSimilarity(const CosineOfNormalizedGradientField *sim)
  :
    _Similarity(sim),
    _dx(sim->Target()->Offset(RegisteredImage::Dx)),
    _dy(sim->Target()->Offset(RegisteredImage::Dy)),
    _dz(sim->Target()->Offset(RegisteredImage::Dz)),
    _Sum(.0), _Cnt(0)
  {}

  EvaluateCosineOfNormalizedGradientFieldSimilarity(const EvaluateCosineOfNormalizedGradientFieldSimilarity &o)
  :
    _Similarity(o._Similarity),
    _dx(o._dx), _dy(o._dy), _dz(o._dz),
    _Sum(o._Sum), _Cnt(o._Cnt)
  {}

  void split(const EvaluateCosineOfNormalizedGradientFieldSimilarity &o)
  {
    _Similarity = o._Similarity;
    _dx = o._dx, _dy = o._dy, _dz = o._dz;
    _Sum = .0, _Cnt =  0;
  }

  void join(const EvaluateCosineOfNormalizedGradientFieldSimilarity &rhs)
  {
    _Sum += rhs._Sum;
    _Cnt += rhs._Cnt;
  }

  template <class TGradient>
  void operator ()(int i, int j, int k, int, const TGradient *dF, const TGradient *dM)
  {

    if (_Similarity->IsForeground(i, j, k)) {
      const int    power = _Similarity->Power();
      const double dF_dn = _Similarity->TargetNormalization();
      const double dM_dn = _Similarity->SourceNormalization();
      const double normt = sqrt(dF[_dx]*dF[_dx] + dF[_dy]*dF[_dy] + dF[_dz]*dF[_dz] + dF_dn*dF_dn);
      const double norms = sqrt(dM[_dx]*dM[_dx] + dM[_dy]*dM[_dy] + dM[_dz]*dM[_dz] + dM_dn*dM_dn);
      const double cos = (dF[_dx]*dM[_dx] + dF[_dy]*dM[_dy] + dF[_dz]*dM[_dz] + dF_dn*dM_dn) / (normt * norms);

      _Sum += pow(cos, power);
      ++_Cnt;
    }
  }
  
  double Value() const
  {
    return (_Cnt == 0 ? .0 : _Sum / _Cnt);
  }
};

// -----------------------------------------------------------------------------
struct EvaluateCosineOfNormalizedGradientFieldSimilarityGradient : public VoxelFunction
{
private:

  const CosineOfNormalizedGradientField *_Similarity;

  int _dx; ///< Offset of 1st order derivative w.r.t x
  int _dy; ///< Offset of 1st order derivative w.r.t y
  int _dz; ///< Offset of 1st order derivative w.r.t z

  int _y;  ///< Offset of output gradient w.r.t y
  int _z;  ///< Offset of output gradient w.r.t z

public:

  EvaluateCosineOfNormalizedGradientFieldSimilarityGradient(const CosineOfNormalizedGradientField *sim)
  :
    _Similarity(sim),
    _dx        (sim->Target()->Offset(RegisteredImage::Dx)),
    _dy        (sim->Target()->Offset(RegisteredImage::Dy)),
    _dz        (sim->Target()->Offset(RegisteredImage::Dz)),
    _y         (sim->NumberOfVoxels()),
    _z         (2 * _y)
  {}

  template <class TGradient, class TReal>
  void operator ()(int i, int j, int k, int, const TGradient *dF, const TGradient *dM, TReal *g)
  {
    if (_Similarity->IsForeground(i, j, k)) {
      const int    power = _Similarity->Power();
      const double dF_dn = _Similarity->TargetNormalization();
      const double dM_dn = _Similarity->SourceNormalization();
      const double normt = sqrt(dF[_dx]*dF[_dx] + dF[_dy]*dF[_dy] + dF[_dz]*dF[_dz] + dF_dn*dF_dn);
      const double norms = sqrt(dM[_dx]*dM[_dx] + dM[_dy]*dM[_dy] + dM[_dz]*dM[_dz] + dM_dn*dM_dn);
      const double cos   = (dF[_dx]*dM[_dx] + dF[_dy]*dM[_dy] + dF[_dz]*dM[_dz] + dF_dn*dM_dn) / (normt * norms);
      const double wt    = 1.0 / (normt * norms);
      const double ws    = cos / (norms * norms);

      g[0]  = wt * dF[_dx] - ws * dM[_dx];
      g[_y] = wt * dF[_dy] - ws * dM[_dy];
      g[_z] = wt * dF[_dz] - ws * dM[_dz];

      // Apply chain rule
      if (power > 1) {
        const double factor = power * pow(cos, power - 1);
        g[0] *= factor, g[_y] *= factor, g[_z] *= factor;
      }
    }
  }
};


} // namespace CosineOfNormalizedGradientFieldUtils
using namespace CosineOfNormalizedGradientFieldUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
CosineOfNormalizedGradientField
::CosineOfNormalizedGradientField(const char *name)
:
  NormalizedGradientFieldSimilarity(name),
  _Power(1)
{
}

// -----------------------------------------------------------------------------
CosineOfNormalizedGradientField
::CosineOfNormalizedGradientField(const CosineOfNormalizedGradientField &other)
:
  NormalizedGradientFieldSimilarity(other),
  _Power(1)
{
}

// -----------------------------------------------------------------------------
CosineOfNormalizedGradientField::~CosineOfNormalizedGradientField()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool CosineOfNormalizedGradientField::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Power") == 0) {
    return FromString(value, _Power) && _Power > 0;
  }
  return NormalizedGradientFieldSimilarity::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList CosineOfNormalizedGradientField::Parameter() const
{
  ParameterList params = NormalizedGradientFieldSimilarity::Parameter();
  InsertWithPrefix(params, "Power", _Power);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double CosineOfNormalizedGradientField::Evaluate()
{
  EvaluateCosineOfNormalizedGradientFieldSimilarity eval(this);
  ParallelForEachVoxel(Domain(), Target(), Source(), eval);
  if (_Power % 2 == 0) return        1.0 - eval.Value();
  else                 return 0.5 * (1.0 - eval.Value());
}

// -----------------------------------------------------------------------------
double CosineOfNormalizedGradientField::RawValue(double value) const
{
  value = NormalizedGradientFieldSimilarity::RawValue(value);
  if (_Power % 2 == 0) value = 1.0 -       value;
  else                 value = 1.0 - 2.0 * value;
  return value;
}

// -----------------------------------------------------------------------------
bool CosineOfNormalizedGradientField
::NonParametricGradient(const RegisteredImage *image, GradientImageType *gradient)
{
  typedef EvaluateCosineOfNormalizedGradientFieldSimilarityGradient EvaluateGradientFunc;

  // Select fixed image
  const RegisteredImage *fixed = ((image == Target()) ? Source() : Target());

  // Evaluate similarity gradient w.r.t gradient of given transformed image
  EvaluateGradientFunc eval(this);
  ParallelForEachVoxel(Domain(), fixed,  image, gradient, eval);
  (*gradient) *= - (_Power % 2 == 0 ? 1.0 : 0.5);

  return true;
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void CosineOfNormalizedGradientField::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  NormalizedGradientFieldSimilarity::WriteDataSets(p, suffix, all);

  const int  sz = 1024;
  char       fname[sz];
  string     _prefix = Prefix(p);
  const char *prefix = _prefix.c_str();

  ImageAttributes attr = Target()->Attributes();
  attr._t = 1;

  GenericImage<VoxelType> sim(attr);
  VoxelType *ptr2vox = sim.Data(0, 0, 0);
  const VoxelType *dtdx = Target()->Data(0, 0, 0, 1);
  const VoxelType *dtdy = Target()->Data(0, 0, 0, 2);
  const VoxelType *dtdz = Target()->Data(0, 0, 0, 3);
  const VoxelType *dsdx = Source()->Data(0, 0, 0, 1);
  const VoxelType *dsdy = Source()->Data(0, 0, 0, 2);
  const VoxelType *dsdz = Source()->Data(0, 0, 0, 3);

  sim += 1;
  for (int idx = 0; idx < NumberOfVoxels(); ++idx, ++ptr2vox, ++dtdx, ++dtdy, ++dtdz, ++dsdx, ++dsdy, ++dsdz) {
    if (IsForeground(idx)) {
      const double dtdn = _TargetNormalization;
      const double dsdn = _SourceNormalization;
      const double normt = sqrt((*dtdx)*(*dtdx) + (*dtdy)*(*dtdy) + (*dtdz)*(*dtdz) + dtdn*dtdn);
      const double norms = sqrt((*dsdx)*(*dsdx) + (*dsdy)*(*dsdy) + (*dsdz)*(*dsdz) + dsdn*dsdn);
      const double cos = ((*dtdx)*(*dsdx) + (*dtdy)*(*dsdy) + (*dtdz)*(*dsdz) + dtdn*dsdn) / (normt * norms);
      *ptr2vox = pow(cos, _Power);
    }
  }

  snprintf(fname, sz, "%sfield%s", prefix, suffix);
  sim.Write(fname);
}


} // namespace mirtk
