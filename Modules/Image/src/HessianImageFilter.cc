/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/HessianImageFilter.h"

#include "mirtk/Matrix.h"
#include "mirtk/ImageAttributes.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
HessianImageFilter<VoxelType>::HessianImageFilter(OutputType type)
{
  _Type           = type;
  _UseVoxelSize   = true;
  _UseOrientation = false;
  _PaddingValue   = -numeric_limits<double>::infinity();
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void HessianImageFilter<VoxelType>::Initialize()
{
  ImageToImage<VoxelType>::Initialize(false);

  if (this->Input()->T() > 1) {
    cerr << this->NameOfClass() << "::Run: Only implemented for images with t = 1" << endl;
    exit(1);
  }

  int n;
  if      (_Type == HESSIAN_MATRIX) n = 9;
  else if (_Type == HESSIAN_VECTOR) n = 6;
  else                              n = 1;
  this->Output()->Initialize(this->Input()->Attributes(), n);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HessianImageFilter<VoxelType>::Run()
{
  double dxx, dxy, dxz, dyy, dyz, dzz, dii, dij, dik, djj, djk, dkk;
  int    x1, y1, z1, x2, y2, z2;

  this->Initialize();

  const BaseImage       *input  = this->Input();
  BaseImage             *output = this->Output();
  const ImageAttributes &attr   = input->Attributes();
  Matrix                 R      = attr.GetWorldToImageOrientation();

  for (int z = 0; z < attr._z; ++z) {
    z1 = z - 1;
    if (z1 < 0) z1 = 0;
    z2 = z + 1;
    if (z2 > attr._z - 1) z2 = attr._z - 1;

    for (int y = 0; y < attr._y; ++y) {
      y1 = y - 1;
      if (y1 < 0) y1 = 0;
      y2 = y + 1;
      if (y2 > attr._y - 1) y2 = attr._y - 1;

      for (int x = 0; x < attr._x; ++x) {
        x1 = x - 1;
        if (x1 < 0) x1 = 0;
        x2 = x + 1;
        if (x2 > attr._x - 1) x2 = attr._x - 1;

        if (x1 != x2 &&
            input->Get(x,  y, z) > _PaddingValue &&
            input->Get(x1, y, z) > _PaddingValue &&
            input->Get(x2, y, z) > _PaddingValue) {
          dxx = (input->Get(x2, y, z) - 2.0 * input->Get(x, y, z) + input->Get(x1, y, z));
        } else {
          dxx = .0;
        }
        if (x1 != x2 &&
            y1 != y2 &&
            input->Get(x1, y1, z) > _PaddingValue &&
            input->Get(x1, y2, z) > _PaddingValue &&
            input->Get(x2, y1, z) > _PaddingValue &&
            input->Get(x2, y2, z) > _PaddingValue) {
          dxy = (input->Get(x2, y2, z) - input->Get(x2, y1, z) - input->Get(x1, y2, z) + input->Get(x1, y1, z)) / ((x2 - x1) * (y2 - y1));
        } else {
          dxy = .0;
        }
        if (x1 != x2 &&
            z1 != z2 &&
            input->Get(x1, y, z1) > _PaddingValue &&
            input->Get(x1, y, z2) > _PaddingValue &&
            input->Get(x2, y, z1) > _PaddingValue &&
            input->Get(x2, y, z2) > _PaddingValue) {
          dxz = (input->Get(x2, y, z2) - input->Get(x2, y, z1) - input->Get(x1, y, z2) + input->Get(x1, y, z1)) / ((x2 - x1) * (z2 - z1));
        } else {
          dxz = .0;
        }

        if (y1 != y2 &&
            input->Get(x, y,  z) > _PaddingValue &&
            input->Get(x, y1, z) > _PaddingValue &&
            input->Get(x, y2, z) > _PaddingValue) {
          dyy = (input->Get(x, y2, z) - 2.0 * input->Get(x, y, z) + input->Get(x, y1, z));
        } else {
          dyy = .0;
        }

        if (y1 != y2 &&
            z1 != z2 &&
            input->Get(x, y1, z1) > _PaddingValue &&
            input->Get(x, y1, z2) > _PaddingValue &&
            input->Get(x, y2, z1) > _PaddingValue &&
            input->Get(x, y2, z2) > _PaddingValue) {
          dyz = (input->Get(x, y2, z2) - input->Get(x, y2, z1) - input->Get(x, y1, z2) + input->Get(x, y1, z1)) / ((y2 - y1) * (z2 - z1));
        } else {
          dyz = .0;
        }

        if (z1 != z2 &&
            input->Get(x, y, z)  > _PaddingValue &&
            input->Get(x, y, z1) > _PaddingValue &&
            input->Get(x, y, z2) > _PaddingValue) {
          dzz = (input->Get(x, y, z2) - 2.0 * input->Get(x, y, z) + input->Get(x, y, z1));
        } else {
          dzz = .0;
        }

        if (_UseVoxelSize) {
          if (input->XSize() > 0.) {
            dxx /= (input->XSize() * input->XSize());
            dxy /= input->XSize();
            dxz /= input->XSize();
          }
          if (input->YSize() > 0.) {
            dxy /= input->YSize();
            dyy /= (input->YSize() * input->YSize());
            dyz /= input->YSize();
          }
          if (input->ZSize() > 0.) {
            dxz /= input->ZSize();
            dyz /= input->ZSize();
            dzz /= (input->ZSize() * input->ZSize());
          }
        }
        if (_UseOrientation) {
          // Using numerator-layout for matrix calculus.
          // http://en.wikipedia.org/wiki/Matrix_calculus#Numerator-layout_notation
          //
          // Expression computed here is transpose(R) * Hessian * R = transpose(Hessian * R) * R
          dii = dxx, dij = dxy, dik = dxz, djj = dyy, djk = dyz, dkk = dzz;
          dxx = R(0, 0) * (R(0, 0) * dii + R(1, 0) * dij + R(2, 0) * dik) + R(1, 0) * (R(0, 0) * dij + R(1, 0) * djj + R(2, 0) * djk) + R(2, 0) * (R(0, 0) * dik + R(1, 0) * djk + R(2, 0) * dkk);
          dxy = R(0, 1) * (R(0, 0) * dii + R(1, 0) * dij + R(2, 0) * dik) + R(1, 1) * (R(0, 0) * dij + R(1, 0) * djj + R(2, 0) * djk) + R(2, 1) * (R(0, 0) * dik + R(1, 0) * djk + R(2, 0) * dkk);
          dxz = R(0, 2) * (R(0, 0) * dii + R(1, 0) * dij + R(2, 0) * dik) + R(1, 2) * (R(0, 0) * dij + R(1, 0) * djj + R(2, 0) * djk) + R(2, 2) * (R(0, 0) * dik + R(1, 0) * djk + R(2, 0) * dkk);
          dyy = R(0, 1) * (R(0, 1) * dii + R(1, 1) * dij + R(2, 1) * dik) + R(1, 1) * (R(0, 1) * dij + R(1, 1) * djj + R(2, 1) * djk) + R(2, 1) * (R(0, 1) * dik + R(1, 1) * djk + R(2, 1) * dkk);
          dyz = R(0, 2) * (R(0, 1) * dii + R(1, 1) * dij + R(2, 1) * dik) + R(1, 2) * (R(0, 1) * dij + R(1, 1) * djj + R(2, 1) * djk) + R(2, 2) * (R(0, 1) * dik + R(1, 1) * djk + R(2, 1) * dkk);
          dzz = R(0, 2) * (R(0, 2) * dii + R(1, 2) * dij + R(2, 2) * dik) + R(1, 2) * (R(0, 2) * dij + R(1, 2) * djj + R(2, 2) * djk) + R(2, 2) * (R(0, 2) * dik + R(1, 2) * djk + R(2, 2) * dkk);
        }

        switch (_Type) {
          case HESSIAN_XX:
            output->PutAsDouble(x, y, z, 0, dxx);
            break;
          case HESSIAN_XY:
            output->PutAsDouble(x, y, z, 0, dxy);
            break;
          case HESSIAN_XZ:
            output->PutAsDouble(x, y, z, 0, dxz);
            break;
          case HESSIAN_YY:
            output->PutAsDouble(x, y, z, 0, dyy);
            break;
          case HESSIAN_YZ:
            output->PutAsDouble(x, y, z, 0, dyz);
            break;
          case HESSIAN_ZZ:
            output->PutAsDouble(x, y, z, 0, dzz);
            break;
          case HESSIAN_VECTOR:
            output->PutAsDouble(x, y, z, 0, dxx);
            output->PutAsDouble(x, y, z, 1, dxy);
            output->PutAsDouble(x, y, z, 2, dxz);
            output->PutAsDouble(x, y, z, 3, dyy);
            output->PutAsDouble(x, y, z, 4, dyz);
            output->PutAsDouble(x, y, z, 5, dzz);
            break;
          case HESSIAN_MATRIX:
            output->PutAsDouble(x, y, z, 0, dxx);
            output->PutAsDouble(x, y, z, 1, dxy);
            output->PutAsDouble(x, y, z, 2, dxz);
            output->PutAsDouble(x, y, z, 3, dxy);
            output->PutAsDouble(x, y, z, 4, dyy);
            output->PutAsDouble(x, y, z, 5, dyz);
            output->PutAsDouble(x, y, z, 6, dxz);
            output->PutAsDouble(x, y, z, 7, dyz);
            output->PutAsDouble(x, y, z, 8, dzz);
            break;
          default:
            cerr << this->NameOfClass() << "::Run: Unknown gradient computation" << endl;
            exit(1);
        }
      }
    }
  }

  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class HessianImageFilter<unsigned char>;
template class HessianImageFilter<short>;
template class HessianImageFilter<float>;
template class HessianImageFilter<double>;


} // namespace mirtk
