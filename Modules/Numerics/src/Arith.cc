/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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


#include "mirtk/Math.h"
#include "mirtk/Arith.h"

using namespace std;


namespace mirtk {


//----------------------------------------------------------------------------
int IAbs(int iValue)
{
  return ( iValue >= 0 ? iValue : -iValue );
}

//----------------------------------------------------------------------------
int ICeil(float fValue)
{
  return int(ceil(fValue));
}

//----------------------------------------------------------------------------
int IFloor(float fValue)
{
  return int(floor(fValue));
}

//----------------------------------------------------------------------------
int ISign(int iValue)
{
  return ( iValue > 0 ? +1 : ( iValue < 0 ? -1 : 0 ) );
}

//----------------------------------------------------------------------------
double Abs(double fValue)
{
  return double(fabs(fValue));
}

//----------------------------------------------------------------------------
double ACos(double fValue)
{
  if ( -1.0 < fValue ) {
    if ( fValue < 1.0 )
      return double(acos(fValue));
    else
      return 0.0;
  } else {
    return pi;
  }
}

//----------------------------------------------------------------------------
double ASin(double fValue)
{
  if ( -1.0 < fValue ) {
    if ( fValue < 1.0 )
      return double(asin(fValue));
    else
      return -pi_half;
  } else {
    return pi_half;
  }
}

//----------------------------------------------------------------------------
double ATan (double fValue)
{
  return double(atan(fValue));
}

//----------------------------------------------------------------------------
double ATan2(double fY, double fX)
{
  return double(atan2(fY,fX));
}

//----------------------------------------------------------------------------
double Ceil(double fValue)
{
  return double(ceil(fValue));
}

//----------------------------------------------------------------------------
double Cos(double fValue)
{
  return double(cos(fValue));
}

//----------------------------------------------------------------------------
double Exp(double fValue)
{
  return double(exp(fValue));
}

//----------------------------------------------------------------------------
double Floor(double fValue)
{
  return double(floor(fValue));
}

//----------------------------------------------------------------------------
double Log(double fValue)
{
  return double(log(fValue));
}

//----------------------------------------------------------------------------
double Pow(double fBase, double fExponent)
{
  return double(pow(fBase, fExponent));
}

//----------------------------------------------------------------------------
double Sign(double fValue)
{
  if ( fValue > 0.0 )
    return 1.0;

  if ( fValue < 0.0 )
    return -1.0;

  return 0.0;
}

//----------------------------------------------------------------------------
double Sin(double fValue)
{
  return double(sin(fValue));
}

//----------------------------------------------------------------------------
double Sqr(double fValue)
{
  return fValue*fValue;
}

//----------------------------------------------------------------------------
double Sqrt(double fValue)
{
  return double(sqrt(fValue));
}

//----------------------------------------------------------------------------
double UnitRandom()
{
  return double(rand())/double(RAND_MAX);
}

//----------------------------------------------------------------------------
double SymmetricRandom()
{
  return 2.0*double(rand())/double(RAND_MAX) - 1.0;
}


} // namespace mirtk
