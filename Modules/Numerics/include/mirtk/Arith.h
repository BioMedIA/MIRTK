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

#ifndef MIRTK_Arith_H
#define MIRTK_Arith_H

namespace mirtk {


int IAbs(int iValue);
int ICeil(float fValue);
int IFloor(float fValue);
int ISign(int iValue);

double Abs(double fValue);
double ACos(double fValue);
double ASin(double fValue);
double ATan(double fValue);
double ATan2(double fY, double fX);
double Ceil(double fValue);
double Cos(double fValue);
double Exp(double fValue);
double Floor(double fValue);
double Log(double fValue);
double Pow(double kBase, double kExponent);
double Sign(double fValue);
double Sin(double fValue);
double Sqr(double fValue);
double Sqrt(double fValue);
double UnitRandom();  // in [0,1]
double SymmetricRandom();  // in [-1,1]


} // namespace mirtk

#endif // MIRTK_Arith_H
