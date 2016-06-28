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

#define _USE_MATH_DEFINES
#include <cmath>

#include "mirtk/Math.h"


namespace mirtk {


const double inf         = numeric_limits<double>::infinity();
const double nan         = numeric_limits<double>::quiet_NaN();
const double NaN         = numeric_limits<double>::quiet_NaN();
const double pi          = double(M_PI);
const double pi_half     = 0.5 * pi;
const double two_pi      = 2.0 * pi;
const double rad_per_deg = pi / 180.0;
const double deg_per_rad = 180.0 / pi;


} // namespace mirtk
