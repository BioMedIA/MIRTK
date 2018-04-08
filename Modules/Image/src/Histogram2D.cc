/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
 * Copyright 2015-2017 Andreas Schuh
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

#include "mirtk/Histogram1D.h"
#include "mirtk/Histogram2D.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/GenericImage.h"

#include "mirtk/CommonExport.h"


namespace mirtk {


// Global "debug" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace Histogram2DUtils {


// -----------------------------------------------------------------------------
template <class T>
class LogTransform
{
  T     *_Bins;
  double _Num;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    double p;
    T *ptr = _Bins + re.begin();
    for (int i = re.begin(); i != re.end(); ++i, ++ptr) {
      p = static_cast<double>(*ptr);
      if (p > .0) p = log(p / _Num);
      else        p = .0;
      (*ptr) = static_cast<T>(p);
    }
  }

  static void Run(Histogram2D<T> *hxy)
  {
    LogTransform body;
    body._Bins = hxy->RawPointer();
    body._Num  = hxy->NumberOfSamples();
    blocked_range<int> i(0, hxy->NumberOfBins());
    parallel_for(i, body);
  }
};


} // namespace Histogram2DUtils
using namespace Histogram2DUtils;

// =============================================================================
// Histogram2D
// =============================================================================

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType>::Histogram2D(const Histogram2D &h)
:
  Object(h)
{
  _min_x   = h._min_x;
  _min_y   = h._min_y;
  _max_x   = h._max_x;
  _max_y   = h._max_y;
  _width_x = h._width_x;
  _width_y = h._width_y;
  _nbins_x = h._nbins_x;
  _nbins_y = h._nbins_y;
  _nsamp   = h._nsamp;
  const int nbins = NumberOfBins();
  if (nbins > 0) {
    Allocate(_bins, _nbins_x, _nbins_y);
    memcpy(RawPointer(), h.RawPointer(), nbins * sizeof(HistogramType));
  } else {
    _bins = nullptr;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType>::Histogram2D(int nbins_x, int nbins_y)
{
  if (nbins_x < 0) nbins_x = 0;
  if (nbins_y < 0) nbins_y = 0;
  _min_x   = 0;
  _min_y   = 0;
  _max_x   = nbins_x;
  _max_y   = nbins_y;
  _width_x = 1;
  _width_y = 1;
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _nsamp   = 0;
  if (nbins_x > 0 && nbins_y > 0) {
    CAllocate(_bins, _nbins_x, _nbins_y);
  } else {
    _bins = nullptr;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType>::Histogram2D(double min_x, double max_x, double width_x,
                                        double min_y, double max_y, double width_y)
{
  _min_x   = min_x;
  _min_y   = min_y;
  _max_x   = max_x;
  _max_y   = max_y;
  _nbins_x = iround((max_x - min_x) / width_x);
  _nbins_y = iround((max_y - min_y) / width_y);
  _width_x = (_max_x - _min_x) / double(_nbins_x);
  _width_y = (_max_y - _min_y) / double(_nbins_y);
  _nsamp = 0;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "Histogram2D<HistogramType>::Histogram2D: Should have at least one bin" << endl;
    exit(1);
  }
  CAllocate(_bins, _nbins_x, _nbins_y);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType>::~Histogram2D()
{
  Deallocate(_bins);
  _nbins_x = 0;
  _nbins_y = 0;
  _min_x   = 0;
  _min_y   = 0;
  _max_x   = 0;
  _max_y   = 0;
  _width_x = 0;
  _width_y = 0;
  _nsamp   = 0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Initialize(double min_x, double max_x, double width_x,
                                            double min_y, double max_y, double width_y)
{
  const int nbins_x = iround((max_x - min_x) / width_x);
  const int nbins_y = iround((max_y - min_y) / width_y);
  if (nbins_x < 1 || nbins_y < 1) {
    cerr << "Histogram2D<HistogramType>::Histogram2D: Should have at least one bin" << endl;
    exit(1);
  }
  if (_nbins_x != nbins_x || nbins_y != _nbins_y) {
    Deallocate(_bins);
    Allocate(_bins, nbins_x, nbins_y);
  }
  _min_x   = min_x;
  _min_y   = min_y;
  _max_x   = max_x;
  _max_y   = max_y;
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _width_x = (_max_x - _min_x) / double(_nbins_x);
  _width_y = (_max_y - _min_y) / double(_nbins_y);
  Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Reset()
{
  memset(RawPointer(), 0, NumberOfBins() * sizeof(HistogramType));
  _nsamp = 0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Reset(const Histogram2D &h)
{
  if (this != &h) {
    if ((_nbins_x != h._nbins_x) || (_nbins_y != h._nbins_y)) {
      Deallocate(_bins);
      Allocate(_bins, h._nbins_x, h._nbins_y);
    }
    _min_x   = h._min_x;
    _min_y   = h._min_y;
    _max_x   = h._max_x;
    _max_y   = h._max_y;
    _width_x = h._width_x;
    _width_y = h._width_y;
    _nbins_x = h._nbins_x;
    _nbins_y = h._nbins_y;
    _nsamp   = h._nsamp;
    memcpy(RawPointer(), h.RawPointer(), NumberOfBins() * sizeof(HistogramType));
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType> &Histogram2D<HistogramType>::Transpose()
{
  swap(_min_x, _min_y);
  swap(_max_x, _max_y);
  swap(_width_x, _width_y);
  swap(_nbins_x, _nbins_y);

  HistogramType **bins = Allocate<HistogramType>(_nbins_x, _nbins_y);
  for (int j = 0; j < _nbins_y; ++j)
  for (int i = 0; i < _nbins_x; ++i) {
    bins[j][i] = _bins[i][j];
  }
  Deallocate(_bins);
  _bins = bins;

  return *this;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType> Histogram2D<HistogramType>::Transposed() const
{
  Histogram2D<HistogramType> h(_nbins_y, _nbins_x);

  h._min_x   = _min_y;
  h._max_x   = _max_y;
  h._width_x = _width_y;
  h._nbins_x = _nbins_y;

  h._min_y   = _min_x;
  h._max_y   = _max_x;
  h._width_y = _width_x;
  h._nbins_y = _nbins_x;

  for (int j = 0; j < _nbins_y; ++j)
  for (int i = 0; i < _nbins_x; ++i) {
    h._bins[i][j] = _bins[j][i];
  }
  h._nsamp = _nsamp;

  return h;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::PutMin(double min_x, double min_y)
{
  _min_x = min_x;
  _min_y = min_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetMin(double *min_x, double *min_y) const
{
  *min_x = _min_x;
  *min_y = _min_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetMin(double &min_x, double &min_y) const
{
  min_x = _min_x;
  min_y = _min_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::PutMax(double max_x, double max_y)
{
  _max_x = max_x;
  _max_y = max_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetMax(double *max_x, double *max_y) const
{
  *max_x = _max_x;
  *max_y = _max_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetMax(double &max_x, double &max_y) const
{
  max_x = _max_x;
  max_y = _max_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::PutWidth(double width_x, double width_y)
{
  if ((_nbins_x > 0) && (_nbins_y > 0)) {
    Deallocate(_bins);
  }
  const int nbins_x = iround((_max_x - _min_x) / width_x);
  const int nbins_y = iround((_max_y - _min_y) / width_y);
  if (_nbins_x != nbins_x || _nbins_y != nbins_y) {
    Deallocate(_bins);
    Allocate(_bins, _nbins_x, _nbins_y);
  }
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "Histogram2D<HistogramType>::PutWidth: Should have at least one bin" << endl;
    exit(1);
  }
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetWidth(double *width_x, double *width_y) const
{
  *width_x = _width_x;
  *width_y = _width_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetWidth(double &width_x, double &width_y) const
{
  width_x = _width_x;
  width_y = _width_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::PutNumberOfBins(int nbins_x, int nbins_y)
{
  if (_nbins_x != nbins_x || _nbins_y != nbins_y) {
    Deallocate(_bins);
    Allocate(_bins, _nbins_x, _nbins_y);
  }
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "Histogram2D<HistogramType>::PutWidth: Should have at least one bin" << endl;
    exit(1);
  }
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::PutNumberOfBinsX(int nbins_x)
{
  this->PutNumberOfBins(nbins_x, this->NumberOfBinsY());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::PutNumberOfBinsY(int nbins_y)
{
  this->PutNumberOfBins(this->NumberOfBinsX(), nbins_y);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetNumberOfBins(int *nbins_x, int *nbins_y) const
{
  *nbins_x = _nbins_x;
  *nbins_y = _nbins_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::GetNumberOfBins(int &nbins_x, int &nbins_y) const
{
  nbins_x = _nbins_x;
  nbins_y = _nbins_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::AddSample(double x, double y, HistogramType n)
{
  if (x < _min_x || x > _max_x || y < _min_y || y > _max_y) return;
  int i = iround(_nbins_x * (x - _min_x - 0.5*_width_x) / (_max_x - _min_x));
  int j = iround(_nbins_y * (y - _min_y - 0.5*_width_y) / (_max_y - _min_y));
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i >= _nbins_x) i = _nbins_x - 1;
  if (j >= _nbins_y) j = _nbins_y - 1;
  _bins[j][i] += n;
  _nsamp      += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::DelSample(double x, double y, HistogramType n)
{
  if (x < _min_x || x > _max_x || y < _min_y || y > _max_y) return;
  int i = iround(_nbins_x * (x - _min_x - 0.5*_width_x) / (_max_x - _min_x));
  int j = iround(_nbins_y * (y - _min_y - 0.5*_width_y) / (_max_y - _min_y));
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i >= _nbins_x) i = _nbins_x - 1;
  if (j >= _nbins_y) j = _nbins_y - 1;
  _bins[j][i] -= n;
  _nsamp      -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::HistogramX(Histogram1D<HistogramType> &hx) const
{
  hx.PutNumberOfBins(_nbins_x);
  hx.Min(_min_x);
  hx.Max(_max_x);
  hx.NumberOfSamples(_nsamp);

  const HistogramType *bin = this->RawPointer();
  for (int j = 0; j < _nbins_y; ++j) {
    HistogramType *out = hx.RawPointer();
    for (int i = 0; i < _nbins_x; ++i, ++bin, ++out) *out += *bin;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::HistogramY(Histogram1D<HistogramType> &hy) const
{
  hy.PutNumberOfBins(_nbins_y);
  hy.Min(_min_y);
  hy.Max(_max_y);
  hy.NumberOfSamples(_nsamp);

  const HistogramType *bin = this->RawPointer();
  HistogramType       *out = hy.RawPointer();
  for (int j = 0; j < _nbins_y; ++j, ++out)
  for (int i = 0; i < _nbins_x; ++i, ++bin) {
    *out += *bin;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Log()
{
  LogTransform<HistogramType>::Run(this);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::MeanX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::MeanX: No samples in Histogram" << endl;
    }
    return 0;
  }
  double tmp, val = .0;
  for (int i = 0; i < _nbins_x; i++) {
    tmp = this->BinToValX(i);
    for (int j = 0; j < _nbins_y; j++) {
      val += _bins[j][i] * tmp;
    }
  }
  return val / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::MeanY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram_2D::MeanY: No samples in Histogram" << endl;
    }
    return 0;
  }
  double tmp, val = .0;
  for (int j = 0; j < _nbins_y; j++) {
    tmp = this->BinToValY(j);
    for (int i = 0; i < _nbins_x; i++) {
      val += _bins[j][i] * tmp;
    }
  }
  return val / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::VarianceX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::VarianceX: No samples in Histogram" << endl;
    }
    return 0;
  }
  double val = .0;
  for (int i = 0; i < _nbins_x; i++) {
    val += this->MarginalProbabilityX(i) * pow(this->BinToValX(i), 2.0);
  }
  return val - pow(this->MeanX(), 2.0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::VarianceY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::VarianceY: No samples in Histogram" << endl;
    }
    return 0;
  }
  double val = .0;
  for (int i = 0; i < _nbins_y; i++) {
    val += this->MarginalProbabilityY(i) * pow(this->BinToValY(i), 2.0);
  }
  return val - pow(this->MeanY(), 2.0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::StandardDeviationX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::StandardDeviationX: No samples in Histogram" << endl;
    }
    return 0;
  }
  return sqrt(this->VarianceX());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::StandardDeviationY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::StandardDeviationY: No samples in Histogram" << endl;
    }
    return 0;
  }
  return sqrt(this->VarianceY());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::Covariance() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::Covariance: No samples in Histogram" << endl;
    }
    return 0;
  }
  double       val    = .0;
  const double mean_x = this->MeanX();
  const double mean_y = this->MeanY();
  for (int j = 0; j < _nbins_y; ++j)
  for (int i = 0; i < _nbins_x; ++i) {
    val += _bins[j][i] * (this->BinToValX(i) - mean_x) * (this->BinToValY(j) - mean_y);
  }
  return val / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::EntropyX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::EntropyX: No samples in Histogram" << endl;
    }
    return 0;
  }
  Histogram1D<HistogramType> hx(0);
  HistogramX(hx);
  return hx.Entropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::EntropyY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::EntropyY: No samples in Histogram" << endl;
    }
    return 0;
  }
  Histogram1D<HistogramType> hy(0);
  HistogramY(hy);
  return hy.Entropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::JointEntropy() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::JointEntropy: No samples in Histogram" << endl;
    }
    return 0;
  }
  // Attention: Parallel summation yielded slightly different results each
  //            time this function was executed. This might be caused by a
  //            different summation of values, which causes different numerical
  //            cancelations. When used for NMI gradient computation, the
  //            registration result could differ from run to run!
  //
  // Additional note from 20 Sep 2017 (Andreas):
  // Changed summation from naive to Kahan summation. This results in identical
  // joint entropy and NMI values for a symmetric SVFFD registration at the first
  // iteration, lowest level, when input images are exchanged. With the naive
  // summation, the joint entropy and therefore NMI values would differ after
  // about 15 significant digits.
  double p, s = 0., c = 0., y, t;
  const HistogramType *bin = _bins[0];
  const int nbins = _nbins_x * _nbins_y;
  for (int i = 0; i != nbins; ++i, ++bin) {
    p = static_cast<double>(*bin);
    if (p > .0) {
      y = p * log(p) - c;
      t = s + y;
      c = (t - s) - y;
      s = t;
    }
  }
  // H = - sum (p/n) log(p/n) = log(n) - sum(p log p) / n
  return log(_nsamp) - s / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::ConditionalMeanXY(int i) const
{
  double m = .0, p = .0;
  for (int j = 0; j < _nbins_x; j++) {
    m += this->JointProbability(j, i) * this->BinToValX(j);
    p += this->JointProbability(j, i);
  }
  return ((p > .0) ? (m / p) : .0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::ConditionalMeanYX(int i) const
{
  double m = .0, p = .0;
  for (int j = 0; j < _nbins_y; j++) {
    m += this->JointProbability(i, j) * this->BinToValY(j);
    p += this->JointProbability(i, j);
  }
  return ((p > .0) ? (m / p) : .0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::ConditionalVarianceXY(int i) const
{
  double s = .0, p = .0;
  const double m = this->ConditionalMeanXY(i);
  for (int j = 0; j < _nbins_x; j++) {
    s += this->JointProbability(j, i) * pow(this->BinToValX(j) - m, 2);
    p += this->JointProbability(j, i);
  }
  return ((p > .0) ? (s / p) : .0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::ConditionalVarianceYX(int i) const
{
  double s = .0, p = .0;
  const double m = this->ConditionalMeanYX(i);
  for (int j = 0; j < _nbins_y; j++) {
    s += this->JointProbability(i, j) * pow(this->BinToValY(j) - m, 2);
    p += this->JointProbability(i, j);
  }
  return ((p > .0) ? (s / p) : .0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::CorrelationRatioXY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::CorrelationRatioXY: No samples in Histogram" << endl;
    }
    return 0;
  }
  double       c = 0;
  const double m = this->MeanX();
  for (int i = 0; i < _nbins_y; i++) {
    c += this->MarginalProbabilityY(i) * pow(this->ConditionalMeanXY(i) - m, 2.0);
  }
  return (c / this->VarianceX());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::CorrelationRatioYX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::CorrelationRatioYX: No samples in Histogram" << endl;
    }
    return 0;
  }
  double       c = 0;
  const double m = this->MeanY();
  for (int i = 0; i < _nbins_x; i++) {
    c += this->MarginalProbabilityX(i) * pow(this->ConditionalMeanYX(i) - m, 2.0);
  }
  return (c / this->VarianceY());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::MutualInformation() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::MutualInformation: No samples in Histogram" << endl;
    }
    return 0;
  }
  return this->EntropyX() + this->EntropyY() - this->JointEntropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::NormalizedMutualInformation() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::NormalizedMutualInformation: No samples in Histogram" << endl;
    }
    return 0;
  }
  return (this->EntropyX() + this->EntropyY()) / this->JointEntropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::CrossCorrelation() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::CrossCorrelation: No samples in Histogram" << endl;
    }
    return 0;
  }
  return abs(this->Covariance() / (sqrt(this->VarianceX()) *
                                   sqrt(this->VarianceY())));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::SumsOfSquaredDifferences() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::SumsOfSquaredDifferences:";
      cerr << " No samples in Histogram" << endl;
    }
    return 0;
  }
  double ssd = .0, val_x, val_y = this->BinToValY(0);
  for (int j = 0; j < _nbins_y; ++j) {
    val_x = this->BinToValX(0);
    for (int i = 0; i < _nbins_x; ++i) {
      ssd   += _bins[j][i] * (val_x - val_y) * (val_x - val_y);
      val_x += _width_x;
    }
    val_y += _width_y;
  }
  return ssd;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::LabelConsistency() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::LabelConsistency: No samples in Histogram" << endl;
    }
    return 0;
  }
  if (_nbins_x != _nbins_y) {
    cerr << "Histogram2D<HistogramType>::LabelConsistency: Histogram must have equal number of bins in X and Y" << endl;
    return 0;
  }
  HistogramType n = 0;
  for (int i = 0; i < _nbins_x; i++) {
    n += _bins[i][i];
  }
  return static_cast<double>(n) / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double Histogram2D<HistogramType>::Kappa() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram2D<HistogramType>::Kappa: No samples in Histogram" << endl;
    }
    return 0;
  }
  if (_nbins_x != _nbins_y) {
    cerr << "Histogram2D<HistogramType>::Kappa: Histogram must have equal number of bins in X and Y" << endl;
    return 0;
  }

  HistogramType *col_sum = new HistogramType[_nbins_x];
  HistogramType *row_sum = new HistogramType[_nbins_x];
  for (int j = 0; j < _nbins_x; j++) {
    col_sum[j] = 0;
    row_sum[j] = 0;
    for (int i = 0; i < _nbins_x; i++) {
      col_sum[j] += _bins[i][j];
      row_sum[j] += _bins[j][i];
    }
  }

  double po = .0, pe = .0;
  for (int j = 0; j < _nbins_x; j++) {
    po += _bins[j][j];
    pe += static_cast<double>(col_sum[j]) * static_cast<double>(row_sum[j]);
  }
  po /= _nsamp;
  pe /= static_cast<double>(_nsamp) * static_cast<double>(_nsamp);

  delete[] row_sum;
  delete[] col_sum;

  return (po - pe) / (1.0 - pe);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Smooth()
{
  if (_nsamp == 0) return;

  // Smoothing kernel
  double kernel[3] = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };

  // Allocate temporary memory
  double value, **tmp = Allocate<double>(_nbins_x, _nbins_y);

  // Smooth along the x-axis
  for (int j = 0; j < _nbins_y; j++)
  for (int i = 0; i < _nbins_x; i++) {
    value = .0;
    for (int k = 0; k < 3; k++) {
      if ((i-1+k >= 0) && (i-1+k < _nbins_x)) {
        value += kernel[k] * static_cast<double>(_bins[j][i-1+k]);
      }
    }
    tmp[j][i] = value;
  }

  // Smooth along the y-axis
  double n = 0., c = 0., y, t;
  for (int i = 0; i < _nbins_x; i++)
  for (int j = 0; j < _nbins_y; j++) {
    value = .0;
    for (int k = 0; k < 3; k++) {
      if ((j-1+k >= 0) && (j-1+k < _nbins_y)) {
        value += kernel[k] * tmp[j-1+k][i];
      }
    }
    _bins[j][i] = static_cast<HistogramType>(value);
    // Kahan summation of number of samples
    y = static_cast<double>(value) - c;
    t = n + y;
    c = (t - n) - y;
    n = t;
  }
  _nsamp = static_cast<HistogramType>(n);

  // Free tmp memory
  Deallocate(tmp);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
Histogram2D<HistogramType> Histogram2D<HistogramType>::Smoothed(bool pad)
{
  Histogram2D<HistogramType> hist;
  if (pad) {
    const int m = 2; // smoothing kernel size - 1
    double mx = m * _width_x;
    double my = m * _width_y;
    hist.Initialize(_min_x - mx, _max_x + mx, _width_x,
                    _min_y - my, _max_y + my, _width_y);
    for (int j = 0; j < _nbins_y; ++j)
    for (int i = 0; i < _nbins_x; ++i) {
      hist._bins[j + m][i + m] = _bins[j][i];
    }
    hist._nsamp = _nsamp;
  } else {
    hist.Reset(*this);
  }
  hist.Smooth();
  return hist;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Read(const char *filename)
{
  char buffer[255];

  ifstream from(filename);
  if (!from) {
    cerr << "Histogram2D<HistogramType>::Read: Can't open file " << filename << endl;
    exit(1);
  }
  if ((_nbins_x > 0) && (_nbins_y > 0)) {
    Deallocate(_bins);
    _nbins_x = 0;
    _nbins_y = 0;
  }

  from >> buffer;
  if (strcmp(buffer, "irtkHistogram2D") != 0) {
    cerr << "Histogram2D<HistogramType>::Read: Invalid format" << endl;
    exit(1);
  }

  // Read no. of bins
  from >> _nbins_x;
  from >> _nbins_y;

  // Read no. of samples
  from >> _nsamp;

  // Read min and max of bins
  from >> _min_x;
  from >> _max_x;
  from >> _min_y;
  from >> _max_y;

  // Read width of bins
  from >> _width_x;
  from >> _width_y;

  Allocate(_bins, _nbins_x, _nbins_y);
  for (int j = 0; j < _nbins_y; j++) {
    for (int i = 0; i < _nbins_x; i++) {
      from >> _bins[j][i];
    }
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Write(const char *filename) const
{
  ofstream to(filename);
  if (!to) {
    cerr << "Histogram2D<HistogramType>::Write: Can't open file " << filename << endl;
    exit(1);
  }
  to << "irtkHistogram2D\n";
  to << _nbins_x << " "
     << _nbins_y << " "
     << _nsamp << " "
     << _min_x << " "
     << _max_x << " "
     << _min_y << " "
     << _max_y << " "
     << _width_x << " "
     << _width_y << "\n";
  for (int j = 0; j < _nbins_y; j++) {
    for (int i = 0; i < _nbins_x; i++) {
      to << _bins[j][i] << " ";
    }
    to << "\n";
  }
  to.close();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::WriteAsImage(const char *filename) const
{
  GenericImage<HistogramType> image(_nbins_x, _nbins_y, 1);
  for (int j = 0; j < _nbins_y; ++j)
  for (int i = 0; i < _nbins_x; ++i) {
    image(i, j, 0) = _bins[j][i];
  }
  image.Write(filename);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void Histogram2D<HistogramType>::Print() const
{
  cout << _nbins_x << " "
       << _nbins_y << " "
       << _nsamp   << " "
       << _min_x   << " "
       << _max_x   << " "
       << _min_y   << " "
       << _max_y   << " "
       << _width_x << " "
       << _width_y
       << "\n";
  for (int j = 0; j < _nbins_y; j++) {
    for (int i = 0; i < _nbins_x; i++) {
      cout << _bins[j][i] << " ";
    }
    cout << "\n";
  }
}

// =============================================================================
// Explicit instantiations
// =============================================================================

template class Histogram2D<int>;
template class Histogram2D<double>;


} // namespace mirtk
