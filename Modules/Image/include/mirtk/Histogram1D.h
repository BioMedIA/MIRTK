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

#ifndef MIRTK_Histogram1D_H
#define MIRTK_Histogram1D_H

#include "mirtk/Object.h"
#include "mirtk/Memory.h"
#include "mirtk/Math.h"


namespace mirtk {


/**
 * Class for 1D histograms.
 */
template <class HistogramType>
class Histogram1D : public Object
{
  mirtkObjectMacro(Histogram1D);

protected:

  /// Number of bins
  int _nbins;

  /// Number of samples
  HistogramType _nsamp;

  /// Min. value for samples, everything below is ignored
  double _min;

  /// Max. value for samples, everything below is ignored
  double _max;

  /// Width of bins
  double _width;

  /// Dynamic memory for bins
  HistogramType *_bins;

public:

  /// Construct a histogram from another histogram
  Histogram1D(const Histogram1D &);

  /// Construct a histogram with 256 bins and samples ranging from 0 to 255
  Histogram1D(int nbins = 256);

  /// Construct a histogram for samples ranging from min to max and width
  Histogram1D(double min, double max, double width);

  /// Read constructor
  Histogram1D(const char *);

  /// Assignment operator
  Histogram1D &operator =(const Histogram1D &);

  /// Destructor
  ~Histogram1D();

  /// Get raw pointer to histogram bins
  HistogramType *RawPointer();

  /// Get raw pointer to histogram bins
  const HistogramType *RawPointer() const;

  /// Clear and reset histogram
  void Reset();

  /// Get number of bins in histogram
  int NumberOfBins() const;

  /// Put number of bins in histogram
  /// \note Unlike PutNumberOfBins, this function does not reset the histogram.
  void NumberOfBins(int);

  /// Put number of bins in histogram
  void PutNumberOfBins(int);

  /// Get minimum value in histogram
  double Min() const;

  /// Set minimum value in histogram
  /// \note Unlike PutMin, this function does not reset the histogram.
  void Min(double);

  /// Get maximum value in histogram
  double Max() const;

  /// Set maximum value in histogram
  /// \note Unlike PutMax, this function does not reset the histogram.
  void Max(double);

  /// Get minimum value in histogram
  double GetMin() const;

  /// Put minimum value in histogram
  void PutMin(double);

  /// Get maximum value in histogram
  double GetMax() const;

  /// Put maximum value in histogram
  void PutMax(double);

  /// Get width of bins in histogram
  double GetWidth() const;

  /// Put width of bins in histogram
  void PutWidth(double);

  /// Get number of samples in histogram
  HistogramType NumberOfSamples() const;

  /// Set number of samples in histogram
  void NumberOfSamples(HistogramType);

  /// Get number of samples in bin(i)
  HistogramType &operator()(int);

  /// Get number of samples in bin(i)
  const HistogramType &operator()(int) const;

  /// Add counts to bin
  void Add(int, HistogramType = 1);

  /// Delete counts from bin
  void Delete(int, HistogramType = 1);

  /// Add sample to bin
  void AddSample(double, HistogramType = 1);

  /// Delete sample from bin
  void DelSample(double, HistogramType = 1);

  /// Convert sample value to continuous bin index
  double ValToRange(double val) const;

  /// Convert sample value to bin index
  int ValToBin(double val) const;

  /// Convert bin index to sample value
  double BinToVal(int bin) const;

  /// Convert bin into probability density distributions
  double BinToPDF(int bin) const;

  /// Convert sample value into probability density distributions
  double ValToPDF(double val) const;

  /// Convert bin into cumulative density distributions
  double BinToCDF(int bin) const;

  /// Convert sample value into cumulative  density distributions
  double ValToCDF(double val) const;

  /// Convert cumulative density distributions to bin value
  int CDFToBin(double p) const;

  /// Convert cumulative density distributions to sample value
  double CDFToVal(double p) const;

  /// Log transform histogram
  void Log();

  /// Smooth histogram
  void Smooth();

  /// Return smallest modal value
  double Mode() const;

  /// Return sorted modal values
  Array<double> Modes() const;

  /// Calculate mean
  double Mean() const;

  /// Calculate variance
  double Variance() const;

  /// Calculate standard deviation
  double StandardDeviation() const;

  /// Calculate entropy
  double Entropy() const;

  /// Read histogram
  void Read(const char *);

  /// Wrirte histogram
  void Write(const char *) const;

  /// Print histogram
  void Print() const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType *Histogram1D<HistogramType>::RawPointer()
{
  return _bins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline const HistogramType *Histogram1D<HistogramType>::RawPointer() const
{
  return _bins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram1D<HistogramType>::NumberOfBins() const
{
  return _nbins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::NumberOfBins(int nbins)
{
  if (_nbins != nbins) {
    Deallocate(_bins);
    Allocate(_bins, nbins);
    _nbins = nbins;
    _width = (_max - _min) / _nbins;
    _nsamp = 0;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::Min() const
{
  return _min;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::Min(double min)
{
  _min   = min;
  _width = (_max - _min) / _nbins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::GetMin() const
{
  return _min;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::Max() const
{
  return _max;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::Max(double max)
{
  _max   = max;
  _width = (_max - _min) / _nbins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::GetMax() const
{
  return _max;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::GetWidth() const
{
  return _width;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType Histogram1D<HistogramType>::NumberOfSamples() const
{
  return _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::NumberOfSamples(HistogramType n)
{
  _nsamp = n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType &Histogram1D<HistogramType>::operator()(int i)
{
  return _bins[i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline const HistogramType &Histogram1D<HistogramType>::operator()(int i) const
{
  return _bins[i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::Add(int i, HistogramType n)
{
  _bins[i] += n;
  _nsamp   += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::Delete(int i, HistogramType n)
{
  _bins[i] -= n;
  _nsamp   -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::AddSample(double x, HistogramType n)
{
  if (x < _min || x > _max) return;
  int index = iround(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index <  0     ) index = 0;
  if (index >= _nbins) index = _nbins - 1;
  _bins[index] += n;
  _nsamp       += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram1D<HistogramType>::DelSample(double x, HistogramType n)
{
  if (x < _min || x > _max) return;
  int index = iround(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index <  0     ) index = 0;
  if (index >= _nbins) index = _nbins - 1;
  _bins[index] -= n;
  _nsamp       -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::ValToRange(double val) const
{
  return _nbins * (val - _min - 0.5 * _width) / (_max - _min);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram1D<HistogramType>::ValToBin(double val) const
{
  const int index = iround(ValToRange(val));
  return (index < 0 ? 0 : (index >= _nbins ? _nbins - 1 : index));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::BinToVal(int i) const
{
  return (i*(_max - _min)/_nbins + _min) + _width/2.0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::BinToPDF(int i) const
{
  return ((_nsamp == .0) ? .0 : _bins[i] / _nsamp);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::ValToPDF(double val) const
{
  return BinToPDF(ValToBin(val));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::BinToCDF(int i) const
{
  if (_nsamp == .0) return .0;
  double s = .0;
  for (int j = 0; j <= i; ++j) {
    s += _bins[j];
  }
  return s / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::ValToCDF(double val) const
{
  return BinToCDF(ValToBin(val));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram1D<HistogramType>::CDFToBin(double p) const
{
  if (p < .0 || p > 1.0) {
    cerr << "Histogram1D<HistogramType>::CDFToBin: Must be between 0 and 1" << endl;
    exit(1);
  }
  int i;
  double sum = .0;
  for (i = 0; i < _nbins; ++i) {
    sum += static_cast<double>(_bins[i]);
    if (sum / _nsamp >= p) return i;
  }
  return _nbins - 1;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram1D<HistogramType>::CDFToVal(double p) const
{
  return BinToVal(CDFToBin(p));
}


} // namespace mirtk

#endif // MIRTK_Histogram1D_H
