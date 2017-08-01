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

#ifndef MIRTK_Histogram2D_H
#define MIRTK_Histogram2D_H

#include "mirtk/Math.h"
#include "mirtk/Object.h"
#include "mirtk/Histogram1D.h"


namespace mirtk {


/**
 * Class for 2D histograms.
 */
template <class HistogramType>
class Histogram2D : public Object
{
  mirtkObjectMacro(Histogram2D);

public:

  typedef HistogramType BinType;

private:

  /// Number of bins in x-direction
  int _nbins_x;

  /// Number of bins in x-direction
  int _nbins_y;

  /// Number of samples
  HistogramType _nsamp;

  /// Min. value for x-samples, everything below is ignored
  double _min_x;

  /// Min. value for y-samples, everything below is ignored
  double _min_y;

  /// Max. value for x-samples, everything above is ignored
  double _max_x;

  /// Max. value for y-samples, everything above is ignored
  double _max_y;

  /// Width of bin in x-direction
  double _width_x;

  /// Width of bin in y-direction
  double _width_y;

  /// Dynamic memory for bins
  HistogramType **_bins;

public:

  /// Construct a histogram from another histogram
  Histogram2D(const Histogram2D &);

  /// Construct a histogram with 256 bins and samples ranging from 0 to 255
  Histogram2D(int nbins_x = 256, int nbins_y = 256);

  /// Construct a histogram for samples ranging from min to max and width
  Histogram2D(double min_x, double max_x, double width_x,
              double min_y, double max_y, double width_y);

  /// Destructor
  ~Histogram2D();

  /// Assignment operator
  Histogram2D &operator =(const Histogram2D &);

  /// Construct a histogram for samples ranging from min to max and width
  void Initialize(double min_x, double max_x, double width_x,
                  double min_y, double max_y, double width_y);

  /// Clear and reset histogram
  void Reset();

  /// Get raw pointer to histogram bins
  HistogramType *RawPointer();

  /// Get raw pointer to histogram bins
  const HistogramType *RawPointer() const;

  /// Clear and copy histogram
  void Reset(const Histogram2D &);

  /// Transpose histogram in place with x and y direction exchanged
  Histogram2D<HistogramType> &Transpose();

  /// Get transposed histogram with x and y direction exchanged
  Histogram2D<HistogramType> Transposed() const;

  /// Get total number of bins
  int  NumberOfBins() const;

  /// Get number of bins in x-direction
  int  NumberOfBinsX() const;

  /// Put number of bins in x-direction
  void PutNumberOfBinsX(int);

  /// Get number of bins in x-direction
  int  NumberOfBinsY() const;

  /// Put number of bins in x-direction
  void PutNumberOfBinsY(int);

  /// Get number of bins in x- and y-direction
  void GetNumberOfBins(int *, int *) const;

  /// Get number of bins in x- and y-direction
  void GetNumberOfBins(int &, int &) const;

  /// Put number of bins in x- and y-direction
  void PutNumberOfBins(int, int);

  /// Get minimum value in histogram in x-direction
  double MinX() const;

  /// Get maximum value in histogram in x-direction
  double MaxX() const;

  /// Get minimum value in histogram in y-direction
  double MinY() const;

  /// Get maximum value in histogram in y-direction
  double MaxY() const;

  /// Get width of bins in histogram in x-direction
  double WidthX() const;

  /// Get width of bins in histogram in y-direction
  double WidthY() const;

  /// Get minimum value in histogram
  void GetMin(double *, double *) const;

  /// Get minimum value in histogram
  void GetMin(double &, double &) const;

  /// Put minimum value in histogram
  void PutMin(double, double);

  /// Get maximum value in histogram
  void GetMax(double *, double *) const;

  /// Get maximum value in histogram
  void GetMax(double &, double &) const;

  /// Put maximum value in histogram
  void PutMax(double, double);

  /// Get width of bins in histogram
  void GetWidth(double *, double *) const;

  /// Get width of bins in histogram
  void GetWidth(double &, double &) const;

  /// Put width of bins in histogram
  void PutWidth(double, double);

  /// Get number of samples in histogram
  HistogramType NumberOfSamples() const;

  /// Set number of samples in histogram
  void NumberOfSamples(HistogramType);

  /// Get number of samples in bin(i, j)
  HistogramType &operator()(int, int);

  /// Get number of samples in bin(i, j)
  const HistogramType &operator()(int, int) const;

  /// Add counts to bins
  void Add(int, int, HistogramType = 1);

  /// Delete counts from bins
  void Delete(int, int, HistogramType = 1);

  /// Add samples
  void AddSample(double, double, HistogramType = 1);

  /// Delete samples
  void DelSample(double, double, HistogramType = 1);

  /// Convert sample value to bin index
  int  ValToBinX(double val) const;

  /// Convert bin index to sample value
  double BinToValX(int bin) const;

  /// Convert sample value to bin index
  int  ValToBinY(double val) const;

  /// Convert bin index sample value
  double BinToValY(int bin) const;

  /// Compute marginal histogram of X
  void HistogramX(Histogram1D<HistogramType> &) const;

  /// Compute marginal histogram of X
  Histogram1D<HistogramType> HistogramX() const;

  /// Compute marginal histogram of Y
  void HistogramY(Histogram1D<HistogramType> &) const;

  /// Compute marginal histogram of Y
  Histogram1D<HistogramType> HistogramY() const;

  /// Log transform histogram
  void Log();

  /// Smooth histogram
  void Smooth();

  /// Get smoothed histogram, optionally with padded boundaries
  Histogram2D<HistogramType> Smoothed(bool = false);

  /// Calculate joint probability p(x, y)
  double JointProbability(int, int) const;

  /// Calculate marginal probability p(x)
  double MarginalProbabilityX(int) const;

  /// Calculate marginal probability p(y)
  double MarginalProbabilityY(int) const;

  /// Calculate conditional probability p(x|y)
  double ConditionalProbabilityXY(int, int) const;

  /// Calculate conditional probability p(y|x)
  double ConditionalProbabilityYX(int, int) const;

  /// Calculate mean
  double MeanX() const;

  /// Calculate mean
  double MeanY() const;

  /// Calculate conditional mean
  double ConditionalMeanXY(int) const;

  /// Calculate conditional mean
  double ConditionalMeanYX(int) const;

  /// Calculate variance
  double VarianceX() const;

  /// Calculate variance
  double VarianceY() const;

  /// Calculate conditional variance
  double ConditionalVarianceXY(int) const;

  /// Calculate conditional variance
  double ConditionalVarianceYX(int) const;

  /// Calculate covariance
  double Covariance() const;

  /// Calculate standard deviation
  double StandardDeviationX() const;

  /// Calculate standard deviation
  double StandardDeviationY() const;

  /// Calculate marginal entropy
  double EntropyX() const;

  /// Calculate marginal entropy
  double EntropyY() const;

  /// Calculate joint entropy
  double JointEntropy() const;

  /// Calculate mutual information
  double MutualInformation() const;

  /// Calculate normalized mutual information
  double NormalizedMutualInformation() const;

  /// Calculate cross correlation
  double CrossCorrelation() const;

  /// Calculate correlation ratio
  double CorrelationRatioXY() const;

  /// Calculate correlation ratio
  double CorrelationRatioYX() const;

  /// Calculate sums of squared differences
  double SumsOfSquaredDifferences() const;

  /// Calcualate label consistency
  double LabelConsistency() const;

  /// Calcualate kappa statistic
  double Kappa() const;

  /// Read histogram
  void Read(const char *);

  /// Write histogram
  void Write(const char *) const;

  /// Write histogram as 2D image
  void WriteAsImage(const char *) const;

  /// Print histogram
  void Print() const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class HistogramType>
inline Histogram2D<HistogramType> &
Histogram2D<HistogramType>::operator =(const Histogram2D &h)
{
  Reset(h);
  return *this;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram2D<HistogramType>::NumberOfBins() const
{
  return _nbins_x * _nbins_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram2D<HistogramType>::NumberOfBinsX() const
{
  return _nbins_x;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram2D<HistogramType>::NumberOfBinsY() const
{
  return _nbins_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType *Histogram2D<HistogramType>::RawPointer()
{
  return _bins[0];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline const HistogramType *Histogram2D<HistogramType>::RawPointer() const
{
  return _bins[0];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType Histogram2D<HistogramType>::NumberOfSamples() const
{
  return _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram2D<HistogramType>::NumberOfSamples(HistogramType n)
{
  _nsamp = n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::MinX() const
{
  return _min_x;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::MaxX() const
{
  return _max_x;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::MinY() const
{
  return _min_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::MaxY() const
{
  return _max_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::WidthX() const
{
  return _width_x;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::WidthY() const
{
  return _width_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType &Histogram2D<HistogramType>::operator()(int i, int j)
{
  return _bins[j][i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline const HistogramType &Histogram2D<HistogramType>::operator()(int i, int j) const
{
  return _bins[j][i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram2D<HistogramType>::Add(int i, int j, HistogramType n)
{
  _bins[j][i] += n;
  _nsamp      += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void Histogram2D<HistogramType>::Delete(int i, int j, HistogramType n)
{
  _bins[j][i] -= n;
  _nsamp      -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram2D<HistogramType>::ValToBinX(double val) const
{
  int index = static_cast<int>(round(_nbins_x * (val - _min_x - 0.5*_width_x) / (_max_x - _min_x)));
  if (index < 0) index = 0;
  if (index > _nbins_x-1) index = _nbins_x - 1;
  return index;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int Histogram2D<HistogramType>::ValToBinY(double val) const
{
  int index = static_cast<int>(round(_nbins_y * (val - _min_y - 0.5*_width_y) / (_max_y - _min_y)));
  if (index < 0) index = 0;
  if (index > _nbins_y-1) index = _nbins_y - 1;
  return index;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::BinToValX(int i) const
{
  return (i*(_max_x - _min_x)/_nbins_x + _min_x) + _width_x/2.0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::BinToValY(int i) const
{
  return (i*(_max_y - _min_y)/_nbins_y + _min_y) + _width_y/2.0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline Histogram1D<HistogramType> Histogram2D<HistogramType>::HistogramX() const
{
  Histogram1D<HistogramType> hx(_nbins_x);
  HistogramX(hx);
  return hx;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline Histogram1D<HistogramType> Histogram2D<HistogramType>::HistogramY() const
{
  Histogram1D<HistogramType> hy(_nbins_y);
  HistogramY(hy);
  return hy;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::JointProbability(int i, int j) const
{
  return _bins[j][i] / (double) _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::MarginalProbabilityX(int i) const
{
  HistogramType n = 0;
  for (int j = 0; j < _nbins_y; j++) {
    n += _bins[j][i];
  }
  return n / (double) _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::MarginalProbabilityY(int i) const
{
  HistogramType n = 0;
  for (int j = 0; j < _nbins_x; j++) {
    n += _bins[i][j];
  }
  return static_cast<double>(n) / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::ConditionalProbabilityXY(int i, int j) const
{
  double p = this->MarginalProbabilityY(j);
  if (p > 0) {
    return this->JointProbability(i, j) / p;
  } else {
    return 0;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double Histogram2D<HistogramType>::ConditionalProbabilityYX(int i, int j) const
{
  double p = this->MarginalProbabilityX(j);
  if (p > 0) {
    return this->JointProbability(j, i) / p;
  } else {
    return 0;
  }
}


} // namespace mirtk

#endif // MIRTK_Histogram2D_H
