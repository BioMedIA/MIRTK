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

#ifndef MIRTK_DataFunctions_H
#define MIRTK_DataFunctions_H

#include "mirtk/DataOp.h"
#include "mirtk/DataStatistics.h"
#include "mirtk/Parallel.h"


namespace mirtk { namespace data { namespace op {


// -----------------------------------------------------------------------------
/// Reset data mask
class ResetMask : public Op
{
private:

  bool _Inverse;

public:

  ResetMask(bool inv = false) : _Inverse(inv) {}

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *, bool *mask = nullptr)
  {
    if (mask != nullptr) {
      memset(mask, _Inverse ? 1 : 0, n * sizeof(bool));
    }
  }
};

// -----------------------------------------------------------------------------
/// Invert mask values
class InvertMask : public Op
{
public:

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *, bool *mask = nullptr)
  {
    if (mask != nullptr) {
      for (int i = 0; i < n; ++i) {
        mask[i] = !mask[i];
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Set value of unmasked data points
class SetInsideValue : public Op
{
private:

  double *_Data;
  bool   *_Mask;
  double  _Value;

public:

  SetInsideValue(double value = .0) : _Value(value) {}

  // Called by TBB's parallel_for, for internal use only!
  void operator()(const blocked_range<int> &re) const
  {
    double *data = _Data + re.begin();
    if (_Mask) {
      bool *mask = _Mask + re.begin();
      for (int i = re.begin(); i != re.end(); ++i, ++data, ++mask) {
        if (*mask == true) *data = _Value;
      }
    } else {
      for (int i = re.begin(); i != re.end(); ++i, ++data) {
        *data = _Value;
      }
    }
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    _Data = data;
    _Mask = mask;
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Set value of masked data points
class SetOutsideValue : public Op
{
private:

  double *_Data;
  bool   *_Mask;
  double  _Value;

public:

  SetOutsideValue(double value = .0) : _Value(value) {}

  // Called by TBB's parallel_for, for internal use only!
  void operator()(const blocked_range<int> &re) const
  {
    double *data = _Data + re.begin();
    bool   *mask = _Mask + re.begin();
    for (int i = re.begin(); i != re.end(); ++i, ++data, ++mask) {
      if (*mask == false) *data = _Value;
    }
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    if (mask) {
      _Data = data;
      _Mask = mask;
      parallel_for(blocked_range<int>(0, n), *this);
    }
  }
};

// -----------------------------------------------------------------------------
/// Base class of element-wise data transformations
class ElementWiseUnaryOp : public Op
{
private:

  double *_Data;
  bool   *_Mask;

public:

  // Called by TBB's parallel_for, for internal use only!
  void operator()(const blocked_range<int> &re) const
  {
    double *data = _Data + re.begin();
    if (_Mask) {
      bool *mask = _Mask + re.begin();
      for (int i = re.begin(); i != re.end(); ++i, ++data, ++mask) {
        if (*mask) *data = this->Op(*data, *mask);
      }
    } else {
      bool mask = true;
      for (int i = re.begin(); i != re.end(); ++i, ++data) {
        *data = this->Op(*data, mask);
      }
    }
  }

  /// Transform data value and/or mask it by setting mask = false
  virtual double Op(double value, bool &) const = 0;

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    _Data = data;
    _Mask = mask;
    // MUST be called in the base class which defines Op!
    //parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Base class of element-wise data transformations
class ElementWiseBinaryOp : public Op
{
  /// Constant value to add
  mirtkPublicAttributeMacro(double, Constant);

  /// File path of second dataset
  mirtkPublicAttributeMacro(string, FileName);

private:

  double *_Data;
  bool   *_Mask;
  double *_Other;

protected:

  /// Constructor
  ElementWiseBinaryOp(double value) : _Constant(value) {}

  /// Constructor
  ElementWiseBinaryOp(const char *fname) : _Constant(.0), _FileName(fname) {}

public:

  // Called by TBB's parallel_for, for internal use only!
  void operator()(const blocked_range<int> &re) const
  {
    double *data = _Data + re.begin();
    if (_Other) {
      double *other = _Other + re.begin();
      if (_Mask) {
        bool *mask = _Mask + re.begin();
        for (int i = re.begin(); i != re.end(); ++i) {
          if (*mask) *data = this->Op(*data, *other, *mask);
          ++data, ++other, ++mask;
        }
      } else {
        bool mask = true;
        for (int i = re.begin(); i != re.end(); ++i) {
          *data = this->Op(*data, *other, mask);
          ++data, ++other;
        }
      }
    } else {
      if (_Mask) {
        bool *mask = _Mask + re.begin();
        for (int i = re.begin(); i != re.end(); ++i) {
          if (*mask) *data = this->Op(*data, _Constant, *mask);
          ++data, ++mask;
        }
      } else {
        bool mask = true;
        for (int i = re.begin(); i != re.end(); ++i) {
          *data = this->Op(*data, _Constant, mask);
          ++data;
        }
      }
    }
  }

  /// Transform data value and/or mask it by setting mask = false
  virtual double Op(double value, double, bool &) const = 0;

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    _Data     = data;
    _Mask     = mask;
    _Other    = NULL;
    if (!_FileName.empty()) {
      if (Read(_FileName.c_str(), _Other) != n) {
        cerr << "Input file " << _FileName << " has different number of data points!" << endl;
        exit(1);
      }
    }
    // MUST be called in the base class which defines Op!
    //parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise absolute value
class Abs : public ElementWiseUnaryOp
{
public:

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &) const
  {
    return abs(value);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise power
class Pow : public ElementWiseUnaryOp
{
  mirtkPublicAttributeMacro(double, Exponent);

public:

  Pow(double exponent) : _Exponent(exponent) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &) const
  {
    return pow(value, _Exponent);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise exponential map
class Exp : public ElementWiseUnaryOp
{
public:

  /// Transform data value and/or mask it by setting mask = false
  virtual double Op(double value, bool &) const
  {
    return exp(value);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise logarithmic map
class Log : public ElementWiseUnaryOp
{
  mirtkPublicAttributeMacro(double, Base);
  mirtkPublicAttributeMacro(double, Threshold);

public:

  Log(double base, double threshold = .01) : _Base(base), _Threshold(threshold) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &) const
  {
    if (value < _Threshold) value = _Threshold;
    return log(value) / log(_Base);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise binary logarithm
class Lb : public ElementWiseUnaryOp
{
  mirtkPublicAttributeMacro(double, Threshold);

public:

  Lb(double threshold = .01) : _Threshold(threshold) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &) const
  {
    if (value < _Threshold) value = _Threshold;
    return log2(value);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise natural logarithm
class Ln : public ElementWiseUnaryOp
{
  mirtkPublicAttributeMacro(double, Threshold);

public:

  Ln(double threshold = .01) : _Threshold(threshold) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &) const
  {
    if (value < _Threshold) value = _Threshold;
    return log(value);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Compute element-wise logarithm to base 10
class Lg : public ElementWiseUnaryOp
{
  mirtkPublicAttributeMacro(double, Threshold);

public:

  Lg(double threshold = .01) : _Threshold(threshold) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &) const
  {
    if (value < _Threshold) value = _Threshold;
    return log10(value);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Element-wise addition
class Add : public ElementWiseBinaryOp
{
public:

  /// Constructor
  Add(double value) : ElementWiseBinaryOp(value) {}

  /// Constructor
  Add(const char *fname) : ElementWiseBinaryOp(fname) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, double constant, bool &) const
  {
    return value + constant;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseBinaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Element-wise subtraction
class Sub : public ElementWiseBinaryOp
{
public:

  /// Constructor
  Sub(double value) : ElementWiseBinaryOp(value) {}

  /// Constructor
  Sub(const char *fname) : ElementWiseBinaryOp(fname) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, double constant, bool &) const
  {
    return value - constant;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseBinaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Element-wise multiplication
class Mul : public ElementWiseBinaryOp
{
public:

  /// Constructor
  Mul(double value) : ElementWiseBinaryOp(value) {}

  /// Constructor
  Mul(const char *fname) : ElementWiseBinaryOp(fname) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, double constant, bool &) const
  {
    return value * constant;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseBinaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Element-wise division
class Div : public ElementWiseBinaryOp
{
public:

  /// Constructor
  Div(double value) : ElementWiseBinaryOp(value) {}

  /// Constructor
  Div(const char *fname) : ElementWiseBinaryOp(fname) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, double constant, bool &) const
  {
    return (constant != .0 ? value / constant : numeric_limits<double>::quiet_NaN());
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseBinaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Element-wise division
class DivWithZero : public ElementWiseBinaryOp
{
public:

  /// Constructor
  DivWithZero(const char *fname) : ElementWiseBinaryOp(fname) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, double constant, bool &) const
  {
    return (constant != .0 ? value / constant : .0);
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseBinaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask values
class Mask : public ElementWiseBinaryOp
{
public:

  /// Constructor
  Mask(double value) : ElementWiseBinaryOp(value) {}

  /// Constructor
  Mask(const char *fname) : ElementWiseBinaryOp(fname) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, double constant, bool &mask) const
  {
    if (( _FileName.empty() && (fequal(value, constant) || (IsNaN(value) && IsNaN(constant)))) ||
        (!_FileName.empty() && constant == .0)) {
      mask = false;
    }
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseBinaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask values below or above a specified lower/upper threshold
class MaskOutsideInterval : public ElementWiseUnaryOp
{
  /// Lower threshold value
  mirtkPublicAttributeMacro(double, LowerThreshold);

  /// Upper threshold value
  mirtkPublicAttributeMacro(double, UpperThreshold);

public:

  /// Constructor
  MaskOutsideInterval(double l, double u)
  :
    _LowerThreshold(l), _UpperThreshold(u)
  {
    if (_UpperThreshold < _LowerThreshold) swap(_LowerThreshold, _UpperThreshold);
  }

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (value < _LowerThreshold || value > _UpperThreshold) mask = false;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask values below, equal, or above a specified lower/upper threshold
class MaskOutsideOpenInterval : public ElementWiseUnaryOp
{
  /// Lower threshold value
  mirtkPublicAttributeMacro(double, LowerThreshold);

  /// Upper threshold value
  mirtkPublicAttributeMacro(double, UpperThreshold);

public:

  /// Constructor
  MaskOutsideOpenInterval(double l, double u)
  :
    _LowerThreshold(l), _UpperThreshold(u)
  {
    if (_UpperThreshold < _LowerThreshold) swap(_LowerThreshold, _UpperThreshold);
  }

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (value <= _LowerThreshold || value >= _UpperThreshold) mask = false;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask values inside closed interval
class MaskInsideInterval : public ElementWiseUnaryOp
{
  /// Lower threshold value
  mirtkPublicAttributeMacro(double, LowerThreshold);

  /// Upper threshold value
  mirtkPublicAttributeMacro(double, UpperThreshold);

public:

  /// Constructor
  MaskInsideInterval(double l, double u)
  :
    _LowerThreshold(l), _UpperThreshold(u)
  {
    if (_UpperThreshold < _LowerThreshold) swap(_LowerThreshold, _UpperThreshold);
  }

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (_LowerThreshold <= value && value <= _UpperThreshold) mask = false;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask values inside open interval
class MaskInsideOpenInterval : public ElementWiseUnaryOp
{
  /// Lower threshold value
  mirtkPublicAttributeMacro(double, LowerThreshold);

  /// Upper threshold value
  mirtkPublicAttributeMacro(double, UpperThreshold);

public:

  /// Constructor
  MaskInsideOpenInterval(double l, double u)
  :
    _LowerThreshold(l), _UpperThreshold(u)
  {
    if (_UpperThreshold < _LowerThreshold) swap(_LowerThreshold, _UpperThreshold);
  }

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (_LowerThreshold < value && value < _UpperThreshold) mask = false;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask even values (e.g., segmentation labels of right hemisphere; cf MAL 2012)
class MaskEvenValues : public ElementWiseUnaryOp
{
public:

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (static_cast<int>(value) % 2 == 0) mask = false;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Mask odd values (e.g., segmentation labels of left hemisphere; cf MAL 2012)
class MaskOddValues : public ElementWiseUnaryOp
{
public:

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (static_cast<int>(value) % 2 == 1) mask = false;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Clamp values below or equal a given threshold value
class LowerThreshold : public ElementWiseUnaryOp
{
  /// Lower threshold value
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  LowerThreshold(double value) : _Threshold(value) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (value <= _Threshold) value = _Threshold;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Clamp values below or equal a given threshold value
class UpperThreshold : public ElementWiseUnaryOp
{
  /// Upper threshold value
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  UpperThreshold(double value) : _Threshold(value) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if (value >= _Threshold) value = _Threshold;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Clamp values below or equal a given threshold value
class Clamp : public ElementWiseUnaryOp
{
  /// Lower threshold value
  mirtkPublicAttributeMacro(double, LowerThreshold);

  /// Upper threshold value
  mirtkPublicAttributeMacro(double, UpperThreshold);

public:

  /// Constructor
  Clamp(double l, double u) : _LowerThreshold(l), _UpperThreshold(u) {}

  /// Transform data value and/or mask data value by setting *mask = false
  virtual double Op(double value, bool &mask) const
  {
    if      (value <= _LowerThreshold) value = _LowerThreshold;
    else if (value >= _UpperThreshold) value = _UpperThreshold;
    return value;
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    ElementWiseUnaryOp::Process(n, data, mask);
    parallel_for(blocked_range<int>(0, n), *this);
  }
};

// -----------------------------------------------------------------------------
/// Rescale values to new range
class Rescale : public Op
{
private:

  double *_Data;
  bool   *_Mask;
  double  _Min;
  double  _Max;
  double  _Slope;
  double  _Intercept;

public:

  Rescale(double min = .0, double max = 1.0) : _Min(min), _Max(max) {}

  // Called by TBB's parallel_for, for internal use only!
  void operator()(const blocked_range<int> &re) const
  {
    double *data = _Data + re.begin();
    for (int i = re.begin(); i != re.end(); ++i, ++data) {
      *data = (*data * _Slope) + _Intercept;
    }
  }

  /// Process given data (not thread-safe!)
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    statistic::Extrema extrema;
    extrema.Evaluate(n, data, mask);
    _Slope     = (_Max - _Min) / (extrema.Max() - extrema.Min());
    _Intercept = _Min - (extrema.Min() * _Slope);
    _Data      = data;
    _Mask      = mask;
    parallel_for(blocked_range<int>(0, n), *this);
  }
};


} } } // namespace mirtk::data::op

#endif // MIRTK_DataFunctions_H
