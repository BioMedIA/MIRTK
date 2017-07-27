/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_DataStatistics_H
#define MIRTK_DataStatistics_H

#include "mirtk/DataOp.h"

#include "mirtk/Math.h"
#include "mirtk/String.h"
#include "mirtk/Stream.h"
#include "mirtk/Array.h"
#include "mirtk/Algorithm.h" // partial_sort


namespace mirtk { namespace data {


// =============================================================================
// Base class
// =============================================================================

// -----------------------------------------------------------------------------
/// Base class of all data statistics
class Statistic : public Op
{
  /// Do not include statistic in output report
  mirtkPublicAttributeMacro(bool, Hidden);

  /// Description of data statistic
  mirtkPublicAttributeMacro(string, Description);

  /// Column names for (CSV) output
  mirtkPublicAttributeMacro(Array<string>, Names);

  /// Values of data statistic
  mirtkReadOnlyAttributeMacro(Array<double>, Values);

protected:

  /// Constructor
  Statistic(int nvalues, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    _Hidden(false),
    _Description(desc ? desc : "Unknown statistic")
  {
    if (nvalues <= 0) nvalues = 1;
    _Names .resize(nvalues);
    _Values.resize(nvalues, NaN);
    if (names) {
      for (int i = 0; i < nvalues; ++i) {
        _Names[i] = names->at(i);
      }
    }
  }

  /// Set value of statistic (first entry of _Values vector)
  void Value(double v)
  {
    _Values[0] = v;
  }

public:

  /// Destructor
  virtual ~Statistic() {}

  /// Hide/Show statistic in output report
  mirtkOnOffMacro(Hidden);

  /// Get value of statistic (first entry of _Values vector)
  ///
  /// @returns Reference to first element of _Values vector such that the
  ///          address of the this element can be passed on to subsequent
  ///          data operations.
  virtual const double &Value() const
  {
    return _Values[0];
  }

  /// Process given data
  virtual void Process(int n, double *data, bool *mask = nullptr)
  {
    this->Evaluate(_Values, n, data, mask);
  }

#if MIRTK_Image_WITH_VTK

  /// Process given vtkDataArray
  virtual void Process(vtkDataArray *data, bool *mask = nullptr)
  {
    const int n = static_cast<int>(data->GetNumberOfTuples() * data->GetNumberOfComponents());
    if (data->GetDataType() == VTK_DOUBLE) {
      this->Process(n, reinterpret_cast<double *>(data->GetVoidPointer(0)), mask);
    } else {
      UniquePtr<double[]> _data(new double[n]);
      double *tuple = _data.get();
      for (vtkIdType i = 0; i < data->GetNumberOfTuples(); ++i) {
        data->GetTuple(i, tuple);
        tuple += data->GetNumberOfComponents();
      }
      this->Process(n, _data.get(), mask);
      // Unlike base class, no need to copy data back to input array
    }
  }

#endif // MIRTK_Image_WITH_VTK

  /// Evaluate statistic for given data
  virtual void Evaluate(Array<double> &, int, const double *, const bool * = nullptr) const = 0;

  /// Print column names of statistic values to output stream
  virtual void PrintHeader(ostream &os = cout, const char *delimiter = ",") const
  {
    for (size_t i = 0; i < _Values.size(); ++i) {
      if (i > 0) os << delimiter;
      if (i < _Names.size()) os << _Names[i];
    }
  }

  /// Print delimited statistic values to output stream
  virtual void PrintValues(ostream &os = cout, int digits = 5, const char *delimiter = ",") const
  {
    streamsize prev_precision = os.precision(digits);
    for (size_t i = 0; i < _Values.size(); ++i) {
      if (i > 0) os << delimiter;
      os << _Values[i];
    }
    os.precision(prev_precision);
  }

  /// Print statistic to output stream as "<desc> = <value>"
  virtual void Print(ostream &os = cout, int digits = 5, const char *prefix = "") const
  {
    streamsize prev_precision = os.precision(digits);
    if (prefix && prefix[0] != '\0') {
      os << prefix << ' ';
      os << char(tolower(_Description[0]));
      os << _Description.substr(1);
    } else {
      os << _Description;
    }
    os << " = ";
    if (_Values.size() != 1) os << "[";
    for (size_t i = 0; i < _Values.size(); ++i) {
      if (i > 0) os << ", ";
      os << _Values[i];
    }
    if (_Values.size() != 1) os << "]";
    os << endl;
    os.precision(prev_precision);
  }

};


// =============================================================================
// Data statistics
// =============================================================================

namespace statistic {


// -----------------------------------------------------------------------------
/// Count number of masked values (mask != false)
class Count : public Statistic
{
public:

  Count(const char *desc = "N", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "N";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    if (mask) {
      int count = 0;
      for (int i = 0; i < n; ++i) {
        if (mask[i]) ++count;
      }
      return count;
    } else {
      return n;
    }
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  /// Print delimited statistic values to output stream
  virtual void PrintValues(ostream &os = cout, int = 5, const char *delimiter = ",") const
  {
    for (size_t i = 0; i < _Values.size(); ++i) {
      if (i > 0) os << delimiter;
      os << int(_Values[i]);
    }
  }

  /// Print statistic to output stream as "<desc> = <value>"
  virtual void Print(ostream &os = cout, int = 5, const char *prefix = "") const
  {
    if (prefix && prefix[0] != '\0') {
      os << prefix << ' ';
      os << char(tolower(_Description[0]));
      os << _Description.substr(1);
    } else {
      os << _Description;
    }
    os << " = ";
    if (_Values.size() != 1) os << "[";
    for (size_t i = 0; i < _Values.size(); ++i) {
      if (i > 0) os << ", ";
      os << int(_Values[i]);
    }
    if (_Values.size() != 1) os << "]";
    os << endl;
  }
};

// -----------------------------------------------------------------------------
/// Get sum of values
class Sum : public Statistic
{
public:

  Sum(const char *desc = "Sum", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Sum";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double sum = 0.;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          sum += static_cast<double>(data[i]);
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        sum += static_cast<double>(data[i]);
      }
    }
    return sum;
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get minimum value
class Min : public Statistic
{
public:

  Min(const char *desc = "Minimum", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Min";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double d, v = inf;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = static_cast<double>(data[i]);
          if (d < v) v = d;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = static_cast<double>(data[i]);
        if (d < v) v = d;
      }
    }
    return (IsInf(v) ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get minimum absolute value
class MinAbs : public Statistic
{
public:

  MinAbs(const char *desc = "Minimum absolute value", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Min abs value";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double d, v = inf;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = abs(static_cast<double>(data[i]));
          if (d < v) v = d;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = abs(static_cast<double>(data[i]));
        if (d < v) v = d;
      }
    }
    return (IsInf(v) ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get maximum value
class Max : public Statistic
{
public:

  Max(const char *desc = "Maximum", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Max";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double d, v = -inf;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = static_cast<double>(data[i]);
          if (d > v) v = d;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = static_cast<double>(data[i]);
        if (d > v) v = d;
      }
    }
    return (IsInf(v) ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get maximum absolute value
class MaxAbs : public Statistic
{
public:

  MaxAbs(const char *desc = "Maximum absolute value", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Max abs value";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double d, v = -inf;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = abs(static_cast<double>(data[i]));
          if (d > v) v = d;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = abs(static_cast<double>(data[i]));
        if (d > v) v = d;
      }
    }
    return (IsInf(v) ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get minimum and maximum values
class Extrema : public Statistic
{
public:

  Extrema(const char *desc = "Extrema", const Array<string> *names = nullptr)
  :
    Statistic(2, desc, names)
  {
    if (!names) {
      _Names[0] = "Min";
      _Names[1] = "Max";
    }
  }

  template <class T>
  static void Calculate(double &v1, double &v2, int n, const T *data, const bool *mask = nullptr)
  {
    double d;
    v1 = +inf;
    v2 = -inf;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = static_cast<double>(data[i]);
          if (d < v1) v1 = d;
          if (d > v2) v2 = d;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = static_cast<double>(data[i]);
        if (d < v1) v1 = d;
        if (d > v2) v2 = d;
      }
    }
    if (v1 > v2) {
      v1 = v2 = NaN;
    }
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(2);
    Calculate(values[0], values[1], n, data, mask);
  }

  double Min() const { return _Values[0]; }
  double Max() const { return _Values[1]; }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray2();
};

// -----------------------------------------------------------------------------
/// Get value range
class Range : public Statistic
{
public:

  Range(const char *desc = "Range", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Range";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double v1, v2;
    Extrema::Calculate(v1, v2, n, data, mask);
    return v2 - v1;
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Median value (50th percentile)
class Median : public Statistic
{
public:

  Median(const char *desc = "Median", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Median";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    Array<T> values;
    if (mask) {
      values.reserve(n);
      for (int i = 0; i < n; ++i) {
        if (mask[i]) values.push_back(data[i]);
      }
    } else {
      values.resize(n);
      for (int i = 0; i < n; ++i) {
        values[i] = data[i];
      }
    }
    if (values.size() <= 0) return NaN;
    return NthElement(values, static_cast<int>(values.size())/2);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Robust mean/average evaluation
class Mean : public Statistic
{
public:

  Mean(const char *desc = "Mean", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) _Names[0] = "Mean";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    int    m = 0;
    double v = 0.;

    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          ++m, v += (static_cast<double>(data[i]) - v) / m;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        v += (static_cast<double>(data[i]) - v) / (i + 1);
      }
      m = n;
    }

    return (m < 1 ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Base class of statistics which compute mean, variance, and/or standard deviation
class MeanVar : public Statistic
{
protected:

  /// Constructor
  MeanVar(int nvalues, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    Statistic(nvalues, desc, names)
  {}

public:

  template <class T>
  static void Calculate(double &mean, double &var, int n, const T *data, const bool *mask = nullptr)
  {
    int m =  0;
    mean = var = 0.;

    double delta, d;
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        ++m;
        d = static_cast<double>(data[i]);
        delta = d - mean;
        mean += delta / m;
        var  += delta * (d - mean);
      }
    }

    if (m < 1) {
      mean = var = NaN;
    } else if (m < 2) {
      var = 0.;
    } else {
      var /= m - 1;
    }
  }
};

// -----------------------------------------------------------------------------
/// Robust variance evaluation
class Var : public MeanVar
{
public:

  Var(const char *desc = "Variance", const Array<string> *names = nullptr)
  :
    MeanVar(1, desc, names)
  {
    if (!names) _Names[0] = "Var";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double mean, var;
    MeanVar::Calculate(mean, var, n, data, mask);
    return var;
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    double mean;
    values.resize(1);
    MeanVar::Calculate(mean, values[0], n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Robust evaluation of standard deviation
class StDev : public MeanVar
{
public:

  StDev(const char *desc = "Standard deviation", const Array<string> *names = nullptr)
  :
    MeanVar(1, desc, names)
  {
    if (!names) _Names[0] = "Sigma";
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double mean, var;
    MeanVar::Calculate(mean, var, n, data, mask);
    return sqrt(var);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    double mean, var;
    MeanVar::Calculate(mean, var, n, data, mask);
    values.resize(1);
    values[0] = sqrt(var);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Robust evaluation of Gaussian mean and standard deviation
class NormalDistribution : public MeanVar
{
public:

  NormalDistribution(const char *desc = "Normal distribution", const Array<string> *names = nullptr)
  :
    MeanVar(2, desc, names)
  {
    if (!names) {
      _Names[0] = "Mean";
      _Names[1] = "Sigma";
    }
  }

  template <class T>
  static void Calculate(double &mean, double &sigma, int n, const T *data, const bool *mask = nullptr)
  {
    MeanVar::Calculate(mean, sigma, n, data, mask);
    sigma = sqrt(sigma);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(2);
    Calculate(values[0], values[1], n, data, mask);
  }

  double Mean()  const { return _Values[0]; }
  double Sigma() const { return _Values[1]; }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray2();
};

// -----------------------------------------------------------------------------
/// Mean absolute difference/deviation around the specified value
class AverageAbsoluteDifference : public Statistic
{
  /// Mean value
  mirtkPublicAttributeMacro(double, Mean);

  /// Address of mean value storage
  ///
  /// This pointer can be set to the memory location of a double value that
  /// will be set by a preceeding data statistic operation such as
  /// data::statistic::Mean and data::statistic::Median. The result of this
  /// statistic calculation thus becomes the input parameter of this operation.
  mirtkPublicAggregateMacro(const double, MeanPointer);

public:

  AverageAbsoluteDifference(double mean, const double *mean_ptr,
                            const char *desc = "Average absolute difference",
                            const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names),
    _Mean(mean),
    _MeanPointer(mean_ptr)
  {
    if (!names) _Names[0] = "MAD";
  }

  AverageAbsoluteDifference(double mean)
  :
    AverageAbsoluteDifference(mean, nullptr)
  {}

  AverageAbsoluteDifference(const double *mean)
  :
    AverageAbsoluteDifference(0., mean)
  {}

  template <class T>
  static double Calculate(double mean, int n, const T *data, const bool *mask = nullptr)
  {
    int    num = 0;
    double mad = 0.;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          mad += abs(data[i] - mean);
          num += 1;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        mad += abs(data[i] - mean);
      }
      num = n;
    }
    return (num > 0 ? mad / num : 0.);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(_MeanPointer ? *_MeanPointer : _Mean, n, data, mask);
  }
};

// -----------------------------------------------------------------------------
/// Mean absolute difference/deviation around the mean
class MeanAbsoluteDifference : public Statistic
{
public:

  MeanAbsoluteDifference(const char *desc = "Mean absolute difference", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) {
      _Names[0] = "MAD mean";
    }
  }

  template <class T>
  static void Calculate(double &mean, double &mad, int n, const T *data, const bool *mask = nullptr)
  {
    int num = 0;
    mean = Mean::Calculate(n, data, mask);
    mad  = 0.;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          mad += abs(data[i] - mean);
          num += 1;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        mad += abs(data[i] - mean);
      }
      num = n;
    }
    if (num > 0) mad /= num;
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double mean, mad;
    Calculate(mean, mad, n, data, mask);
    return mad;
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Mean absolute difference/deviation around the median
class MedianAbsoluteDifference : public Statistic
{
public:

  MedianAbsoluteDifference(const char *desc = "Median absolute difference", const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names)
  {
    if (!names) {
      _Names[0] = "MAD median";
    }
  }

  template <class T>
  static void Calculate(double &median, double &mad, int n, const T *data, const bool *mask = nullptr)
  {
    int num = 0;
    median = Median::Calculate(n, data, mask);
    mad    = 0.;
    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          mad += abs(data[i] - median);
          num += 1;
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        mad += abs(data[i] - median);
      }
      num = n;
    }
    if (num > 0) mad /= num;
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = nullptr)
  {
    double median, mad;
    Calculate(median, mad, n, data, mask);
    return mad;
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Percentile calculation
class Percentile : public Statistic
{
  /// Percentage of values that are lower than the percentile to look for
  mirtkReadOnlyAttributeMacro(int, P);

public:

  Percentile(int p, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names), _P(p)
  {
    if (!desc) {
      _Description = ToString(p);
      if      (p == 1) _Description += "st";
      else if (p == 2) _Description += "nd";
      else if (p == 3) _Description += "rd";
      else             _Description += "th";
      _Description += " percentile";
    }
    if (!names) {
      _Names[0] = ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "%";
    }
  }

  template <class T>
  static double CalculateGivenSortedData(int p, double &rank, int n, const T *data)
  {
    if (n == 0) {
      rank = NaN;
      return NaN;
    }

    // Compute percentile rank
    rank = (double(p) / 100.) * double(n + 1);

    // Split rank into integer and decimal components
    int    k = int(rank);
    double d = rank - k;

    // Compute percentile value according to NIST method
    // (cf. http://en.wikipedia.org/wiki/Percentile#Definition_of_the_NIST_method )
    if (k == 0) {
      return Min::Calculate(n, data);
    } else if (k >= n) {
      return Max::Calculate(n, data);
    } else {
      return data[k - 1] + d * (data[k] - data[k - 1]);
    }
  }

  template <class T>
  static double CalculateGivenSortedData(int p, int n, const T *data)
  {
    double rank;
    return CalculateGivenSortedData(p, rank, n, data);
  }

  template <class T>
  static double Calculate(int p, double &rank, int n, const T *data, const bool *mask = nullptr)
  {
    // Determine number of unmasked values
    int m = n;
    if (mask) {
      m = 0;
      for (int i = 0; i < n; ++i) {
        if (mask[i]) ++m;
      }
    }

    if (m == 0) {
      rank = NaN;
      return NaN;
    }

    // Compute percentile rank
    rank = (double(p) / 100.) * double(m + 1);

    // Split rank into integer and decimal components
    int    k = int(rank);
    double d = rank - k;

    // Compute percentile value according to NIST method
    // (cf. http://en.wikipedia.org/wiki/Percentile#Definition_of_the_NIST_method )
    if (k == 0) {
      return Min::Calculate(n, data, mask);
    } else if (k >= m) {
      return Max::Calculate(n, data, mask);
    } else {
      // Copy unmasked data only
      double *v = new double[m];
      if (mask) {
        m = 0;
        for (int i = 0; i < n; ++i) {
          if (mask[i]) v[m++] = static_cast<double>(data[i]);
        }
      } else {
        for (int i = 0; i < n; ++i) {
          v[i] = static_cast<double>(data[i]);
        }
      }
      // Sort data below k + 1
      partial_sort(v, v + k + 1, v + m);
      // Compute percentile value
      double percentile = v[k - 1] + d * (v[k] - v[k - 1]);
      // Free (partially) sorted data copy
      delete[] v;
      return percentile;
    }
  }

  template <class T>
  static double Calculate(int p, int n, const T *data, const bool *mask = nullptr)
  {
    double rank;
    return Calculate(p, rank, n, data, mask);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(_P, n, data, mask);
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, double &rank, vtkDataArray *data, const bool *mask = nullptr)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, rank, n, static_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, rank, n, static_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, rank, n, static_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, rank, n, static_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    rank = NaN;
    return NaN;
  }

  static double Calculate(int p, vtkDataArray *data, const bool *mask = nullptr)
  {
    double rank;
    return Calculate(p, rank, data, mask);
  }

#endif // MIRTK_Image_WITH_VTK
};

// -----------------------------------------------------------------------------
/// Absolute value percentile calculation
class AbsPercentile : public Statistic
{
  /// Percentage of values that are lower than the percentile to look for
  mirtkReadOnlyAttributeMacro(int, P);

public:

  AbsPercentile(int p, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    Statistic(1, desc, names), _P(p)
  {
    if (!desc) {
      _Description = ToString(p);
      if      (p == 1) _Description += "st";
      else if (p == 2) _Description += "nd";
      else if (p == 3) _Description += "rd";
      else             _Description += "th";
      _Description += " absolute value percentile";
    }
    if (!names) {
      _Names[0] = ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "% (abs)";
    }
  }

  template <class T>
  static double Calculate(int p, double &rank, int n, const T *data, const bool *mask = nullptr)
  {
    // Determine number of unmasked values
    int m = n;
    if (mask) {
      m = 0;
      for (int i = 0; i < n; ++i) {
        if (mask[i]) ++m;
      }
    }

    if (m == 0) {
      rank = NaN;
      return NaN;
    }

    // Compute percentile rank
    rank = (double(p) / 100.) * double(m + 1);

    // Split rank into integer and decimal components
    int    k = int(rank);
    double d = rank - k;

    // Compute percentile value according to NIST method
    // (cf. http://en.wikipedia.org/wiki/Percentile#Definition_of_the_NIST_method )
    if (k == 0) {
      return MinAbs::Calculate(n, data, mask);
    } else if (k >= m) {
      return MaxAbs::Calculate(n, data, mask);
    } else {
      // Copy unmasked data only
      double *v = new double[m];
      if (mask) {
        m = 0;
        for (int i = 0; i < n; ++i) {
          if (mask[i]) v[m++] = abs(static_cast<double>(data[i]));
        }
      } else {
        for (int i = 0; i < n; ++i) {
          v[i] = abs(static_cast<double>(data[i]));
        }
      }
      // Sort data below k + 1
      partial_sort(v, v + k + 1, v + m);
      // Compute percentile value
      double percentile = v[k - 1] + d * (v[k] - v[k - 1]);
      // Free (partially) sorted data copy
      delete[] v;
      return percentile;
    }
  }

  template <class T>
  static double Calculate(int p, int n, const T *data, const bool *mask = nullptr)
  {
    double rank;
    return Calculate(p, rank, n, data, mask);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(_P, n, data, mask);
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, double &rank, vtkDataArray *data, const bool *mask = nullptr)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, rank, n, static_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, rank, n, static_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, rank, n, static_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, rank, n, static_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    rank = numeric_limits<double>::quiet_NaN();
    return numeric_limits<double>::quiet_NaN();
  }

  static double Calculate(int p, vtkDataArray *data, const bool *mask = nullptr)
  {
    double rank;
    return Calculate(p, rank, data, mask);
  }

#endif
};

// -----------------------------------------------------------------------------
/// Lower percentile mean calculation
class LowerPercentileMean : public Percentile
{
public:

  LowerPercentileMean(int p, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    Percentile(p, desc, names)
  {
    if (!desc) {
      _Description  = "Mean below ";
      _Description += ToString(p);
      if      (p == 1) _Description += "st";
      else if (p == 2) _Description += "nd";
      else if (p == 3) _Description += "rd";
      else             _Description += "th";
      _Description += " percentile";
    }
    if (!names) {
      _Names[0]  = "Mean <";
      _Names[0] += ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "%";
    }
  }

  template <class T>
  static double Calculate(int p, int n, const T *data, const bool *mask = nullptr)
  {
    const double threshold = Percentile::Calculate(p, n, data, mask);

    int    m = 0;
    double v = 0., d;

    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = static_cast<double>(data[i]);
          if (d <= threshold) {
            ++m;
            v += (d - v) / m;
          }
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = static_cast<double>(data[i]);
        if (d <= threshold) {
          ++m;
          v += (d - v) / m;
        }
      }
    }

    return (m < 1 ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(_P, n, data, mask);
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, vtkDataArray *data, const bool *mask = nullptr)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, n, static_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, n, static_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, n, static_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, n, static_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    return NaN;
  }

#endif // MIRTK_Image_WITH_VTK
};

// -----------------------------------------------------------------------------
/// Upper percentile mean calculation
class UpperPercentileMean : public Percentile
{
public:
  
  UpperPercentileMean(int p, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    Percentile(p, desc, names)
  {
    if (!desc) {
      _Description  = "Mean above ";
      _Description += ToString(p);
      if      (p == 1) _Description += "st";
      else if (p == 2) _Description += "nd";
      else if (p == 3) _Description += "rd";
      else             _Description += "th";
      _Description += " percentile";
    }
    if (!names) {
      _Names[0]  = "Mean >";
      _Names[0] += ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "%";
    }
  }

  template <class T>
  static double Calculate(int p, int n, const T *data, const bool *mask = nullptr)
  {
    const double threshold = Percentile::Calculate(p, n, data, mask);

    int    m = 0;
    double v = 0., d;

    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = static_cast<double>(data[i]);
          if (d >= threshold) {
            ++m;
            v += (d - v) / m;
          }
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = static_cast<double>(data[i]);
        if (d >= threshold) {
          ++m;
          v += (d - v) / m;
        }
      }
    }

    return (m < 1 ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(_P, n, data, mask);
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, vtkDataArray *data, const bool *mask = nullptr)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, n, static_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, n, static_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, n, static_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, n, static_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    return NaN;
  }

#endif // MIRTK_Image_WITH_VTK
};

// -----------------------------------------------------------------------------
/// Robust mean calculation, ignoring values below and above a certain percentile
class RobustMean : public Percentile
{
public:
  
  RobustMean(int p, const char *desc = nullptr, const Array<string> *names = nullptr)
  :
    Percentile(p, desc, names)
  {
    if (!desc) {
      _Description  = "Mean excl. ";
      _Description += ToString(p);
      if      (p == 1) _Description += "st";
      else if (p == 2) _Description += "nd";
      else if (p == 3) _Description += "rd";
      else             _Description += "th";
      _Description += " percentile";
    }
    if (!names) {
      _Names[0]  = "Mean <>";
      _Names[0] += ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "%";
    }
  }

  template <class T>
  static double Calculate(int p, int n, const T *data, const bool *mask = nullptr)
  {
    const double min = Percentile::Calculate(      p, n, data, mask);
    const double max = Percentile::Calculate(100 - p, n, data, mask);

    int    m = 0;
    double v = 0., d;

    if (mask) {
      for (int i = 0; i < n; ++i) {
        if (mask[i]) {
          d = static_cast<double>(data[i]);
          if (min <= d && d <= max) {
            ++m;
            v += (d - v) / m;
          }
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        d = static_cast<double>(data[i]);
        if (min <= d && d <= max) {
          ++m;
          v += (d - v) / m;
        }
      }
    }

    return (m < 1 ? NaN : v);
  }

  void Evaluate(Array<double> &values, int n, const double *data, const bool *mask = nullptr) const
  {
    values.resize(1);
    values[0] = Calculate(_P, n, data, mask);
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, vtkDataArray *data, const bool *mask = nullptr)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, n, static_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, n, static_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, n, static_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, n, static_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    return NaN;
  }

#endif // MIRTK_Image_WITH_VTK
};


} } } // namespace mirtk::data::statistic

#endif // MIRTK_DataStatistics_H

