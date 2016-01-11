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

#ifndef MIRTK_DataStatistics_H
#define MIRTK_DataStatistics_H

#include <mirtkDataOp.h>

#include <mirtkMath.h>
#include <mirtkString.h>
#include <mirtkStream.h>
#include <mirtkArray.h>
#include <mirtkAlgorithm.h> // partial_sort


namespace mirtk { namespace data { namespace statistic {


// -----------------------------------------------------------------------------
/// Base class of all data statistics
class Statistic : public Op
{
  /// Description of data statistic
  mirtkPublicAttributeMacro(string, Description);

  /// Column names for (CSV) output
  mirtkPublicAttributeMacro(Array<string>, Names);

  /// Values of data statistic
  mirtkReadOnlyAttributeMacro(Array<double>, Values);

protected:

  /// Constructor
  Statistic(const char *desc = NULL, const Array<string> *names = NULL)
  :
    _Description(desc ? desc : "Unknown statistic"),
    _Values(1, numeric_limits<double>::quiet_NaN())
  {
    if (names) _Names = *names;
  }

  /// Set value of statistic (first entry of _Values vector)
  void Value(double v)
  {
    _Values[0] = v;
  }

public:

  /// Destructor
  virtual ~Statistic() {}

  /// Get value of statistic (first entry of _Values vector)
  virtual double Value() const
  {
    return _Values[0];
  }

  /// Process given data
  virtual void Process(int n, double *data, bool *mask = NULL)
  {
    this->Evaluate(n, data, mask);
  }

#if MIRTK_Image_WITH_VTK

  /// Process given vtkDataArray
  virtual void Process(vtkDataArray *data, bool *mask = NULL)
  {
    const int n = static_cast<int>(data->GetNumberOfTuples() * data->GetNumberOfComponents());
    if (data->GetDataType() == VTK_DOUBLE) {
      this->Process(n, reinterpret_cast<double *>(data->GetVoidPointer(0)), mask);
    } else {
      double *d = new double[n];
      double *tuple = d;
      for (vtkIdType i = 0; i < data->GetNumberOfTuples(); ++i) {
        data->GetTuple(i, tuple);
        tuple += data->GetNumberOfComponents();
      }
      this->Process(n, d, mask);
      // Unlike base class, no need to copy data back to input array
      delete[] d;
    }
  }

#endif // MIRTK_Image_WITH_VTK

  /// Evaluate statistic for given data
  virtual void Evaluate(int, const double *, const bool * = NULL) = 0;

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

// -----------------------------------------------------------------------------
/// Get minimum value
class Min : public Statistic
{
public:

  Min(const char *desc = "Minimum", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Min";
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    int m = 0;
    double d, v = numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        ++m;
        d = static_cast<double>(data[i]);
        if (d < v) v = d;
      }
    }
    return (m > 0 ? v : numeric_limits<double>::quiet_NaN());
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(n, data, mask));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get minimum absolute value
class MinAbs : public Statistic
{
public:

  MinAbs(const char *desc = "Minimum absolute value", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Min abs value";
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    int m = 0;
    double d, v = numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        ++m;
        d = abs(static_cast<double>(data[i]));
        if (d < v) v = d;
      }
    }
    return (m > 0 ? v : numeric_limits<double>::quiet_NaN());
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(n, data, mask));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get maximum value
class Max : public Statistic
{
public:

  Max(const char *desc = "Maximum", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Max";
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    int m = 0;
    double d, v = -numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        ++m;
        d = static_cast<double>(data[i]);
        if (d > v) v = d;
      }
    }
    return (m > 0 ? v : numeric_limits<double>::quiet_NaN());
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(n, data, mask));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get maximum absolute value
class MaxAbs : public Statistic
{
public:

  MaxAbs(const char *desc = "Maximum absolute value", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Max abs value";
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    int m = 0;
    double d, v = -numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        ++m;
        d = abs(static_cast<double>(data[i]));
        if (d > v) v = d;
      }
    }
    return (m > 0 ? v : numeric_limits<double>::quiet_NaN());
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(n, data, mask));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Get minimum and maximum values
class Extrema : public Statistic
{
public:

  Extrema(const char *desc = "Extrema", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(2);
      _Names[0] = "Min";
      _Names[1] = "Max";
    }
  }

  template <class T>
  static void Calculate(double &v1, double &v2, int n, const T *data, const bool *mask = NULL)
  {
    double d;
    v1 = +numeric_limits<double>::infinity();
    v2 = -numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        d = static_cast<double>(data[i]);
        if (d < v1) v1 = d;
        if (d > v2) v2 = d;
      }
    }
    if (v1 > v2) {
      v1 = v2 = numeric_limits<double>::quiet_NaN();
    }
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    _Values.resize(2);
    Calculate(_Values[0], _Values[1], n, data, mask);
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

  Range(const char *desc = "Range", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Range";
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    double v1, v2;
    Extrema::Calculate(v1, v2, n, data, mask);
    return v2 - v1;
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(n, data, mask));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Robust mean/average evaluation
class Mean : public Statistic
{
public:

  Mean(const char *desc = "Mean", const Array<string> *names = NULL)
  :
    Statistic(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Mean";
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    int    m =  0;
    double v = .0;

    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        ++m;
        v += (static_cast<double>(data[i]) - v) / m;
      }
    }

    return (m < 1 ? numeric_limits<double>::quiet_NaN() : v);
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(n, data, mask));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
};

// -----------------------------------------------------------------------------
/// Robust variance evaluation
class Var : public Statistic
{
  /// Mean value computed as a side-product of variance evaluation
  mirtkReadOnlyAttributeMacro(double, Mean);

public:

  Var(const char *desc = "Variance", const Array<string> *names = NULL)
  :
    Statistic(desc, names), _Mean(.0)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Var";
    }
  }

  template <class T>
  static void Calculate(double &mean, double &var, int n, const T *data, const bool *mask = NULL)
  {
    int m =  0;
    mean = var = .0;

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
      mean = var = numeric_limits<double>::quiet_NaN();
    } else if (m < 2) {
      var = .0;
    } else {
      var /= m - 1;
    }
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    double mean, var;
    Calculate(mean, var, n, data, mask);
    return var;
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Calculate(_Mean, _Values[0], n, data, mask);
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
  mirtkCalculateVtkDataArray2();
};

// -----------------------------------------------------------------------------
/// Robust evaluation of standard deviation
class StDev : public Var
{
public:

  StDev(const char *desc = "Standard deviation", const Array<string> *names = NULL)
  :
    Var(desc, names)
  {
    if (!names) {
      _Names.resize(1);
      _Names[0] = "Sigma";
    }
  }

  template <class T>
  static void Calculate(double &mean, double &sigma, int n, const T *data, const bool *mask = NULL)
  {
    Var::Calculate(mean, sigma, n, data, mask);
    sigma = sqrt(sigma);
  }

  template <class T>
  static double Calculate(int n, const T *data, const bool *mask = NULL)
  {
    double mean, sigma;
    Calculate(mean, sigma, n, data, mask);
    return sigma;
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Var::Evaluate(n, data, mask);
    Value(sqrt(this->Value()));
  }

  // Add support for vtkDataArray argument
  mirtkCalculateVtkDataArray1();
  mirtkCalculateVtkDataArray2();
};

// -----------------------------------------------------------------------------
/// Robust evaluation of Gaussian mean and standard deviation
class NormalDistribution : public StDev
{
public:

  NormalDistribution(const char *desc = "Normal distribution", const Array<string> *names = NULL)
  :
    StDev(desc, names)
  {
    if (!names) {
      _Names.resize(2);
      _Names[0] = "Mean";
      _Names[1] = "Sigma";
    }
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    StDev::Evaluate(n, data, mask);
    _Values.resize(2);
    _Values[1] = _Values[0];
    _Values[0] = _Mean;
  }

  double Mean()  const { return _Values[0]; }
  double Sigma() const { return _Values[1]; }
};

// -----------------------------------------------------------------------------
/// Percentile calculation
class Percentile : public Statistic
{
  /// Percentage of values that are lower than the percentile to look for
  mirtkReadOnlyAttributeMacro(int, P);

  /// Percentile rank computed as a side product
  mirtkReadOnlyAttributeMacro(double, Rank);

public:

  Percentile(int p, const char *desc = NULL, const Array<string> *names = NULL)
  :
    _P(p), _Rank(-1)
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
      _Names.resize(1);
      _Names[0] = ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "%";
    }

  }

  template <class T>
  static double Calculate(int p, double &rank, int n, const T *data, const bool *mask = NULL)
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
      rank = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }

    // Compute percentile rank
    rank = (double(p) / 100.0) * double(m + 1);

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
  static double Calculate(int p, int n, const T *data, const bool *mask = NULL)
  {
    double rank;
    return Calculate(p, rank, n, data, mask);
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(_P, _Rank, n, data, mask));
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, double &rank, vtkDataArray *data, const bool *mask = NULL)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, rank, n, reinterpret_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, rank, n, reinterpret_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, rank, n, reinterpret_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, rank, n, reinterpret_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    rank = numeric_limits<double>::quiet_NaN();
    return numeric_limits<double>::quiet_NaN();
  }

  static double Calculate(int p, vtkDataArray *data, const bool *mask = NULL)
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

  /// Percentile rank computed as a side product
  mirtkReadOnlyAttributeMacro(double, Rank);

public:

  AbsPercentile(int p, const char *desc = NULL, const Array<string> *names = NULL)
  :
    _P(p), _Rank(-1)
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
      _Names.resize(1);
      _Names[0] = ToString(p);
      if      (p == 1) _Names[0] += "st";
      else if (p == 2) _Names[0] += "nd";
      else if (p == 3) _Names[0] += "rd";
      else             _Names[0] += "th";
      _Names[0] += "% (abs)";
    }
  }

  template <class T>
  static double Calculate(int p, double &rank, int n, const T *data, const bool *mask = NULL)
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
      rank = numeric_limits<double>::quiet_NaN();
      return numeric_limits<double>::quiet_NaN();
    }

    // Compute percentile rank
    rank = (double(p) / 100.0) * double(m + 1);

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
  static double Calculate(int p, int n, const T *data, const bool *mask = NULL)
  {
    double rank;
    return Calculate(p, rank, n, data, mask);
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(_P, _Rank, n, data, mask));
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, double &rank, vtkDataArray *data, const bool *mask = NULL)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, rank, n, reinterpret_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, rank, n, reinterpret_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, rank, n, reinterpret_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, rank, n, reinterpret_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    rank = numeric_limits<double>::quiet_NaN();
    return numeric_limits<double>::quiet_NaN();
  }

  static double Calculate(int p, vtkDataArray *data, const bool *mask = NULL)
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
  /// Mean of values below the P-th percentile
  mirtkReadOnlyAttributeMacro(double, Mean);

public:

  LowerPercentileMean(int p, const char *desc = NULL, const Array<string> *names = NULL)
  :
    Percentile(p), _Mean(numeric_limits<double>::quiet_NaN())
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
      _Names.resize(1);
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
  static double Calculate(int p, int n, const T *data, const bool *mask = NULL)
  {
    const double threshold = Percentile::Calculate(p, n, data, mask);

    int    m =  0;
    double v = .0, d;
    
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        d = static_cast<double>(data[i]);
        if (d <= threshold) {
          ++m;
          v += (d - v) / m;
        }
      }
    }

    if (m < 1) return numeric_limits<double>::quiet_NaN();
    return v;
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(_P, n, data, mask));
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, vtkDataArray *data, const bool *mask = NULL)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, n, reinterpret_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, n, reinterpret_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, n, reinterpret_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, n, reinterpret_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    return numeric_limits<double>::quiet_NaN();
  }

#endif // MIRTK_Image_WITH_VTK
};

// -----------------------------------------------------------------------------
/// Upper percentile mean calculation
class UpperPercentileMean : public Percentile
{
  /// Mean of values above the P-th percentile
  mirtkReadOnlyAttributeMacro(double, Mean);
  
public:
  
  UpperPercentileMean(int p, const char *desc = NULL, const Array<string> *names = NULL)
  :
    Percentile(p), _Mean(numeric_limits<double>::quiet_NaN())
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
      _Names.resize(1);
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
  static double Calculate(int p, int n, const T *data, const bool *mask = NULL)
  {
    const double threshold = Percentile::Calculate(p, n, data, mask);

    int    m =  0;
    double v = .0, d;
    
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        d = static_cast<double>(data[i]);
        if (d >= threshold) {
          ++m;
          v += (d - v) / m;
        }
      }
    }

    if (m < 1) return numeric_limits<double>::quiet_NaN();
    return v;
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(_P, n, data, mask));
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, vtkDataArray *data, const bool *mask = NULL)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, n, reinterpret_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, n, reinterpret_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, n, reinterpret_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, n, reinterpret_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    return numeric_limits<double>::quiet_NaN();
  }

#endif // MIRTK_Image_WITH_VTK
};

// -----------------------------------------------------------------------------
/// Robust mean calculation, ignoring values below and above a certain percentile
class RobustMean : public Percentile
{
  /// Mean of values between the P-th and (100-P)-th percentile
  mirtkReadOnlyAttributeMacro(double, Mean);
  
public:
  
  RobustMean(int p, const char *desc = NULL, const Array<string> *names = NULL)
  :
    Percentile(p), _Mean(numeric_limits<double>::quiet_NaN())
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
      _Names.resize(1);
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
  static double Calculate(int p, int n, const T *data, const bool *mask = NULL)
  {
    const double min = Percentile::Calculate(        p, n, data, mask);
    const double max = Percentile::Calculate(100.0 - p, n, data, mask);

    int    m =  0;
    double v = .0, d;
    
    for (int i = 0; i < n; ++i) {
      if (!mask || mask[i]) {
        d = static_cast<double>(data[i]);
        if (min <= d && d <= max) {
          ++m;
          v += (d - v) / m;
        }
      }
    }

    if (m < 1) return numeric_limits<double>::quiet_NaN();
    return v;
  }

  void Evaluate(int n, const double *data, const bool *mask = NULL)
  {
    Value(Calculate(_P, n, data, mask));
  }

  // Add support for vtkDataArray
#if MIRTK_Image_WITH_VTK

  static double Calculate(int p, vtkDataArray *data, const bool *mask = NULL)
  {
    const int   n   = static_cast<int>(data->GetNumberOfTuples());
    const void *ptr = data->GetVoidPointer(0);
    switch (data->GetDataType()) {
      case VTK_SHORT:  return Calculate(p, n, reinterpret_cast<const short  *>(ptr), mask);
      case VTK_INT:    return Calculate(p, n, reinterpret_cast<const int    *>(ptr), mask);
      case VTK_FLOAT:  return Calculate(p, n, reinterpret_cast<const float  *>(ptr), mask);
      case VTK_DOUBLE: return Calculate(p, n, reinterpret_cast<const double *>(ptr), mask);
      default:
        cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl;
        exit(1);
    }
    return numeric_limits<double>::quiet_NaN();
  }

#endif // MIRTK_Image_WITH_VTK
};


} } } // namespace mirtk::data::statistic

#endif // MIRTK_DataStatistics_H

