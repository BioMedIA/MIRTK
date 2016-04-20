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

#include "mirtk/Vector.h"

#include "mirtk/Assert.h"
#include "mirtk/Memory.h"
#include "mirtk/Array.h"
#include "mirtk/Indent.h"
#include "mirtk/Cfstream.h"

#include "mirtk/NumericsConfig.h"
#if MIRTK_Numerics_WITH_MATLAB
#  include "mirtk/Matlab.h"
#endif


namespace mirtk {


// -----------------------------------------------------------------------------
void Vector::PermuteRows(Array<int> idx)
{
  mirtkAssert(idx.size() <= static_cast<size_t>(_rows), "valid permutation");
  for (int r1 = 0; r1 < static_cast<int>(idx.size()); ++r1) {
    int r2 = idx[r1];
    if (r2 == r1) continue;
    swap(_vector[r1], _vector[r2]);
    for (int r = r1 + 1; r < _rows; ++r) {
      if (idx[r] == r1) {
        swap(idx[r], idx[r1]);
        break;
      }
    }
  }
}

// -----------------------------------------------------------------------------
ostream &operator <<(ostream &os, const Vector &v)
{
  os << "irtkVector " << max(0, v._rows) << endl;
  if (v._rows > 0) {
    const bool swapped = (GetByteOrder() == LittleEndian);
    if (swapped) {
      swap64((char *)v._vector, (char *)v._vector, v._rows);
    }
    os.write((char *)v._vector, v._rows * sizeof(double));
    if (swapped) {
      swap64((char *)v._vector, (char *)v._vector, v._rows);
    }
  }
  return os;
}

// -----------------------------------------------------------------------------
istream &operator >> (istream &is, Vector &v)
{
  // Read file type identifier
  char keyword[11];
  is.read(keyword, 11);
  if (strncmp(keyword, "irtkVector ", 11) != 0) {
    cerr << "Vector: File does not appear to be a binary (M)IRTK Vector file" << endl;
    exit(1);
  }

  // Read vector size
  int rows = 0;
  is >> rows;
  if (!is || rows < 0) {
    cerr << "Vector: Failed to read vector size from binary file" << endl;
    exit(1);
  }

  // Skip remaining characters (e.g., comment) on first line
  char buffer[256];
  while (is.get() != '\n')  {
    is.get(buffer, 256);
  }

  // Read vector data
  v.Initialize(rows);
  if (rows > 0) {
    is.read((char *)v._vector, rows * sizeof(double));
    if (is.fail()) {
      cerr << "Vector: File contains fewer values than expected: #rows = " << rows << endl;
      exit(1);
    }
    if (GetByteOrder() == LittleEndian) {
      swap64((char *)v._vector, (char *)v._vector, v._rows);
    }
  }

  return is;
}

// -----------------------------------------------------------------------------
Cofstream &operator <<(Cofstream& to, const Vector &v)
{
  to.WriteAsChar("irtkVector", 11);
  const int rows = max(0, v._rows);
  to.WriteAsInt(&rows, 1);
  if (v._rows > 0) to.WriteAsDouble(v._vector, v._rows);
  return to;
}

// -----------------------------------------------------------------------------
Cifstream &operator >>(Cifstream &from, Vector &v)
{
  char keyword[11];
  from.ReadAsChar(keyword, 11);
  keyword[10] = '\0';
  if (strcmp(keyword, "irtkVector") != 0) {
    cerr << "Vector: File does not appear to be a binary (M)IRTK Vector file" << endl;
    exit(1);
  }

  int rows = 0;
  if (!from.ReadAsInt(&rows, 1) || rows < 0) {
    cerr << "Vector: Failed to read vector size from binary file" << endl;
    exit(1);
  }

  v.Initialize(rows);
  if (rows > 0 && !from.ReadAsDouble(v._vector, rows)) {
    cerr << "Vector: File contains fewer values than expected: #rows = " << rows << endl;
    exit(1);
  }

  return from;
}

// -----------------------------------------------------------------------------
void Vector::Print(Indent indent) const
{
  cout << indent << "Vector " << _rows << endl;
  ++indent;
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  for (int i = 0; i < _rows; i++) {
    cout << indent << setw(15) << _vector[i] << endl;
  }
  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

// -----------------------------------------------------------------------------
void Vector::Read(const char *filename)
{
  // Open file stream
  ifstream from(filename, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "Vector::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read vector
  from >> *this;
}

// -----------------------------------------------------------------------------
void Vector::Write(const char *filename) const
{
  // Open file stream
  ofstream to(filename, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "Vector::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write vector
  to << *this;
}

#if MIRTK_Numerics_WITH_MATLAB

// -----------------------------------------------------------------------------
inline mxArray *VectorToMxArray(const Vector &v)
{
  Matlab::Initialize();
  mxArray *m = mxCreateDoubleMatrix(v.Rows(), 1, mxREAL);
  memcpy(mxGetPr(m), v.RawPointer(), v.Rows() * sizeof(double));
  return m;
}

// -----------------------------------------------------------------------------
bool Vector::WriteMAT(const char *fname, const char *varname) const
{
  Matlab::Initialize();
  MATFile *fp = matOpen(fname, "w");
  if (fp == NULL) return false;
  mxArray *m = VectorToMxArray(*this);
  if (matPutVariable(fp, varname, m) != 0) {
    mxDestroyArray(m);
    matClose(fp);
    return false;
  }
  mxDestroyArray(m);
  return (matClose(fp) == 0);
}

#endif // MIRTK_Numerics_WITH_MATLAB


} // namespace mirtk
