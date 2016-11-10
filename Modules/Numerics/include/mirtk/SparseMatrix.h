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

#ifndef MIRTK_SparseMatrix_H
#define MIRTK_SparseMatrix_H

#include "mirtk/Config.h"
#include "mirtk/Object.h"

#include "mirtk/Assert.h"
#include "mirtk/Memory.h"
#include "mirtk/Vector.h"
#include "mirtk/Matrix.h"
#include "mirtk/Math.h"
#include "mirtk/Pair.h"
#include "mirtk/Array.h"
#include "mirtk/Algorithm.h"
#include "mirtk/Profiling.h"

#include "mirtk/NumericsConfig.h"
#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)
#  include "mirtk/Matlab.h"
#endif

#include <iterator> // distance


namespace mirtk {


/**
 * Sparse matrix with generic non-zero value type
 */
template <class TEntry>
class GenericSparseMatrix : public Object
{
  mirtkObjectMacro(SparseMatrix);
 
  // ---------------------------------------------------------------------------
  // Types
public:

  /// Enumeration of implemented storage layouts, i.e., either
  /// compressed row storage (CRS) or compressed column storage (CCS)
  enum StorageLayout { CRS, CCS };

  /// Type of matrix entry values
  typedef TEntry EntryType;

  /// Type of non-zero entry
  typedef Pair<int, TEntry> Entry;

  /// List of non-zero entries
  typedef Array<Entry> Entries;

  // ---------------------------------------------------------------------------
  // Attributes
 
  /// Storage layout
  mirtkReadOnlyAttributeMacro(enum StorageLayout, Layout);
 
  /// Number of rows
  mirtkReadOnlyAttributeMacro(int, Rows);

  /// Number of columns
  mirtkReadOnlyAttributeMacro(int, Cols);

  /// Number of non-zero entries
  mirtkReadOnlyAttributeMacro(int, NNZ);
 
  /// Number of allocated elements
  mirtkAttributeMacro(int, Size);

  /// Maximum number of unused entries, i.e., (_Size - _NNZ)
  mirtkPublicAttributeMacro(int, MaxNumberOfUnusedEntries);

protected:

  /// CRS: Index of first non-zero entry in i-th row
  /// CCS: Row indices of non-zero entries
  int *_Row;

  /// CRS: Column indices of non-zero entries
  /// CCS: Index of first non-zero entry in i-th column
  int *_Col;

  /// Non-zero entries
  EntryType *_Data;

  /// CRS: Indices of non-zero entries for each column
  /// CCS: Indices of non-zero entries for each row
  Array<int> *_Index;

  /// Copy other matrix, used by copy constructor and assignment operator only
  template <class TOtherEntry>
  void CopyAttributes(const GenericSparseMatrix<TOtherEntry> &);

  /// Check non-zero entries, sort them, remove zero entries, and sum up duplicates
  void CheckEntries(Entries &) const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Default constructor
  GenericSparseMatrix(StorageLayout = CCS);

  /// Construct sparse m x n matrix
  GenericSparseMatrix(int, int = 0, StorageLayout = CCS);

  /// Construct sparse m x n matrix with specified number of non-zero entries
  GenericSparseMatrix(int, int, int, StorageLayout = CCS);

  /// Copy constructor
  template <class TOtherEntry>
  explicit GenericSparseMatrix(const GenericSparseMatrix<TOtherEntry> &);

  /// Assignment operator
  template <class TOtherEntry>
  GenericSparseMatrix &operator =(const GenericSparseMatrix<TOtherEntry> &);

  /// Initialize sparse matrix
  void Initialize(int, int = 0, int = 0);

  /// Initialize sparse matrix given compressed rows (CRS) or columns (CCS)
  void Initialize(int, int, Array<Entries> &, bool = false);

  /// Initialize sparse matrix given compressed rows (CRS) or columns (CCS)
  void Initialize(int, int, Entries [], bool = false);

#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)

  /// Initialize sparse matrix from MATLAB array
  ///
  /// If the input MATLAB array is a sparse matrix, this function is most
  /// efficient when also this sparse matrix uses the CCS layout.
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_MATLAB is 1.
  void Initialize(const mxArray *);

#endif

  /// Destructor
  virtual ~GenericSparseMatrix();

  // ---------------------------------------------------------------------------
  // Entries

  /// Change storage layout
  void Layout(StorageLayout);

  /// Get raw access to index and data arrays
  int GetRawData(int *&, int *&, TEntry *&) const;

  /// Get raw access to non-zero values
  TEntry *RawPointer(int = 0);

  /// Get raw access to non-zero values
  const TEntry *RawPointer(int = 0) const;

  /// Reserve more memory
  void Reserve(int);

  /// Remove any zero entries
  void RemoveZeros();

  /// Number of non-zero entries in specified row
  int RowNNZ(int) const;

  /// Number of non-zero entries in specified column
  int ColNNZ(int) const;

  /// Number of non-zero entries in specified submatrix
  int SubNNZ(int, int, int, int) const;

  /// Precompute indices to speed up successive row and column operations.
  /// Not required for operations that only act on entries of a given row in
  /// case of the CRS layout or columns in case of CCS layout. Column (row)
  /// operations on a matrix in CRS (CCS) layout are significantly more
  /// efficient after a call to this function, which however requires O(NNZ)
  /// more memory.
  ///
  /// \note The precomputed index is automatically cleared whenever a new
  ///       non-zero entry is inserted or an existing one is removed.
  void Index();

  /// Free memory needed for precomputed indices
  void ClearIndex();

  /// Get reference to (non-zero) matrix entry
  ///
  /// \attention This operator always inserts a new entry if none exists.
  ///            Use Put to set an entry to zero and remove the previous
  ///            non-zero entry from the sparse matrix. Do not use to iterate
  ///            over all matrix entries as otherwise a new entry is added for
  ///            each matrix row and column index pair! Use Get instead.
  ///
  /// \remarks This function is not thread-safe if non-zero entry does not exist.
  EntryType &operator ()(int, int = -1);

  /// Set value
  /// \remarks This function is not thread-safe if non-zero entry does not exist.
  void Put(int, int, EntryType);

  /// Get value
  EntryType Get(int, int = -1) const;

  /// Get non-zero entries of specified row
  /// \note This function is most efficient when the CRS layout is used.
  void GetRow(int, Entries &) const;

  /// Set non-zero entries of specified row
  ///
  /// \note More efficient then Put when the CRS layout is used.
  ///       Implicitly sorts the input entries.
  ///
  /// \remarks This function is not thread-safe.
  void Row(int, Entries &, bool = false);

  /// Get non-zero entries of specified row
  /// \note This function is most efficient when the CRS layout is used.
  Entries Row(int) const;

  /// Get non-zero entries of specified column
  /// \note This function is most efficient when the CCS layout is used.
  void GetCol(int, Entries &) const;

  /// Set non-zero entries of specified column
  ///
  /// \note More efficient then Put when the CCS layout is used.
  ///       Implicitly sorts the input entries.
  ///
  /// \remarks This function is not thread-safe.
  void Col(int, Entries &, bool = false);

  /// Get non-zero entries of specified column
  /// \note This function is most efficient when the CCS layout is used.
  Entries Col(int) const;

  /// Get non-zero entries of specified column
  /// \note This function is most efficient when the CCS layout is used.
  void GetColumn(int, Entries &) const;

  /// Set non-zero entries of specified column
  /// \note More efficient then Put when the CCS layout is used.
  ///       Implicitly sorts the input entries.
  void Column(int, Entries &, bool = false);

  /// Get non-zero entries of specified column
  /// \note This function is most efficient when the CCS layout is used.
  Entries Column(int) const;

  /// Get diagonal values
  void GetDiag(Vector &d) const;

  /// Get diagonal values
  Vector Diag() const;

  /// Set diagonal to specified (non-zero) values
  void Diag(const Vector &d);

  /// Set diagonal to specified (non-zero) value
  void Diag(TEntry d);

  /// Get non-zero entries of specified submatrix
  GenericSparseMatrix Sub(int, int, int, int) const;

  /// Set non-zero entries of specified submatrix
  void Sub(int, int, const GenericSparseMatrix &);

  /// Free all non-zero entries
  virtual void Clear();

  /// Set all non-zero entries to zero (i.e., remove)
  void Zero();

  // ---------------------------------------------------------------------------
  // Common operations

  /// Transpose matrix
  ///
  /// \param[in] keep_layout Whether to keep the storage layout. If \c false,
  ///                        transposing a CRS matrix is achieved by simply
  ///                        changing to CCS and vice versa.
  void Transpose(bool keep_layout = false);

  /// Calculate sum of all entries in specified row
  /// \note This function is most efficient when the CRS layout is used.
  EntryType RowSum(int) const;

  /// Calculate sum of all entries in specified column
  /// \note This function is most efficient when the CCS layout is used.
  EntryType ColSum(int) const;

  /// Calculate sum of all entries in specified column
  /// \note This function is most efficient when the CCS layout is used.
  EntryType ColumnSum(int) const;

  /// Multiply matrix by vector, used to interface with ARPACK
  /// \note This function is most efficient when the CRS layout is used.
  void MultAv(EntryType [], EntryType []) const;

  /// Multiply by a scalar in-place
  GenericSparseMatrix &operator *=(EntryType);

  /// Divide by a scalar in-place
  GenericSparseMatrix &operator /=(EntryType);

  /// Multiply by a scalar
  GenericSparseMatrix operator *(EntryType) const;

  /// Divide by a scalar
  GenericSparseMatrix operator /(EntryType) const;

  /// Multiply specified row by a scalar
  /// \note This function is most efficient when the CRS layout is used.
  GenericSparseMatrix &ScaleRow(int, EntryType);

  /// Multiply specified column by a scalar
  /// \note This function is most efficient when the CCS layout is used.
  GenericSparseMatrix &ScaleCol(int, EntryType);

  /// Multiply specified column by a scalar
  /// \note This function is most efficient when the CCS layout is used.
  GenericSparseMatrix &ScaleColumn(int, EntryType);

  /// Whether this matrix is symmetric
  bool IsSymmetric() const;

  /// Make square matrix symmetric by adding its transpose and divide by 2
  ///
  /// \param[in] extent Whether to copy the entry of existing entries to their
  ///                   respective transpose entries without dividing by two.
  ///                   This option can be used to create a symmetric matrix
  ///                   from a lower or upper triangular matrix by simply
  ///                   adding the transpose of the triangular matrix to this
  ///                   (while leaving the diagonal entries unmodified).
  void MakeSymmetric(bool extent = false);

  // ---------------------------------------------------------------------------
  // Eigen decomposition

  /// Largest/Smallest eigenvalues
  ///
  /// Uses an iterative eigen solver to compute the \p k eigenvalues
  /// which have largest or smallest magnitude, respectively, or are closest in
  /// magnitude to the specified \p sigma value. Uses the Implicitly Restarted
  /// Arnoldi Method implemented by ARPACK. The LU factorization used for the
  /// shift-and-invert mode is computed using UMFPACK.
  ///
  /// \param[out] v     Converged eigenvalues.
  /// \param[in]  k     Number of requested eigenvalues.
  /// \param[in]  sigma Which eigenvalues of the spectrum, i.e., 'LM', 'LA', 'SM', 'SA',
  ///                   or numeric value sigma used for shift-and-invert mode.
  /// \param[in]  p     Number of Lanczos basis vectors.
  /// \param[in]  tol   Ritz estimate residual <= tol*norm(this).
  /// \param[in]  maxit Maximum number of iterations.
  /// \param[in]  v0    Starting vector. By default randomly generated.
  ///
  /// \returns Number of converged eigenvalues.
  ///
  /// \note Only implemented for real double precision sparse matrices.
  ///       Use only when MIRTK_Numerics_WITH_eigs is 1.
  int Eigenvalues(Vector &v, int k, const char *sigma = "LM",
                  int p = 0, double tol = .0, int maxit = 0, Vector *v0 = NULL) const;

  /// Eigenvectors of largest/smallest eigenvalues
  ///
  /// Uses an iterative eigen solver to compute the eigenvectors corresponding
  /// to the \p k eigenvalues which have largest or smallest magnitude,
  /// respectively, or are closest in magnitude to the specified \p sigma value.
  /// Uses the Implicitly Restarted Arnoldi Method (IRAM) implemented by ARPACK.
  /// The LU factorization for the shift-and-invert mode is computed using UMFPACK.
  ///
  /// \param[out]    E     Matrix of eigenvectors in columns.
  /// \param[in]     k     Number of requested eigenvalues.
  /// \param[in]     sigma Which eigenvalues of the spectrum, i.e., 'LM', 'LA', 'SM', 'SA',
  ///                      or numeric value sigma used for shift-and-invert mode.
  /// \param[in]     p     Number of Lanczos basis vectors.
  /// \param[in]     tol   Ritz estimate residual <= tol*norm(this).
  /// \param[in]     maxit Maximum number of iterations.
  /// \param[in,out] v0    Input:  Starting vector. By default randomly generated.
  ///                      Output: Final residual vector.
  ///
  /// \returns Number of converged eigenvalues.
  ///
  /// \note Only implemented for real double precision sparse matrices.
  ///       Use only when MIRTK_Numerics_WITH_eigs is 1.
  int Eigenvectors(Matrix &E, int k, const char *sigma = "LM",
                   int p = 0, double tol = .0, int maxit = 0, Vector *v0 = NULL) const;

  /// Largest/Smallest eigenvalues and -vectors
  ///
  /// Uses an iterative eigen solver to compute the \p k eigenvalues
  /// which have largest or smallest magnitude, respectively, or are closest in
  /// magnitude to the specified \p sigma value. Uses the Implicitly Restarted
  /// Arnoldi Method implemented by ARPACK. The LU factorization used for the
  /// shift-and-invert mode is computed using UMFPACK.
  ///
  /// \param[out] E     Matrix of eigenvectors in columns.
  /// \param[out] v     Converged eigenvalues.
  /// \param[in]  k     Number of requested eigenvalues.
  /// \param[in]  sigma Which eigenvalues of the spectrum, i.e., 'LM', 'LA', 'SM', 'SA',
  ///                   or numeric value sigma used for shift-and-invert mode.
  /// \param[in]  p     Number of Lanczos basis vectors.
  /// \param[in]  tol   Ritz estimate residual <= tol*norm(this).
  /// \param[in]  maxit Maximum number of iterations.
  /// \param[in]  v0    Starting vector. By default randomly generated.
  ///
  /// \returns Number of converged eigenvalues.
  ///
  /// \note Only implemented for real double precision sparse matrices.
  ///       Use only when MIRTK_Numerics_WITH_eigs is 1.
  int Eigenvectors(Matrix &E, Vector &v, int k, const char *sigma = "LM",
                   int p = 0, double tol = .0, int maxit = 0, Vector *v0 = NULL) const;

  // ---------------------------------------------------------------------------
  // I/O
#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)

  /// Create mxArray from sparse matrix
  ///
  /// This function is most efficient when the CCS layout is used.
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_MATLAB is 1.
  ///
  /// \returns Array created using mxCreateSparse.
  ///          Must be deleted by caller using mxDestroyArray.
  mxArray *MxArray() const;

#endif

  /// Read sparse matrix from MAT-file
  ///
  /// This function is most efficient when the CCS layout is used.
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_MATLAB is 1.
  bool ReadMAT(const char *, const char * = "A");

  /// Write sparse matrix to MAT-file
  ///
  /// This function is most efficient when the CCS layout is used.
  ///
  /// \note Use this function only when MIRTK_Numerics_WITH_MATLAB is 1.
  bool WriteMAT(const char *, const char * = "A") const;

  /// Write matrix to MATLAB script
  ///
  /// This function is most efficient when the CCS layout is used.
  void WriteMFile(const char *, const char * = "A") const;
};

////////////////////////////////////////////////////////////////////////////////
// Template definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::CheckEntries(Entries &entries) const
{
  const int n = (_Layout == CRS ? _Cols : _Rows);
  // Remove zero entries and check indices
  for (typename Entries::iterator i = entries.begin(); i != entries.end(); ++i) {
    if (i->second == TEntry(0)) {
      typename Entries::iterator j = i; --i;
      entries.erase(j);
    } else if (i->first < 0 || i->first >= n) {
      cerr << "GenericSparseMatrix::CheckEntries: Invalid "
                << (_Layout == CRS ? "column" : "row") << " index: "
                << i->first << endl;
      exit(1);
    }
  }
  // Sort non-zero entries
  sort(entries.begin(), entries.end());
  // Sum up duplicate entries
  for (typename Entries::iterator i = entries.begin(); i != entries.end(); ++i) {
    typename Entries::iterator j = i + 1;
    while (j != entries.end() && i->first == j->first) {
      i->second += j->second, ++j;
    }
    entries.erase(i + 1, i + std::distance(i, j));
  }
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TEntry> template <class TOtherEntry>
void GenericSparseMatrix<TEntry>
::CopyAttributes(const GenericSparseMatrix<TOtherEntry> &other)
{
  // Free this matrix
  Deallocate(_Row);
  Deallocate(_Col);
  Deallocate(_Data);
  Deallocate(_Index);
  // Copy attributes of other matrix
  _Layout = other._Layout;
  _Rows   = other._Rows;
  _Cols   = other._Cols;
  _NNZ    = other._NNZ;
  _Size   = other._NNZ; // disregard unused entries
  _MaxNumberOfUnusedEntries = other._MaxNumberOfUnusedEntries;
  // Copy non-zero entries of other matrix
  if (_Layout == CRS) {
    _Row  = Allocate<int>   (_Rows + 1);
    _Col  = Allocate<int>   (_Size);
    _Data = Allocate<TEntry>(_Size);
    for (int r = 0; r <= _Rows; ++r) _Row[r] = other._Row[r];
    for (int i = 0; i <= _NNZ;  ++i) {
      _Col [i] = other._Col [i];
      _Data[i] = other._Data[i];
    }
  } else {
    _Row  = Allocate<int>   (_Size);
    _Col  = Allocate<int>   (_Cols + 1);
    _Data = Allocate<TEntry>(_Size);
    for (int c = 0; c <= _Cols; ++c) _Col[c] = other._Col[c];
    for (int i = 0; i <= _NNZ;  ++i) {
      _Row [i] = other._Row [i];
      _Data[i] = other._Data[i];
    }
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry>::GenericSparseMatrix(StorageLayout layout)
:
  _Layout(layout),
  _Rows(0),
  _Cols(0),
  _NNZ(0),
  _Size(0),
  _MaxNumberOfUnusedEntries(100),
  _Row(NULL),
  _Col(NULL),
  _Data(NULL),
  _Index(NULL)
{
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry>
::GenericSparseMatrix(int m, int n, StorageLayout layout)
:
  _Layout(layout),
  _Rows(m),
  _Cols(n ? n : m),
  _NNZ(0),
  _Size(0),
  _MaxNumberOfUnusedEntries(100),
  _Row(NULL),
  _Col(NULL),
  _Data(NULL),
  _Index(NULL)
{
  if (_Layout == CRS) _Row = CAllocate<int>(_Rows + 1);
  else                _Col = CAllocate<int>(_Cols + 1);
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry>
::GenericSparseMatrix(int m, int n, int nnz, StorageLayout layout)
:
  _Layout(layout),
  _Rows(m),
  _Cols(n ? n : m),
  _NNZ(0),
  _Size(nnz),
  _MaxNumberOfUnusedEntries(100),
  _Row(NULL),
  _Col(NULL),
  _Data(NULL),
  _Index(NULL)
{
  if (_Layout == CRS) {
    _Row  = CAllocate<int>(_Rows + 1);
    _Col  =  Allocate<int>(nnz);
  } else {
    _Row =  Allocate<int>(nnz);
    _Col = CAllocate<int>(_Cols + 1);
  }
  _Data = Allocate<TEntry>(nnz);
}

// -----------------------------------------------------------------------------
template <class TEntry> template <class TOtherEntry>
GenericSparseMatrix<TEntry>
::GenericSparseMatrix(const GenericSparseMatrix<TOtherEntry> &other)
:
  _Row(NULL), _Col(NULL), _Data(NULL), _Index(NULL)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
template <class TEntry> template <class TOtherEntry>
GenericSparseMatrix<TEntry> &GenericSparseMatrix<TEntry>
::operator =(const GenericSparseMatrix<TOtherEntry> &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Initialize(int m, int n, int sz)
{
  if (n == 0) n = m;
  _Rows = m;
  _Cols = n;
  _NNZ  = 0;
  _Size = 0;
  Deallocate(_Row);
  Deallocate(_Col);
  Deallocate(_Data);
  Deallocate(_Index);
  if (_Rows > 0 && _Cols > 0) {
    if (_Layout == CRS) _Row = CAllocate<int>(_Rows + 1);
    else                _Col = CAllocate<int>(_Cols + 1);
  } else {
    _Rows = _Cols = 0;
  }
  if (sz > 0) {
    if (_Layout == CRS) _Col = Allocate<int>(sz);
    else                _Row = Allocate<int>(sz);
    _Data = Allocate<TEntry>(sz);
    _Size = sz;
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>
::Initialize(int m, int n, Array<Entries> &entries, bool as_is)
{
  // Free sparse matrix
  Deallocate(_Row);
  Deallocate(_Col);
  Deallocate(_Data);
  Deallocate(_Index);
  // Set attributes
  _Rows = m, _Cols = n;
  if (_Layout == CCS) swap(m, n);
  // Check/sort input entries
  if (!as_is) {
    for (int r = 0; r < m; ++r) CheckEntries(entries[r]);
  }
  _NNZ  = 0;
  for (int r = 0; r < m; ++r) {
    _NNZ += static_cast<int>(entries[r].size());
  }
  // Allocate memory
  _Row  = Allocate<int>(m + 1);
  _Col  = Allocate<int>(_NNZ);
  _Data = Allocate<EntryType>(_NNZ);
  _Size = _NNZ;
  // Copy non-zero entries
  int i = 0;
  for (int r = 0; r < m; ++r) {
    _Row[r] = i;
    typename Entries::const_iterator entry;
    for (entry = entries[r].begin(); entry != entries[r].end(); ++entry, ++i) {
      _Col [i] = entry->first;
      _Data[i] = entry->second;
    }
  }
  _Row[m] = _NNZ;
  // Swap arrays in case of CCS
  if (_Layout == CCS) swap(_Row, _Col);
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>
::Initialize(int m, int n, Entries entries[], bool as_is)
{
  // Free sparse matrix
  Deallocate(_Row);
  Deallocate(_Col);
  Deallocate(_Data);
  Deallocate(_Index);
  // Set attributes
  _Rows = m, _Cols = n;
  if (_Layout == CCS) swap(m, n);
  if (!as_is) {
    for (int r = 0; r < m; ++r) CheckEntries(entries[r]);
  }
  _NNZ = 0;
  for (int r = 0; r < m; ++r) {
    _NNZ += static_cast<int>(entries[r].size());
  }
  // Allocate memory
  _Row  = Allocate<int>(m + 1);
  _Col  = Allocate<int>(_NNZ);
  _Data = Allocate<EntryType>(_NNZ);
  _Size = _NNZ;
  // Copy non-zero entries
  int i = 0;
  for (int r = 0; r < m; ++r) {
    _Row[r] = i;
    typename Entries::const_iterator entry;
    for (entry = entries[r].begin(); entry != entries[r].end(); ++entry, ++i) {
      _Col [i] = entry->first;
      _Data[i] = entry->second;
    }
  }
  _Row[m] = _NNZ;
  // Swap arrays in case of CCS
  if (_Layout == CCS) swap(_Row, _Col);
}

// -----------------------------------------------------------------------------
#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)
template <class TEntry>
void GenericSparseMatrix<TEntry>::Initialize(const mxArray *pm)
{
  // Free sparse matrix
  Deallocate(_Row);
  Deallocate(_Col);
  Deallocate(_Data);
  Deallocate(_Index);
  // Set attributes
  _Rows = static_cast<int>(mxGetM(pm));
  _Cols = static_cast<int>(mxGetN(pm));
  // Copy matrix elements
  if (mxIsSparse(pm)) {
    const StorageLayout layout = _Layout;
    // Get pointers to sparse matrix data
    const mwIndex *ir = mxGetIr(pm);
    const mwIndex *jc = mxGetJc(pm);
    const double  *pr = mxGetPr(pm);
    // Allocate memory
    _Layout = CCS;
    _Size = _NNZ = jc[_Cols];
    _Col  = Allocate<int>(_Cols + 1);
    _Row  = Allocate<int>(_Size);
    _Data = Allocate<EntryType>(_Size);
    // Copy sparse matrix entries
    for (int c = 0; c <= _Cols; ++c) {
      _Col[c] = jc[c];
    }
    for (int i = 0; i < _NNZ; ++i) {
      _Row [i] = ir[i];
      _Data[i] = pr[i];
    }
    // Change layout to CRS if necessary
    Layout(layout);
  } else {
    cerr << "GenericSparseMatrix::Initialize: Initialization from non-sparse mxArray not implemented" << endl;
    exit(1);
  }
}
#endif // MIRTK_Numerics_WITH_MATLAB

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry>::~GenericSparseMatrix()
{
  Deallocate(_Row);
  Deallocate(_Col);
  Deallocate(_Data);
  Deallocate(_Index);
}

// =============================================================================
// Entries
// =============================================================================

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Layout(StorageLayout layout)
{
  if (layout != CRS && layout != CCS) {
    cerr << "GenericSparseMatrix::Layout: Unknown storage layout: " << layout << endl;
    exit(1);
  }
  if (layout != _Layout) {
    if (_NNZ == 0) swap(_Row, _Col);
    else           Transpose(true);
    _Layout = layout;
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>::GetRawData(int *&row, int *&col, TEntry *&data) const
{
  row  = _Row;
  col  = _Col;
  data = _Data;
  return _NNZ;
}

// -----------------------------------------------------------------------------
template <class TEntry>
TEntry *GenericSparseMatrix<TEntry>::RawPointer(int i)
{
  return &_Data[i];
}

// -----------------------------------------------------------------------------
template <class TEntry>
const TEntry *GenericSparseMatrix<TEntry>::RawPointer(int i) const
{
  return &_Data[i];
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Reserve(int sz)
{
  if (sz < _NNZ) sz = _NNZ;
  int    *&old_idx = (_Layout == CRS ? _Col : _Row);
  int    *new_idx  = Allocate<int>   (sz);
  TEntry *new_data = Allocate<TEntry>(sz);
  memcpy(new_idx,  old_idx, _NNZ * sizeof(int));
  memcpy(new_data, _Data,   _NNZ * sizeof(TEntry));
  Deallocate(old_idx);
  Deallocate(_Data);
  old_idx = new_idx;
  _Data   = new_data;
  _Size   = sz;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::RemoveZeros()
{
  // CCS is equivalent to CRS of transpose matrix
  int rows, *col, *row;
  if (_Layout == CCS) rows = _Cols, row = _Col, col = _Row;
  else                rows = _Rows, row = _Row, col = _Col;
  // Find zero entries and overwrite them by following non-zero entries
  int i = 0; // new position
  int j = 0; // old position
  for (int r = 0; r < rows; ++r) {
    row[r] = i;
    while (j < row[r+1]) {
      if (_Data[j] == .0) {
        --_NNZ;
      } else {
        if (i < j) col[i] = col[j], _Data[i] = _Data[j];
        ++i;
      }
      ++j;
    }
  }
  row[rows] = _NNZ;
  // Free memory if number of unused entries exceeds max buffer size
  if (_Size > _NNZ + _MaxNumberOfUnusedEntries) Reserve(_NNZ);
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>::RowNNZ(int r) const
{
  int nnz = 0;
  if (_Layout == CRS) {
    nnz = _Row[r+1] - _Row[r];
  } else {
    for (int c = 0; c < _Cols; ++c) {
      for (int i = _Col[c]; i != _Col[c+1]; ++i) {
        if (_Row[i] == r) ++nnz;
      }
    }
  }
  return nnz;
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>::ColNNZ(int c) const
{
  int nnz = 0;
  if (_Layout == CCS) {
    nnz = _Col[c+1] - _Col[c];
  } else {
    for (int r = 0; r < _Rows; ++r) {
      for (int i = _Row[r]; i != _Row[r+1]; ++i) {
        if (_Col[i] == c) ++nnz;
      }
    }
  }
  return nnz;
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>::SubNNZ(int r1, int c1, int r2, int c2) const
{
  int i, j, nnz = 0;
  if (_Layout == CRS) {
    for (int r = r1; r <= r2; ++r) {
      i = _Row[r];
      while (i != _Row[r+1] && _Col[i] <  c1) ++i;
      j = i;
      while (j != _Row[r+1] && _Col[j] <= c2) ++j;
      nnz += j - i;
    }
  } else {
    for (int c = c1; c <= c2; ++c) {
      i = _Col[c];
      while (i != _Col[c+1] && _Row[i] <  r1) ++i;
      j = i;
      while (j != _Col[c+1] && _Row[j] <= r2) ++j;
      nnz += j - i;
    }
  }
  return nnz;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Index()
{
  if (!_Index) {
    MIRTK_START_TIMING();
    int i, j;
    if (_Layout == CRS) {
      _Index = Allocate<Array<int> >(_Cols);
      Array<int> count(_Cols, 0);
      for (i = 0; i < _NNZ;  ++i) {
        ++count[_Col[i]];
      }
      for (j = 0; j < _Cols; ++j) {
        _Index[j].resize(count[j]);
        count [j] = 0;
      }
      for (i = 0; i < _NNZ;  ++i) {
        j = _Col[i];
        _Index[j][count[j]++] = i;
      }
    } else {
      _Index = Allocate<Array<int> >(_Rows);
      Array<int> count(_Rows, 0);
      for (i = 0; i < _NNZ;  ++i) {
        ++count[_Row[i]];
      }
      for (j = 0; j < _Rows; ++j) {
        _Index[j].resize(count[j]);
        count [j] = 0;
      }
      for (i = 0; i < _NNZ;  ++i) {
        j = _Row[i];
        _Index[j][count[j]++] = i;
      }
    }
    MIRTK_DEBUG_TIMING(10, "building sparse matrix index");
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::ClearIndex()
{
  Deallocate(_Index);
}

// -----------------------------------------------------------------------------
template <class TEntry>
TEntry &GenericSparseMatrix<TEntry>::operator ()(int r, int c)
{
  if (c < 0) c = r;
  // CCS is equivalent to CRS of transpose matrix
  int rows, *col, *row;
  if (_Layout == CCS) {
    swap(r, c);
    row  = _Col;
    col  = _Row;
    rows = _Cols;
  } else {
    row  = _Row;
    col  = _Col;
    rows = _Rows;
  }
  // Find existing non-zero entry or insert position
  int i;
  for (i = row[r]; i != row[r+1]; ++i) {
    if (col[i] >  c) break;
    if (col[i] == c) return _Data[i];
  }
  // Invalidate precomputed _Index
  Deallocate(_Index);
  // Reserve memory for new number of non-zero elements
  if (_NNZ == _Size) {
    Reserve(_Size + _MaxNumberOfUnusedEntries);
    if (_Layout == CCS) row = _Col, col = _Row;
    else                row = _Row, col = _Col;
  }
  // Insert new non-zero element at i-th position
  for (int j = _NNZ; j > i; --j) {
    col  [j] = col  [j-1];
    _Data[j] = _Data[j-1];
  }
  col  [i] = c;
  _Data[i] = TEntry(0);
  for (int j = r + 1; j < rows; ++j) ++row[j];
  row[rows] = ++_NNZ;
  return _Data[i];
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Put(int r, int c, TEntry v)
{
  // Remove non-zero entry if set to zero
  if (v == TEntry(0)) {
    // CCS is equivalent to CRS of transpose matrix
    int rows, *col, *row;
    if (_Layout == CCS) {
      swap(r, c);
      row  = _Col;
      col  = _Row;
      rows = _Cols;
    } else {
      row  = _Row;
      col  = _Col;
      rows = _Rows;
    }
    // Remove existing non-zero entry
    for (int i = row[r]; i != row[r+1]; ++i) {
      if (col[i] <  c) continue;
      if (col[i] == c) {
        Deallocate(_Index);
        for (int j = i + 1; j < _NNZ; ++j) {
          col  [j-1] = col  [j];
          _Data[j-1] = _Data[j];
        }
        for (int j = r + 1; j < rows; ++j) --row[j];
        row[rows] = (--_NNZ);
        if (_Size > _NNZ + _MaxNumberOfUnusedEntries) {
          Reserve(_NNZ + _MaxNumberOfUnusedEntries);
        }
      }
      break;
    }
  // Otherwise, set/add non-zero entry
  } else {
    this->operator ()(r, c) = v;
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
TEntry GenericSparseMatrix<TEntry>::Get(int r, int c) const
{
  if (c < 0) c = r;
  if (_Layout == CRS) {
    for (int i = _Row[r]; i != _Row[r+1]; ++i) {
      if (_Col[i] == c) return _Data[i];
      if (_Col[i] >  c) break;
    }
  } else {
    for (int i = _Col[c]; i != _Col[c+1]; ++i) {
      if (_Row[i] == r) return _Data[i];
      if (_Row[i] >  r) break;
    }
  }
  return TEntry(0);
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Row(int r, Entries &row, bool as_is)
{
  if (!as_is) {
    // Remove zero entries
    for (typename Entries::iterator i = row.begin(); i != row.end(); ++i) {
      if (i->second == TEntry(0)) {
        typename Entries::iterator j = i; --i;
        row.erase(j);
      }
    }
    // Sort non-zero entries
    sort(row.begin(), row.end());
  }
  // Set non-zero entries of r-th row
  if (_Layout == CRS) {
    // Number of additional (positive) or fewer (negative) row entries
    int dnnz = static_cast<int>(row.size()) - (_Row[r+1] - _Row[r]);
    // Allocate more memory if needed
    if (_NNZ + dnnz > _Size) Reserve(_NNZ + dnnz);
    // Increase/Decrease counter of non-zero entries
    _NNZ += dnnz;
    // Move non-zero entries of rows > r to leave just enough space for
    // the new non-zero entries of the r-th row
    if (dnnz != 0) {
      for (int j = r + 1; j < _Rows; ++j) _Row[j] += dnnz;
    }
    _Row[_Rows] = _NNZ;
    if (dnnz < 0) {
      for (int i = _Row[r+1], j = i - dnnz; i < _NNZ; ++i, ++j) {
        _Col [i] = _Col [j];
        _Data[i] = _Data[j];
      }
    } else if (dnnz > 0) {
      for (int i = _NNZ - 1, j = i - dnnz; i >= _Row[r+1]; --i, --j) {
        _Col [i] = _Col [j];
        _Data[i] = _Data[j];
      }
    }
    // Set non-zero entries of r-th row
    bool index_invalid = (dnnz != 0);
    typename Entries::const_iterator entry = row.begin();
    for (int i = _Row[r]; i != _Row[r+1]; ++i, ++entry) {
      if (_Col[i] != entry->first) {
         _Col[i] = entry->first;
        index_invalid = true;
      }
      _Data[i] = entry->second;
    }
    // Delete invalid precomputed index
    if (index_invalid) Deallocate(_Index);
    // Free memory if number of unused entries exceeds max buffer size
    if (_NNZ < _Size - _MaxNumberOfUnusedEntries) Reserve(_NNZ);
  } else {
    // Preallocate more memory if needed
    if (_NNZ + static_cast<int>(row.size()) > _Size) {
      Reserve(_NNZ + static_cast<int>(row.size()));
    }
    // Insert new non-zero row entries
    typename Entries::const_iterator entry = row.begin();
    for (; entry != row.end(); ++entry) {
      this->operator ()(r, entry->first) = entry->second;
    }
    // Free memory if number of unused entries exceeds max buffer size
    if (_Size > _NNZ + _MaxNumberOfUnusedEntries) Reserve(_NNZ);
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::GetRow(int r, Entries &row) const
{
  if (_Layout == CRS) {
    row.resize(_Row[r+1] - _Row[r]);
    for (int i = _Row[r], c = 0; i != _Row[r+1]; ++i, ++c) {
      row[c] = MakePair(_Col[i], _Data[i]);
    }
  } else {
    if (_Index) {
      int i, c = 0;
      row.resize(_Index[r].size());
      for (size_t j = 0; j < _Index[r].size(); ++j) {
        i = _Index[r][j];
        while (i >= _Col[c+1]) ++c;
        row[j] = MakePair(c, _Data[i]);
      }
    } else {
      row.clear();
      for (int c = 0; c < _Cols; ++c) {
        for (int i = _Col[c]; i != _Col[c+1]; ++i) {
          if (_Row[i] <  r) continue;
          if (_Row[i] == r) row.push_back(MakePair(c, _Data[i]));
          break;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
typename GenericSparseMatrix<TEntry>::Entries
GenericSparseMatrix<TEntry>::Row(int r) const
{
  Entries row;
  GetRow(r, row);
  return row;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Col(int c, Entries &col, bool as_is)
{
  if (!as_is) {
    // Remove zero entries
    for (typename Entries::iterator i = col.begin(); i != col.end(); ++i) {
      if (i->second == TEntry(0)) {
        typename Entries::iterator j = i; --i;
        col.erase(j);
      }
    }
    // Sort non-zero entries
    sort(col.begin(), col.end());
  }
  // Set non-zero entries of c-th column
  if (_Layout == CCS) {
    // Number of additional (positive) or fewer (negative) column entries
    int dnnz = static_cast<int>(col.size()) - (_Col[c+1] - _Col[c]);
    // Allocate more memory if needed
    if (_NNZ + dnnz > _Size) Reserve(_NNZ + dnnz);
    // Increase/Decrease counter of non-zero entries
    _NNZ += dnnz;
    // Move non-zero entries of columns > c to leave just enough space for
    // the new non-zero entries of the c-th column
    if (dnnz != 0) {
      for (int j = c + 1; j < _Cols; ++j) _Col[j] += dnnz;
    }
    _Col[_Cols] = _NNZ;
    if (dnnz < 0) {
      for (int i = _Col[c+1], j = i - dnnz; i < _NNZ; ++i, ++j) {
        _Row [i] = _Row [j];
        _Data[i] = _Data[j];
      }
    } else if (dnnz > 0) {
      for (int i = _NNZ - 1, j = i - dnnz; i >= _Col[c+1]; --i, --j) {
        _Row [i] = _Row [j];
        _Data[i] = _Data[j];
      }
    }
    // Set non-zero entries of r-th row
    bool index_invalid = (dnnz != 0);
    typename Entries::const_iterator entry = col.begin();
    for (int i = _Col[c]; i != _Col[c+1]; ++i, ++entry) {
      if (_Row[i] != entry->first) {
        _Row[i] = entry->first;
        index_invalid = true;
      }
      _Data[i] = entry->second;
    }
    // Delete invalid precomputed _Index
    if (index_invalid) Deallocate(_Index);
    // Free memory if number of unused entries exceeds max buffer size
    if (_NNZ < _Size - _MaxNumberOfUnusedEntries) Reserve(_NNZ);
  } else {
    // Preallocate more memory if needed
    if (_NNZ + static_cast<int>(col.size()) > _Size) {
      Reserve(_NNZ + static_cast<int>(col.size()));
    }
    // Insert new non-zero row entries
    typename Entries::const_iterator entry = col.begin();
    for (; entry != col.end(); ++entry) {
      this->operator ()(entry->first, c) = entry->second;
    }
    // Free memory if number of unused entries exceeds max buffer size
    if (_Size > _NNZ + _MaxNumberOfUnusedEntries) Reserve(_NNZ);
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Column(int c, Entries &col, bool as_is)
{
  Col(c, col, as_is);
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::GetCol(int c, Entries &col) const
{
  if (_Layout == CRS) {
    if (_Index) {
      int i, r = 0;
      col.resize(_Index[c].size());
      for (size_t j = 0; j < _Index[c].size(); ++j) {
        i = _Index[c][j];
        while (i >= _Row[r+1]) ++r;
        col[j] = MakePair(r, _Data[i]);
      }
    } else {
      col.clear();
      for (int r = 0; r < _Rows; ++r) {
        for (int i = _Row[r]; i != _Row[r+1]; ++i) {
          if (_Col[i] <  c) continue;
          if (_Col[i] == c) col.push_back(MakePair(r, _Data[i]));
          break;
        }
      }
    }
  } else {
    col.resize(_Col[c+1] - _Col[c]);
    for (int i = _Col[c], r = 0; i != _Col[c+1]; ++i, ++r) {
      col[r] = MakePair(_Row[i], _Data[i]);
    }
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::GetColumn(int c, Entries &col) const
{
  GetCol(c, col);
}

// -----------------------------------------------------------------------------
template <class TEntry>
typename GenericSparseMatrix<TEntry>::Entries
GenericSparseMatrix<TEntry>::Col(int c) const
{
  Entries col;
  GetCol(c, col);
  return col;
}

// -----------------------------------------------------------------------------
template <class TEntry>
typename GenericSparseMatrix<TEntry>::Entries
GenericSparseMatrix<TEntry>::Column(int c) const
{
  return Col(c);
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::GetDiag(Vector &d) const
{
  if (_Rows != _Cols) {
    cerr << "GenericSparseMatrix::GetDiag: Matrix must be square" << endl;
    exit(1);
  }
  d.Initialize(_Rows);
  for (int i = 0; i < _Rows; ++i) d(i) = this->Get(i, i);
}

// -----------------------------------------------------------------------------
template <class TEntry>
Vector GenericSparseMatrix<TEntry>::Diag() const
{
  Vector d;
  GetDiag(d);
  return d;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Diag(const Vector &d)
{
  if (_Rows != _Cols) {
    cerr << "GenericSparseMatrix::Diag: Matrix must be square" << endl;
    exit(1);
  }
  if (d.Rows() != _Rows) {
    cerr << "GenericSparseMatrix::Diag: Invalid vector size" << endl;
    exit(1);
  }
  // CCS is equivalent to CRS of transpose matrix
  int rows, *col, *row;
  if (_Layout == CCS) rows = _Cols, row = _Col, col = _Row;
  else                rows = _Rows, row = _Row, col = _Col;
  // Count number of missing diagonal entries
  int i, j, r, dnnz = 0;
  for (r = 0; r < rows; ++r) {
    if (d(r) != .0) {
      i = row[r];
      while (i < row[r+1] && col[i] < r) ++i;
      if (i == row[r+1] || col[i] != r) ++dnnz;
    }
  }
  // Reserve memory for missing diagonal entries
  if (_NNZ + dnnz > _Size) {
    Reserve(_NNZ + dnnz);
    if (_Layout == CCS) row = _Col, col = _Row;
    else                row = _Row, col = _Col;
  }
  // Delete invalid precomputed _Index
  if (dnnz != 0) Deallocate(_Index);
  // Insert/set diagonal entries
  _NNZ += dnnz;
  row[rows] = _NNZ;    // new nnz
  i = _NNZ - 1 - dnnz; // old position
  j = _NNZ - 1;        // new position
  r = rows - 1;        // iterate entries in reverse order
  while (r >= 0 && i < j) {
    // Move columns above diagonal to make space for new entries
    for (; i >= row[r] && col[i] > r; --i, --j) {
      col[j] = col[i], _Data[j] = _Data[i];
    }
    // Set/insert diagonal entry
    if (i >= row[r] && col[i] == r) {
      col[j] = r, _Data[j] = static_cast<TEntry>(d(r)), --i, --j;
    } else if (d(r) != .0) {
      col[j] = r, _Data[j] = static_cast<TEntry>(d(r)), --j;
    }
    // Move columns below diagonal to make space for new entries
    if (i < j) {
      for (; i >= row[r]; --i, --j) {
        col[j] = col[i], _Data[j] = _Data[i];
      }
    } else {
      i = j = row[r] - 1;
    }
    // Update row start index
    row[r] = j + 1;
    // Next row above this one
    --r;
  }
  // ...just copy remaining entries and modify the diagonal value in each row
  mirtkAssert(i == j, "no more missing diagonal entries");
  while (r >= 0) {
    if (d(r) != .0) {
      // Leave columns above diagonal untouched
      while (i >= row[r] && col[i] > r) --i;
      // Modify diagonal entry of this row
      mirtkAssert(i >= row[r] && col[i] == r, "diagonal entry must exist");
      _Data[i] = static_cast<TEntry>(d(r));
      // Leave columns below diagonal untouched
    }
    i = row[r] - 1;
    // Next row above this one
    --r;
  }
  // Discard previous diagonal entries which are now zero
  RemoveZeros();
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Diag(TEntry d)
{
  Diag(Vector(_Rows, static_cast<double>(d)));
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry>
GenericSparseMatrix<TEntry>::Sub(int r1, int c1, int r2, int c2) const
{
  GenericSparseMatrix<TEntry> sub(r2 - r1 + 1, c2 - c1 + 1, _Layout);
  cerr << "GenericSparseMatrix::Sub (getter): Not implemented" << endl;
  exit(1);
  return sub;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>
::Sub(int r1, int c1, const GenericSparseMatrix &sub)
{
  MIRTK_START_TIMING();
  int r2 = r1 + sub.Rows() - 1;
  int c2 = c1 + sub.Cols() - 1;
  if (r2 >= _Rows || c2 >= _Cols) {
    cerr << "GenericSparseMatrix::Sub (setter): " << sub.Rows() << " x " << sub.Cols()
              << " submatrix starting at (" << r1 << ", " << c1 << ") exceeds matrix index range [0, "
              << _Rows-1 << "] x [0, " << _Cols-1 << "]" << endl;
    exit(1);
  }
  if (_Layout != sub.Layout()) {
    // TODO: Accept also submatrix in different storage layout
    cerr << "GenericSparseMatrix::Sub (setter): Submatrix must have same layout as this matrix" << endl;
    exit(1);
  }
  // Determine new number of non-zero entries
  const int dif_nnz = sub._NNZ - this->SubNNZ(r1, c1, r2, c2);
  const int new_nnz = _NNZ + dif_nnz;
  // Allocate more memory if needed
  if (new_nnz > _Size) Reserve(new_nnz);
  // CCS is equivalent to CRS of transpose matrix
  int rows, *col, *row;
  const int *sub_col, *sub_row;
  if (_Layout == CCS) {
    swap(r1, c1);
    swap(r2, c2);
    row     = _Col;
    col     = _Row;
    rows    = _Cols;
    sub_row = sub._Col;
    sub_col = sub._Row;
  } else {
    row     = _Row;
    col     = _Col;
    rows    = _Rows;
    sub_row = sub._Row;
    sub_col = sub._Col;
  }
  // New position of new or existing i-th non-zero entry
  int i, k, k1, k2, j = new_nnz - 1;
  // Move non-zero entries below submatrix to the end
  for (i = _NNZ - 1; i >= row[r2 + 1]; --i, --j) {
    col  [j] = col  [i];
    _Data[j] = _Data[i];
  }
  for (int r = r2 + 1; r < rows; ++r) row[r] += dif_nnz;
  row[rows] = _NNZ = new_nnz;
  // Insert submatrix rows in reverse order
  for (int r = r2; r >= r1; --r) {
    // Existing non-zero entries with column index c > c2
    for (i = row[r+1] - 1; i >= row[r] && col[i] > c2; --i, --j) {
      col  [j] = col  [i];
      _Data[j] = _Data[i];
    }
    // Skip existing non-zero entries with column index c1 <= c <= c2
    while (i >= row[r] && col[i] >= c1) --i;
    // ...and replace them by new submatrix entries instead
    k1 = sub_row[r - r1];
    k2 = sub_row[r - r1 + 1] - 1;
    for (k = k2; k >= k1; --k, --j) {
      col  [j] = c1 + sub_col[k];
      _Data[j] = sub._Data[k];
    }
    // Existing non-zero entries with column index c < c1
    for (; i >= row[r]; --i, --j) {
      col  [j] = col  [i];
      _Data[j] = _Data[i];
    }
    // Row start index
    row[r] = j + 1;
  }
  MIRTK_DEBUG_TIMING(7, "GenericSparseMatrix::Sub (setter)");
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Clear()
{
  ClearIndex();
  if (_Layout == CRS) {
    if (_Row) memset(_Row, 0, (_Rows + 1) * sizeof(int));
    Deallocate(_Col);
  } else {
    Deallocate(_Row);
    if (_Col) memset(_Col, 0, (_Cols + 1) * sizeof(int));
  }
  Deallocate(_Data);
  _Rows = _Cols = _Size = _NNZ = 0;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Zero()
{
  _NNZ = 0;
  if (_Layout == CRS) {
    if (_Row) memset(_Row, 0, (_Rows + 1) * sizeof(int));
  } else {
    if (_Col) memset(_Col, 0, (_Cols + 1) * sizeof(int));
  }
  Deallocate(_Index);
}

// =============================================================================
// Common operations
// =============================================================================

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::Transpose(bool keep_layout)
{
  if (keep_layout) {
    Index(); // makes Col/Row retrieval efficient
    Array<Entries> entries;
    if (_Layout == CRS) {
      entries.resize(_Cols);
      for (int c = 0; c < _Cols; ++c) entries[c] = Col(c);
      _Layout = CCS; // *after* all columns are retrieved
    } else {
      entries.resize(_Rows);
      for (int r = 0; r < _Rows; ++r) entries[r] = Row(r);
      _Layout = CRS; // *after* all rows are retrieved
    }
    Initialize(_Cols, _Rows, entries);
  } else {
    if (_Layout == CRS) _Layout = CCS;
    else                _Layout = CRS;
    swap(_Rows, _Cols);
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
TEntry GenericSparseMatrix<TEntry>::RowSum(int r) const
{
  TEntry s = TEntry(0);
  if (_Layout == CRS) {
    for (int i = _Row[r]; i != _Row[r+1]; ++i) s += _Data[i];
  } else {
    if (_Index) {
      for (size_t c = 0; c < _Index[r].size(); ++c) {
        s += _Data[_Index[r][c]];
      }
    } else {
      for (int c = 0; c < _Cols; ++c) {
        for (int i = _Col[c]; i != _Col[c+1]; ++i) {
          if (_Row[i] <  r) continue;
          if (_Row[i] == r) s += _Data[i];
          break;
        }
      }
    }
  }
  return s;
}

// -----------------------------------------------------------------------------
template <class TEntry>
TEntry GenericSparseMatrix<TEntry>::ColSum(int c) const
{
  TEntry s = TEntry(0);
  if (_Layout == CRS) {
    if (_Index) {
      for (size_t r = 0; r < _Index[c].size(); ++r) {
        s += _Data[_Index[c][r]];
      }
    } else {
      for (int r = 0; r < _Rows; ++r) {
        for (int i = _Row[r]; i != _Row[r+1]; ++i) {
          if (_Col[i] <  c) continue;
          if (_Col[i] == c) s += _Data[i];
          break;
        }
      }
    }
  } else {
    for (int i = _Col[c]; i != _Col[c+1]; ++i) s += _Data[i];
  }
  return s;
}

// -----------------------------------------------------------------------------
template <class TEntry>
TEntry GenericSparseMatrix<TEntry>::ColumnSum(int c) const
{
  return ColSum(c);
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::MultAv(TEntry v[], TEntry w[]) const
{
  const TEntry zero(0);
  if (_Layout == CRS) {
    for (int r = 0; r < _Rows; ++r) {
      w[r] = zero;
      for (int i = _Row[r]; i != _Row[r+1]; ++i) w[r] += _Data[i] * v[_Col[i]];
    }
  } else {
    if (_Index) {
      int i;
      for (int r = 0; r < _Rows; ++r) {
        w[r] = zero;
        for (size_t c = 0; c < _Index[r].size(); ++c) {
          i = _Index[r][c];
          w[r] += _Data[i] * v[_Col[i]];
        }
      }
    } else {
      for (int r = 0; r < _Rows; ++r) w[r] = zero;
      for (int c = 0; c < _Cols; ++c) {
        for (int i = _Col[c]; i != _Col[c+1]; ++i) {
          w[_Row[c]] += _Data[i] * v[c];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> &GenericSparseMatrix<TEntry>::operator *=(TEntry s)
{
  for (int i = 0; i < _NNZ; ++i) _Data[i] *= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> &GenericSparseMatrix<TEntry>::operator /=(TEntry s)
{
  for (int i = 0; i < _NNZ; ++i) _Data[i] /= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> GenericSparseMatrix<TEntry>::operator *(TEntry s) const
{
  GenericSparseMatrix<TEntry> m(*this);
  return m *= s;
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> GenericSparseMatrix<TEntry>::operator /(TEntry s) const
{
  GenericSparseMatrix<TEntry> m(*this);
  return m /= s;
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> &GenericSparseMatrix<TEntry>::ScaleRow(int r, TEntry s)
{
  if (_Layout == CRS) {
    for (int i = _Row[r]; i != _Row[r+1]; ++i) _Data[i] *= s;
  } else {
    if (_Index) {
      for (size_t c = 0; c < _Index[r].size(); ++c) {
        _Data[_Index[r][c]] *= s;
      }
    } else {
      for (int c = 0; c < _Cols; ++c) {
        for (int i = _Col[c]; i != _Col[c+1]; ++i) {
          if (_Row[i] <  r) continue;
          if (_Row[i] == r) _Data[i] *= s;
          break;
        }
      }
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> &GenericSparseMatrix<TEntry>::ScaleColumn(int c, TEntry s)
{
  if (_Layout == CRS) {
    if (_Index) {
      for (size_t r = 0; r < _Index[c].size(); ++r) {
        _Data[_Index[c][r]] *= s;
      }
    } else {
      for (int r = 0; r < _Rows; ++r) {
        for (int i = _Row[r]; i != _Row[r+1]; ++i) {
          if (_Col[i] <  c) continue;
          if (_Col[i] == c) _Data[i] *= s;
          break;
        }
      }
    }
  } else {
    for (int i = _Col[c]; i != _Col[c+1]; ++i) _Data[i] *= s;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class TEntry>
GenericSparseMatrix<TEntry> &GenericSparseMatrix<TEntry>::ScaleCol(int c, TEntry s)
{
  return ScaleColumn(c, s);
}

// -----------------------------------------------------------------------------
template <class TEntry>
bool GenericSparseMatrix<TEntry>::IsSymmetric() const
{
  if (_Rows != _Cols) return false;
  // CCS is equivalent to CRS of transpose matrix
  int rows, *col, *row;
  if (_Layout == CCS) rows = _Cols, row = _Col, col = _Row;
  else                rows = _Rows, row = _Row, col = _Col;
  // Compare transposed upper triangular submatrix to lower triangular submatrix
  int i, iend, j, jend, r, c;
  for (r = 0; r < rows-1; ++r) {
    i    = row[r];
    iend = row[r+1];
    while (i < iend && col[i] < r) ++i;
    while (i < iend) {
      c    = col[i];
      j    = row[c];
      jend = row[c+1];
      while (j < jend && col[j] < r) ++j;
      if (j == jend || col[j] != r || !fequal(_Data[i], _Data[j])) {
        return false;
      }
      ++i;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>::MakeSymmetric(bool extent)
{
  if (_Rows != _Cols) {
    cerr << "GenericSparseMatrix::MakeSymmetric: Matrix must be square" << endl;
    exit(1);
  }
  // CCS is equivalent to CRS of transpose matrix
  int rows, *col, *row;
  if (_Layout == CCS) rows = _Cols, row = _Col, col = _Row;
  else                rows = _Rows, row = _Row, col = _Col;
  // Count number of missing transposed entries
  int i, j, k, r, c, dnnz = 0;
  for (r = 0; r < rows; ++r) {
    for (i = row[r]; i < row[r+1]; ++i) {
      c = col[i];
      j = row[c];
      while (j < row[c+1] && col[j] < r) ++j;
      if (j < row[c+1] && col[j] != r) ++dnnz;
    }
  }
  // Reserve memory for missing transposed entries
  if (_NNZ + dnnz > _Size) {
    Reserve(_NNZ + dnnz);
    if (_Layout == CCS) row = _Col, col = _Row;
    else                row = _Row, col = _Col;
  }
  // Add values above main diagonal to those below it and divide by two
  // Note that this way new elements are only added to rows below the current
  // one and therefore the current iteration bounds are not affected by it.
  for (r = 0; r < rows-1; ++r) {
    i = row[r];
    while (i < row[r+1] && col[i] < r) ++i;
    while (i < row[r+1]) {
      c = col[i];
      j = row[c];
      while (j < row[c+1] && col[j] < r) ++j;
      if (j < row[c+1] && col[j] == r) {
        _Data[i] = _Data[j] = (_Data[i] + _Data[j]) / 2;
      } else {
        if (!extent) _Data[i] = _Data[i] / 2;
        for (k = _NNZ; k > j; --k) {
          col  [k] = col  [k-1];
          _Data[k] = _Data[k-1];
        }
        col  [j] = r;
        _Data[j] = _Data[i];
        for (k = c + 1; k < rows; ++k) ++row[k];
        row[rows] = ++_NNZ;
      }
      ++i;
    }
  }
}

// =============================================================================
// I/O
// =============================================================================

#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)

// -----------------------------------------------------------------------------
template <class TEntry>
mxArray *GenericSparseMatrix<TEntry>::MxArray() const
{
  Matlab::Initialize();
  mxArray *sa = mxCreateSparse(_Rows, _Cols, _NNZ, mxREAL);
  mwIndex *ir = mxGetIr(sa);
  mwIndex *jc = mxGetJc(sa);
  double  *pr = mxGetPr(sa);
  if (_Layout == CCS) {
    for (int c = 0; c < _Cols; ++c) {
      jc[c] = _Col[c];
    }
    for (int i = 0; i < _NNZ; ++i) {
      ir[i] = _Row [i];
      pr[i] = _Data[i];
    }
  } else {
    Entries col;
    int i = 0;
    for (int c = 0; c < _Cols; ++c) {
      jc[c] = i;
      GetCol(c, col);
      typename Entries::const_iterator entry;
      for (entry = col.begin(); entry != col.end(); ++entry, ++i) {
        ir[i] = entry->first;
        pr[i] = entry->second;
      }
    }
  }
  jc[_Cols] = _NNZ;
  return sa;
}

#endif // MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)

// -----------------------------------------------------------------------------
template <class TEntry>
bool GenericSparseMatrix<TEntry>::ReadMAT(const char *fname, const char *varname)
{
#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)
  Matlab::Initialize();
  MATFile *fp = matOpen(fname, "r");
  if (fp == NULL) return false;
  mxArray *pm = matGetVariable(fp, varname);
  if (pm == NULL) {
    matClose(fp);
    return false;
  }
  Initialize(pm);
  mxDestroyArray(pm);
  return (matClose(fp) == 0);
#else
  cerr << "GenericSparseMatrix::ReadMAT: Must be compiled WITH_MATLAB enabled" << endl;
  return false;
#endif
}

// -----------------------------------------------------------------------------
template <class TEntry>
bool GenericSparseMatrix<TEntry>::WriteMAT(const char *fname, const char *varname) const
{
#if MIRTK_Numerics_WITH_MATLAB && defined(HAVE_MATLAB)
  Matlab::Initialize();
  MATFile *fp = matOpen(fname, "w");
  if (fp == NULL) return false;
  mxArray *m = MxArray();
  if (matPutVariable(fp, varname, m) != 0) {
    mxDestroyArray(m);
    matClose(fp);
    return false;
  }
  mxDestroyArray(m);
  return (matClose(fp) == 0);
#else
  cerr << "GenericSparseMatrix::WriteMAT: Must be compiled WITH_MATLAB enabled" << endl;
  return false;
#endif
}

// -----------------------------------------------------------------------------
template <class TEntry>
void GenericSparseMatrix<TEntry>
::WriteMFile(const char *fname, const char *varname) const
{
  ofstream of(fname);
  // Row indices
  of << varname << "_r = [";
  if (_Layout == CRS) {
    for (int r = 1; r <= _Rows; ++r) {
      for (int i = _Row[r-1]; i != _Row[r]; ++i) {
        of << r << " ";
      }
    }
  } else {
    for (int i = 0; i < _NNZ; ++i) {
      of << (_Row[i] + 1) << " ";
    }
  }
  of << "];" << endl;
  // Column indices
  of << varname << "_c = [";
  if (_Layout == CCS) {
    for (int c = 1; c <= _Cols; ++c) {
      for (int i = _Col[c-1]; i != _Col[c]; ++i) {
        of << c << " ";
      }
    }
  } else {
    for (int i = 0; i < _NNZ; ++i) {
      of << (_Col[i] + 1) << " ";
    }
  }
  of << "];" << endl;
  // Non-zero values
  of << varname << "_v = [";
  for (int i = 0; i < _NNZ; ++i) {
    of << _Data[i] << " ";
  }
  of << "];" << endl;
  of << varname << " = sparse(" << varname << "_r, " << varname << "_c, " << varname << "_v, ";
  of << _Rows << ", " << _Cols << ", " << _NNZ << ");" << endl;
}

////////////////////////////////////////////////////////////////////////////////
// Sparse matrix types
////////////////////////////////////////////////////////////////////////////////

/// Typically used template instantiations
typedef GenericSparseMatrix<float>  SparseFloatMatrix;
typedef GenericSparseMatrix<double> SparseDoubleMatrix;

/// Sparse matrix with default entry value type
#if MIRTK_USE_FLOAT_BY_DEFAULT
typedef SparseFloatMatrix SparseMatrix;
#else
typedef SparseDoubleMatrix SparseMatrix;
#endif


} // namespace mirtk

#endif // MIRTK_SparseMatrix_H
