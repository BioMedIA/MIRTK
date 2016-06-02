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

#include "mirtk/SparseMatrix.h"

#include "mirtk/NumericsConfig.h"
#if MIRTK_Numerics_WITH_eigs
#  include "mirtk/Arpack.h"
#  include "mirtk/Umfpack.h"
#  include "boost/random/mersenne_twister.hpp"
#  include "boost/random/uniform_01.hpp"
#endif // MIRTK_Numerics_WITH_eigs


namespace mirtk {


// =============================================================================
// Eigen decomposition
// =============================================================================

// -----------------------------------------------------------------------------
int eigs(const GenericSparseMatrix<double> &A, Matrix *E, Vector &v,
         int k, const char *eigs_sigma, int p, double tol, int maxit, Vector *v0)
{
  int nconv = 0; // Number of converged eigenvalues

#if MIRTK_Numerics_WITH_eigs
  if (A.Cols() != A.Rows()) {
    cerr << "eigs: Matrix must be square" << endl;
    exit(1);
  }
  int n = A.Rows();

  // Check input and derive (default) ARPACK parameters
  bool issym = A.IsSymmetric();
  if (p     <=  0) p     = min(max(2 * k + (issym ? 0 : 1), 20), n);
  if (tol   <= .0) tol   = numeric_limits<double>::epsilon();
  if (maxit <=  0) maxit = max(300.0, ceil(2 * n / max(p, 1)));
  double *resid = Allocate<double>(n);

  if (v0 && v0->Rows() != 0) {
    if (v0->Rows() != n) {
      cerr << "eigs: Initial vector v0 must have " << n << " rows" << endl;
      exit(1);
    }
    for (int i = 0; i < n; ++i) resid[i] = v0->Get(i);
  } else {
    boost::mt19937            gen;
    boost::uniform_01<double> dist;
    for (int i = 0; i < n; ++i) resid[i] = dist(gen);
  }

  char   which[3] = { "LM" };
  double sigma    = .0;
  int    mode     = 1;

  if (strlen(eigs_sigma) == 2 && (!isdigit(eigs_sigma[0]) || !isdigit(eigs_sigma[1]))) {
    which[0] = toupper(eigs_sigma[0]);
    which[1] = toupper(eigs_sigma[1]);
    if (strcmp(which, "SM") == 0) {
      which[0] = 'L';
      mode     = 3;
    }
    if (strcmp(which, "LM") != 0 && strcmp(which, "SM") != 0 &&
        strcmp(which, "LA") != 0 && strcmp(which, "SA") != 0) {
      cerr << "eigs: Invalid sigma string: " << sigma << endl;
      exit(1);
    }
  } else {
    if (!FromString(eigs_sigma, sigma)) {
      cerr << "eigs: Invalid sigma string or value: " << sigma << endl;
      exit(1);
    }
    mode = 3;
  }

  // LU factorization of (A - sigma I) in shift-and-invert mode
  GenericSparseMatrix<double> *AsI = NULL;
  void                        *Alu = NULL;

  int    *Ai = NULL, *Ap = NULL, status;
  double *Ax = NULL;

  if (mode == 3) {
    void *symbolic = NULL;
    if (sigma == .0) {
      if (A.Layout() == GenericSparseMatrix<double>::CCS) {
        AsI = const_cast<GenericSparseMatrix<double> *>(&A);
      } else {
        AsI = new GenericSparseMatrix<double>(A);
        AsI->Layout(GenericSparseMatrix<double>::CCS);
      }
      AsI->GetRawData(Ai, Ap, Ax);
      status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
      if (status != UMFPACK_OK) {
        cerr << "GenericSparseMatrix::Eigenvalues: UMFPACK symbolic analysis failed (status="
                  << status << "): " << umfpack_status_message(status) << endl;
        exit(1);
      }
      status = umfpack_di_numeric(Ap, Ai, Ax, symbolic, &Alu, NULL, NULL);
      if (status != UMFPACK_OK) {
        cerr << "GenericSparseMatrix::Eigenvalues: UMFPACK factorization failed (status="
                  << status << "): " << umfpack_status_message(status) << endl;
        exit(1);
      }
    } else {
      AsI = new GenericSparseMatrix<double>(A);
      Vector d = AsI->Diag();
      for (int i = 0; i < n; ++i) d(i) -= sigma;
      AsI->Diag(d);
      AsI->Layout(GenericSparseMatrix<double>::CCS);
      AsI->GetRawData(Ai, Ap, Ax);
      status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &symbolic, NULL, NULL);
      if (status != UMFPACK_OK) {
        cerr << "GenericSparseMatrix::Eigenvalues: UMFPACK symbolic analysis failed (status="
                  << status << "): " << umfpack_status_message(status) << endl;
        exit(1);
      }
      status = umfpack_di_numeric (Ap, Ai, Ax, symbolic, &Alu, NULL, NULL);
      if (status != UMFPACK_OK) {
        cerr << "GenericSparseMatrix::Eigenvalues: UMFPACK factorization failed (status="
                  << status << "): " << umfpack_status_message(status) << endl;
        exit(1);
      }
    }
    umfpack_di_free_symbolic(&symbolic);
  }

  // Prepare working memory
  int    lworkl = (issym ? (p * (p + 8)) : (3 * p * (p + 2)));
  double *workl = Allocate<double>(lworkl);
  double *workd = Allocate<double>(3 * n);
  double *basis = Allocate<double>(n * p);
  char    bmat  = (mode == 3) ? 'G' : 'I';
  int ipntr[14] = {0}, iparam[11] = {0};
  int info  = 1; // info  != 0 ensures that our initial resid is used
  iparam[0] = 1; // ishift = 1 ensures we are never asked to handle ido=3
  iparam[2] = maxit;
  iparam[3] = 1; // block size must be one (according to ARPACK++)
  iparam[6] = mode;

  // Iterate until ARPACK's reverse communication parameter tells us to stop
  int ido = 0;
  while (ido != 99) {
    if (issym) {
      ARPACK_F77NAME(dsaupd)(&ido, &bmat, &n, which, &k, &tol, resid, &p, basis,
                             &n, iparam, ipntr, workd, workl, &lworkl, &info);
    } else {
      ARPACK_F77NAME(dnaupd)(&ido, &bmat, &n, which, &k, &tol, resid, &p, basis,
                             &n, iparam, ipntr, workd, workl, &lworkl, &info);
    }
    if (info < 0) {
      // TODO: Map error code to error string
      cerr << "GenericSparseMatrix::Eigenvalues: ARPACK d?aupd error (info=" << info << ")" << endl;
      exit(1);
    }
    switch (ido) {
      case -1: case 1:
        // Compute y = OP * x where
        if (mode == 3) {
          // y = inv(A - sigma I) x <=> (A - sigma I) y = x
          status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax,
                                    workd + ipntr[1] - 1, workd + ipntr[0] - 1,
                                    Alu, NULL, NULL);
          if (status != UMFPACK_OK) {
            cerr << "GenericSparseMatrix::Eigenvalues: UMFPACK solver failed (status="
                      << status << "): " << umfpack_status_message(status) << endl;
            exit(1);
          }
        } else {
          // y = A x
          A.MultAv(workd + ipntr[0] - 1, workd + ipntr[1] - 1);
        }
        break;
      case 2:
        // Compute y = B x, here: B = I
        memcpy(workd + ipntr[1] - 1, workd + ipntr[0] - 1, n * sizeof(double));
        break;
      case 99:
        // ARPACK converged
        break;
      default:
        cerr << "GenericSparseMatrix::Eigenvalues: Unexpected ARPACK action ido=" << ido << endl;
        exit(1);
    }
  }

  // Free no longer needed factorization
  Ap = NULL, Ai = NULL, Ax = NULL;
  umfpack_di_free_numeric(&Alu), Alu = NULL;
  if (AsI != &A) Delete(AsI);

  // Optionally return final residual vector
  if (v0) {
    v0->Resize(n);
    for (int i = 0; i < n; ++i) v0->Put(i, resid[i]);
  }

  // Number of converged Ritz values
  nconv = min(k, iparam[4]);

  // Compute eigenvalues (and -vectors) given the found Arnolid basis
  int rvec = static_cast<int>(E != NULL);

  v.Initialize(nconv);
  if (rvec == 1) E->Initialize(n, nconv);

  int    *select = Allocate<int   >(p);
  double *d      = Allocate<double>(k);

  if (issym) {
    char HowMny[] = "A";
    ARPACK_F77NAME(dseupd)(&rvec, HowMny, select, d, basis, &n, &sigma, &bmat, &n,
                           which, &k, &tol, resid, &p, basis, &n, iparam,
                           ipntr, workd, workl, &lworkl, &info);
    if (info != 0) {
      // TODO: Map error code to error string
      cerr << "GenericSparseMatrix::Eigenvalues: ARPACK dseupd error (info=" << info << ")" << endl;
      exit(1);
    }
    if (strcmp(which, "LM") == 0 || strcmp(which, "LA") == 0) {
      for (int i = nconv-1, c = 0; i >= 0; --i, ++c) {
        v(c) = d[i];
        if (rvec == 1) {
          for (int r = 0; r < n; ++r) {
            E->Put(r, c, basis[i*n + r]);
          }
        }
      }
    }
    if ((strcmp(which, "SM") == 0 || strcmp(which, "SA") == 0) && rvec == 0) {
      for (int i = nconv-1, c = 0; i >= 0; --i, ++c) {
        v(c) = d[i];
      }
    }
  } else {
    char   HowMny[] = "A";
    double sigmai   = .0;
    double *di      = Allocate<double>(k);
    double *workev  = Allocate<double>(3 * p);
    ARPACK_F77NAME(dneupd)(&rvec, HowMny, select, d, di, basis, &n, &sigma, &sigmai,
                           workev, &bmat, &n, which, &k, &tol, resid, &p,
                           basis, &n, iparam, ipntr, workd, workl, &lworkl, &info);
    Deallocate(di);
    Deallocate(workev);
    if (info != 0) {
      // TODO: Map error code to error string
      cerr << "GenericSparseMatrix::Eigenvalues: ARPACK dneupd error (info=" << info << ")" << endl;
      exit(1);
    }
    if (rvec) {
      for (int c = 0; c < nconv; ++c) {
        v(c) = d[c];
        for (int r = 0; r < n; ++r) {
          E->Put(r, c, basis[c*n + r]);
        }
      }
    } else if (nconv == k) {
      for (int i = k+1, c = 0; i > 0; --i, ++c) {
        v(c) = d[i];
      }
    } else {
      for (int i = nconv-1, c = 0; i >= 0; --i, ++c) {
        v(c) = d[i];
      }
    }
  }

  // Clean up
  Deallocate(d);
  Deallocate(select);
  Deallocate(resid);
  Deallocate(workd);
  Deallocate(workl);
  Deallocate(basis);

#else // MIRTK_Numerics_WITH_eigs
  cerr << "eigs: Only available when Numerics module was built with ARPACK and UMFPACK" << endl;
  exit(1);
#endif // MIRTK_Numerics_WITH_eigs

  return nconv;
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>
::Eigenvalues(Vector &, int, const char *, int, double, int, Vector *) const
{
  cerr << "GenericSparseMatrix::Eigenvalues: Only implemented for sparse double precision matrices" << endl;
  exit(1);
}

template <>
int GenericSparseMatrix<double>
::Eigenvalues(Vector &v, int k, const char *sigma, int p, double tol, int maxit, Vector *v0) const
{
  return eigs(*this, NULL, v, k, sigma, p, tol, maxit, v0);
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>
::Eigenvectors(Matrix &, int, const char *, int, double, int, Vector *) const
{
  cerr << "GenericSparseMatrix::Eigenvectors: Only implemented for sparse double precision matrices" << endl;
  exit(1);
}

template <>
int GenericSparseMatrix<double>
::Eigenvectors(Matrix &E, int k, const char *sigma, int p, double tol, int maxit, Vector *v0) const
{
  Vector v;
  return eigs(*this, &E, v, k, sigma, p, tol, maxit, v0);
}

// -----------------------------------------------------------------------------
template <class TEntry>
int GenericSparseMatrix<TEntry>
::Eigenvectors(Matrix &, Vector &, int, const char *, int, double, int, Vector *) const
{
  cerr << "GenericSparseMatrix::Eigenpairs: Only implemented for sparse double precision matrices" << endl;
  exit(1);
}

template <>
int GenericSparseMatrix<double>
::Eigenvectors(Matrix &E, Vector &v, int k, const char *sigma, int p, double tol, int maxit, Vector *v0) const
{
  return eigs(*this, &E, v, k, sigma, p, tol, maxit, v0);
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GenericSparseMatrix<float>;
template class GenericSparseMatrix<double>;


} // namespace mirtk
