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

/**
 * \file  mirtk/Arpack.h
 * \brief Interface to ARPACK FORTRAN routines.
 *
 * \attention Include this header in internal files such as .cc translation units only!
 */

#ifndef MIRTK_Arpack_H
#define MIRTK_Arpack_H


#ifndef ARPACK_F77NAME
#  define ARPACK_F77NAME(x) x##_
#endif

extern "C"
{


// ----------------------------------------------------------------------------
// debug "common" statement

struct { 
  int logfil, ndigit, mgetv0;
  int msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
  int mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
  int mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
} ARPACK_F77NAME(debug);

// -----------------------------------------------------------------------------
// double precision symmetric routines

void ARPACK_F77NAME(dsaupd)(int *ido, char *bmat, int *n, char *which,
                            int *nev, double *tol, double *resid,
                            int *ncv, double *V, int *ldv,
                            int *iparam, int *ipntr, double *workd,
                            double *workl, int *lworkl, int *info);

void ARPACK_F77NAME(dseupd)(int *rvec, char *HowMny, int *select,
                            double *d, double *Z, int *ldz,
                            double *sigma, char *bmat, int *n,
                            char *which, int *nev, double *tol,
                            double *resid, int *ncv, double *V,
                            int *ldv, int *iparam, int *ipntr,
                            double *workd, double *workl,
                            int *lworkl, int *info);

// -----------------------------------------------------------------------------
// double precision nonsymmetric routines

void ARPACK_F77NAME(dnaupd)(int *ido, char *bmat, int *n, char *which,
                            int *nev, double *tol, double *resid,
                            int *ncv, double *V, int *ldv,
                            int *iparam, int *ipntr, double *workd,
                            double *workl, int *lworkl, int *info);

void ARPACK_F77NAME(dneupd)(int *rvec, char *HowMny, int *select,
                            double *dr, double *di, double *Z,
                            int *ldz, double *sigmar,
                            double *sigmai, double *workev,
                            char *bmat, int *n, char *which,
                            int *nev, double *tol, double *resid,
                            int *ncv, double *V, int *ldv,
                            int *iparam, int *ipntr,
                            double *workd, double *workl,
                            int *lworkl, int *info);

// -----------------------------------------------------------------------------
// single precision symmetric routines

void ARPACK_F77NAME(ssaupd)(int *ido, char *bmat, int *n, char *which,
                            int *nev, float *tol, float *resid,
                            int *ncv, float *V, int *ldv,
                            int *iparam, int *ipntr, float *workd,
                            float *workl, int *lworkl, int *info);

void ARPACK_F77NAME(sseupd)(int *rvec, char *HowMny, int *select,
                            float *d, float *Z, int *ldz,
                            float *sigma, char *bmat, int *n,
                            char *which, int *nev, float *tol,
                            float *resid, int *ncv, float *V,
                            int *ldv, int *iparam, int *ipntr,
                            float *workd, float *workl,
                            int *lworkl, int *info);

// -----------------------------------------------------------------------------
// single precision nonsymmetric routines

void ARPACK_F77NAME(snaupd)(int *ido, char *bmat, int *n, char *which,
                            int *nev, float *tol, float *resid,
                            int *ncv, float *V, int *ldv,
                            int *iparam, int *ipntr, float *workd,
                            float *workl, int *lworkl, int *info);

void ARPACK_F77NAME(sneupd)(int *rvec, char *HowMny, int *select,
                            float *dr, float *di, float *Z,
                            int *ldz, float *sigmar,
                            float *sigmai, float *workev, char *bmat,
                            int *n, char *which, int *nev,
                            float *tol, float *resid, int *ncv,
                            float *V, int *ldv, int *iparam,
                            int *ipntr, float *workd, float *workl,
                            int *lworkl, int *info);


} // extern "C"


#endif // MIRTK_Arpack_H
