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

#include "mirtk/GenericRegistrationDebugger.h"

#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Event.h"
#include "mirtk/Point.h"
#include "mirtk/Matrix.h"
#include "mirtk/GenericImage.h"
#include "mirtk/GenericRegistrationFilter.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/ImageSimilarity.h"

#ifdef HAVE_VTK
#  include "mirtk/Vtk.h"
#  include "vtkPoints.h"
#  include "vtkPointData.h"
#  include "vtkShortArray.h"
#  include "vtkFloatArray.h"
#  include "vtkStructuredGrid.h"
#  include "vtkXMLStructuredGridWriter.h"
#endif // HAVE_VTK

#ifdef HAVE_MIRTK_PointSet
#  include "vtkPointSet.h"
#  include "mirtk/PointSetIO.h"
#endif // HAVE_MIRTK_PointSet

#include "mirtk/CommonExport.h"

#include <cstdio>


namespace mirtk {


// global "debug" flag (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


// -----------------------------------------------------------------------------
void CopyString(char *out, size_t sz, const string &str)
{
#ifdef WINDOWS
  strcpy_s(out, sz, str.c_str());
#else
  strcpy(out, str.c_str());
#endif
}

// -----------------------------------------------------------------------------
template <class TReal>
void WriteGradient(const char *fname, FreeFormTransformation *ffd, int l, const TReal *g)
{
  GenericImage<TReal> gradient(ffd->Attributes(), 3);
  int xdof, ydof, zdof;
  for (int k = 0; k < ffd->Z(); ++k)
  for (int j = 0; j < ffd->Y(); ++j)
  for (int i = 0; i < ffd->X(); ++i) {
    ffd->IndexToDOFs(ffd->LatticeToIndex(i, j, k, l), xdof, ydof, zdof);
    gradient(i, j, k, 0) = g[xdof];
    gradient(i, j, k, 1) = g[ydof];
    gradient(i, j, k, 2) = g[zdof];
  }
  gradient.Write(fname);
}

// -----------------------------------------------------------------------------
static void WriteAsVTKDataSet(const char *fname, FreeFormTransformation *ffd, const double *g = NULL)
{
#ifdef HAVE_VTK
  vtkSmartPointer<vtkPoints>         pos  = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkShortArray>     stat = vtkSmartPointer<vtkShortArray>::New();
  vtkSmartPointer<vtkFloatArray>     coef = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray>     disp = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray>     grad = (g ? vtkSmartPointer<vtkFloatArray>::New() : NULL);
  vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

  pos->SetNumberOfPoints(ffd->NumberOfCPs());

  stat->SetName("status");
  stat->SetNumberOfComponents(1);
  stat->SetNumberOfTuples(ffd->NumberOfCPs());

  coef->SetName("coefficient");
  coef->SetNumberOfComponents(3);
  coef->SetNumberOfTuples(ffd->NumberOfCPs());

  disp->SetName("displacement");
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(ffd->NumberOfCPs());

  if (grad) {
    grad->SetName("gradient");
    grad->SetNumberOfComponents(3);
    grad->SetNumberOfTuples(ffd->NumberOfCPs());
  }

  int    i, j, k;
  double x1, y1, z1, x2, y2, z2;
  for (int cp = 0; cp < ffd->NumberOfCPs(); ++cp) {
    ffd->IndexToLattice(cp, i, j, k);
    x1 = i, y1 = j, z1 = k;
    ffd->LatticeToWorld(x1, y1, z1);
    pos->SetPoint(cp, x1, y1, z1);
    x2 = x1, y2 = y1, z2 = z1;
    ffd->Transform(x2, y2, z2);
    disp->SetTuple3(cp, x2 - x1, y2 - y1, z2 - z1);
    stat->SetTuple1(cp, ffd->IsActive(cp));
    ffd->Get(i, j, k, x2, y2, z2);
    coef->SetTuple3(cp, x2, y2, z2);
    if (grad) {
      ffd->IndexToDOFs(cp, i, j, k);
      grad->SetTuple3(cp, g[i], g[j], g[k]);
    }
  }

  grid->SetDimensions(ffd->X(), ffd->Y(), ffd->Z());
  grid->SetPoints(pos);
  grid->GetPointData()->SetScalars(stat);
  grid->GetPointData()->SetVectors(coef);
  grid->GetPointData()->AddArray(disp);
  if (grad) grid->GetPointData()->AddArray(grad);

  vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(fname);
  writer->SetCompressorTypeToZLib();
  SetVTKInput(writer, grid);
  writer->Update();
#endif // HAVE_VTK
}

// -----------------------------------------------------------------------------
void WriteTransformation(const char                *fname,
                         HomogeneousTransformation *lin,
                         const Point               &target_offset,
                         const Point               &source_offset)
{
  const Matrix mat = lin->GetMatrix();
  Matrix pre (4, 4);
  Matrix post(4, 4);
  pre .Ident();
  post.Ident();
  pre (0, 3) = - target_offset._x;
  pre (1, 3) = - target_offset._y;
  pre (2, 3) = - target_offset._z;
  post(0, 3) = + source_offset._x;
  post(1, 3) = + source_offset._y;
  post(2, 3) = + source_offset._z;
  lin->PutMatrix(post * mat * pre);
  lin->Write(fname);
  lin->PutMatrix(mat);
}

// -----------------------------------------------------------------------------
GenericRegistrationDebugger::GenericRegistrationDebugger(const char *prefix)
:
  _Prefix      (prefix),
  _LevelPrefix (true),
  _Registration(NULL)
{
}

// -----------------------------------------------------------------------------
GenericRegistrationDebugger::~GenericRegistrationDebugger()
{
}

// -----------------------------------------------------------------------------
void GenericRegistrationDebugger::HandleEvent(Observable *obj, Event event, const void *data)
{
  GenericRegistrationFilter * const r = _Registration;

  const int sz = 256;
  char prefix[sz];
  char suffix[sz];
  char fname [sz];

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Initialize/update debugger, set common file name prefix/suffix
  switch (event) {

    // -------------------------------------------------------------------------
    // Attach/detach debugger
    case RegisteredEvent:
      _Registration = dynamic_cast<GenericRegistrationFilter *>(obj);
      if (!_Registration) {
        cerr << "GenericRegistrationDebugger::HandleEvent: Cannot attach debugger to object which is not of type GenericRegistrationFilter" << endl;
        exit(1);
      }
      _Iteration = 0;
      break;
    case UnregisteredEvent:
      _Registration = NULL;
      break;

    // -------------------------------------------------------------------------
    // Start/end
    case StartEvent:
      _Level = reinterpret_cast<const struct Iteration *>(data)->Count() + 1; // Iter()
      if (_LevelPrefix) _Iteration = 0;
      // Get pointers to similarity terms
      _Similarity.clear();
      for (int i = 0; i < r->_Energy.NumberOfTerms(); ++i) {
        ImageSimilarity *similarity = dynamic_cast<ImageSimilarity *>(r->_Energy.Term(i));
        if (similarity) _Similarity.push_back(similarity);
      }
      // Do not add a break statement here!
    case EndEvent:
      snprintf(prefix, sz, "%slevel_%d_", _Prefix.c_str(), _Level);
      suffix[0] = '\0';
      break;

    // -------------------------------------------------------------------------
    // Iteration
    case IterationStartEvent:
    case IterationEvent:
      ++_Iteration;
      return; // No data to write yet
    case LineSearchIterationStartEvent:
      _LineIteration = reinterpret_cast<const struct Iteration *>(data)->Iter();
      return; // No data to write yet

    case LineSearchStartEvent:
      if (_LevelPrefix) {
        snprintf(prefix, sz, "%slevel_%d_",         _Prefix.c_str(), _Level);
        snprintf(suffix, sz, "_%03d",               _Iteration);
      } else {
        CopyString(prefix, sz, _Prefix);
        snprintf(suffix, sz, "_%03d",               _Iteration);
      }
      break;
    case AcceptedStepEvent:
      if (_LevelPrefix) {
        snprintf(prefix, sz, "%slevel_%d_",         _Prefix.c_str(), _Level);
        snprintf(suffix, sz, "_%03d_%03d_accepted", _Iteration, _LineIteration);
      } else {
        CopyString(prefix, sz, _Prefix);
        snprintf(suffix, sz, "_%03d_%03d_accepted", _Iteration, _LineIteration);
      }
      break;
    case RejectedStepEvent:
      if (_LevelPrefix) {
        snprintf(prefix, sz, "%slevel_%d_",         _Prefix.c_str(), _Level);
        snprintf(suffix, sz, "_%03d_%03d_rejected", _Iteration, _LineIteration);
      } else {
        CopyString(prefix, sz, _Prefix);
        snprintf(suffix, sz, "_%03d_%03d_rejected", _Iteration, _LineIteration);
      }
      break;

    // -------------------------------------------------------------------------
    // Ignored event
    default: return;
  }

  MultiLevelTransformation  *mffd = NULL;
  FreeFormTransformation    *ffd  = NULL;
  HomogeneousTransformation *lin  = NULL;

  if (r) {
    (mffd = dynamic_cast<MultiLevelTransformation  *>(r->_Transformation)) ||
    (ffd  = dynamic_cast<FreeFormTransformation    *>(r->_Transformation)) ||
    (lin  = dynamic_cast<HomogeneousTransformation *>(r->_Transformation));
  }
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (mffd->LocalTransformationIsActive(l)) {
        if (ffd) {
          ffd = NULL;
          break;
        }
        ffd = mffd->GetLocalTransformation(l);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Write debug information
  switch (event) {

    // -------------------------------------------------------------------------
    // Write initial state
    case StartEvent: {

      // Write input images and their derivatives
      for (size_t i = 0; i < r->_Image[r->_CurrentLevel].size(); ++i) {
        snprintf(fname, sz, "%simage_%02zu", prefix, i+1);
        r->_Image[r->_CurrentLevel][i].Write(fname);
        if (debug >= 2) {
          BaseImage *gradient = NULL;
          BaseImage *hessian  = NULL;
          for (size_t j = 0; j < _Similarity.size(); ++j) {
            if (_Similarity[j]->Target()->InputImage() == &r->_Image[r->_CurrentLevel][i]) {
              if (_Similarity[j]->Target()->PrecomputeDerivatives()) {
                gradient = _Similarity[j]->Target()->InputGradient();
                hessian  = _Similarity[j]->Target()->InputHessian();
              }
              break;
            }
            if (_Similarity[j]->Source()->InputImage() == &r->_Image[r->_CurrentLevel][i]) {
              if (_Similarity[j]->Source()->PrecomputeDerivatives()) {
                gradient = _Similarity[j]->Source()->InputGradient();
                hessian  = _Similarity[j]->Source()->InputHessian();
              }
              break;
            }
          }
          if (gradient) {
            snprintf(fname, sz, "%simage_%02zu_gradient", prefix, i+1);
            gradient->Write(fname);
          }
          if (hessian) {
            snprintf(fname, sz, "%simage_%02zu_hessian", prefix, i+1);
            hessian->Write(fname);
          }
        }
      }

      // Write input domain mask
      if (r->_Mask[r->_CurrentLevel]) {
        snprintf(fname, sz, "%smask", prefix);
        r->_Mask[r->_CurrentLevel]->Write(fname);
      }

      // Write input point set
      #ifdef HAVE_MIRTK_PointSet
        for (size_t i = 0; i < r->_PointSet[r->_CurrentLevel].size(); ++i) {
          vtkPointSet *pointset = r->_PointSet[r->_CurrentLevel][i];
          snprintf(fname, sz, "%spointset_%02zu%s", prefix, i+1, DefaultExtension(pointset));
          WritePointSet(fname, pointset);
        }
      #endif // HAVE_MIRTK_PointSet

    } break;

    // -------------------------------------------------------------------------
    // Write intermediate results after each gradient step
    case LineSearchStartEvent: {

      // Energy gradient vector
      const double * const gradient = reinterpret_cast<const LineSearchStep *>(data)->_Direction;

      // Write input and other debug output of energy terms
      r->_Energy.WriteDataSets(prefix, suffix, _Iteration == 1);

      if (debug >= 3) {

        // Write non-parametric gradient of data fidelity terms
        r->_Energy.WriteGradient(prefix, suffix);

        // Write energy gradient(s) w.r.t control points
        if (ffd) {
          if (ffd->T() > 1) {
            for (int l = 0; l < ffd->T(); ++l) {
              snprintf(fname, sz, "%senergy_gradient_t%02d%s", prefix, l+1, suffix);
              WriteGradient(fname, ffd, l, gradient);
            }
          } else {
            snprintf(fname, sz, "%senergy_gradient%s", prefix, suffix);
            WriteGradient(fname, ffd, 0, gradient);
          }
        } else if (mffd) {
          const double *g = gradient;
          for (int i = 0; i < mffd->NumberOfLevels(); ++i) {
            if (!mffd->LocalTransformationIsActive(i)) continue;
            ffd = mffd->GetLocalTransformation(i);
            if (ffd->T() > 1) {
              for (int l = 0; l < ffd->T(); ++l) {
                snprintf(fname, sz, "%senergy_gradient_wrt_ffd_%d_t%02d%s", prefix, i+1, l+1, suffix);
                WriteGradient(fname, ffd, l, gradient);
              }
            } else {
              snprintf(fname, sz, "%senergy_gradient_wrt_ffd_%d_%s", prefix, i+1, suffix);
              WriteGradient(fname, ffd, 0, g);
            }
            g += ffd->NumberOfDOFs();
          }
          ffd = NULL;
        } else if (lin) {
          snprintf(fname, sz, "%senergy_gradient%s.txt", prefix, suffix);
          ofstream of(fname);
          for (int dof = 0; dof < r->_Energy.NumberOfDOFs(); ++dof) {
            of << gradient[dof] << "\n";
          }
          of.close();
        }

      }

      // Write current transformation estimate
      snprintf(fname, sz, "%stransformation%s.dof.gz", prefix, suffix);
      if (lin) {
        WriteTransformation(fname, lin, r->_TargetOffset, r->_SourceOffset);
      } else {
        r->_Transformation->Write(fname);
        if (ffd && r->_Input.empty() && r->NumberOfPointSets() > 0 && debug >= 4) {
          snprintf(fname, sz, "%stransformation%s.vtp", prefix, suffix);
          WriteAsVTKDataSet(fname, ffd, gradient);
        }
      }
    } break;

    case AcceptedStepEvent:
    case RejectedStepEvent: {
      if (debug >= 5) {

        // Write updated input of data fidelity terms
        r->_Energy.WriteDataSets(prefix, suffix, false);

        // Write current transformation estimate
        snprintf(fname, sz, "%stransformation%s.dof.gz", prefix, suffix);
        if (lin) WriteTransformation(fname, lin, r->_TargetOffset, r->_SourceOffset);
        else     r->_Transformation->Write(fname);

      }
    } break;

    // -------------------------------------------------------------------------
    // Unhandled event
    default: break;
  }
}


} // namespace mirtk
