/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/GenericRegistrationLogger.h"

#include "mirtk/Math.h"
#include "mirtk/Event.h"
#include "mirtk/Matrix.h"
#include "mirtk/Terminal.h"
#include "mirtk/EnergyTerm.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/GenericRegistrationFilter.h"
#include "mirtk/ImageSimilarity.h"

#include "mirtk/CommonExport.h"


namespace mirtk {


// Global "debug_time" flag (cf. mirtkProfiling.cc)
MIRTK_Common_EXPORT extern int debug_time;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace GenericRegistrationLoggerUtils {


// -----------------------------------------------------------------------------
ostream &PrintNumber(ostream &os, double value)
{
  if (value != .0 && (abs(value) < 1.e-5 || abs(value) >= 1.e5)) {
    os << scientific << setprecision(5) << setw(16) << value; // e-0x
  } else os << fixed << setprecision(5) << setw(12) << value << "    ";
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintNormalizedNumber(ostream &os, double value)
{
  os << fixed << setprecision(5) << setw(8) << value;
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintWeight(ostream &os, double weight, int nterms)
{
  const int w = iround(abs(weight) * 100.0);
  if      (w ==   0)               os << "< 1";
  else if (w == 100 && nterms > 1) os << ">99";
  else                             os << fixed << setw(3) << w;
  return os;
}


} // namespace GenericRegistrationLoggerUtils
using namespace GenericRegistrationLoggerUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
GenericRegistrationLogger::GenericRegistrationLogger(ostream *stream)
:
  _AlwaysReportRawValue(true),
  _Verbosity  (0),
  _Stream     (stream),
  _Color      (stdout_color),
  _FlushBuffer(true)
{
}

// -----------------------------------------------------------------------------
GenericRegistrationLogger::~GenericRegistrationLogger()
{
}

// =============================================================================
// Logging
// =============================================================================

// -----------------------------------------------------------------------------
void GenericRegistrationLogger::HandleEvent(Observable *obj, Event event, const void *data)
{
  if (!_Stream) return;
  ostream &os = *_Stream;

  // Change/Remember stream output format
  const streamsize w = os.width(0);
  const streamsize p = os.precision(6);
  const ios::fmtflags io_flags = os.flags();

  // Registration filter
  GenericRegistrationFilter *reg = dynamic_cast<GenericRegistrationFilter *>(obj);
  if (reg == NULL) {
    cerr << "GenericRegistrationLogger::HandleEvent: Cannot log events of object which is not of type GenericRegistrationFilter" << endl;
    exit(1);
  }

  // Typed pointers to event data (which of these is valid depends on the type of event)
  const Iteration      *iter = reinterpret_cast<const Iteration      *>(data);
  const LineSearchStep *step = reinterpret_cast<const LineSearchStep *>(data);
  const char           *msg  = reinterpret_cast<const char           *>(data);

  // Note: endl forces flush of stream buffer! Use \n instead and only flush buffer
  //       at end of this function if _FlushBuffer flag is set.

  switch (event) {

    // Before GenericRegistrationFilter::Initialize
    case InitEvent: {
      if (_Verbosity > 0) {
        if (iter->Iter() == reg->NumberOfLevels() && reg->InitialGuess()) {
          const HomogeneousTransformation *lin;
          lin = dynamic_cast<const HomogeneousTransformation *>(reg->InitialGuess());
          if (lin) {
            os << "\nInitial guess:\n\n";
            lin->Print(os, Indent(2, 2));
          }
        }
        os << "\n";
      }
      os << "\nResolution level " << iter->Iter() << "\n\n";
    } break;

    // Start of optimization, after GenericRegistrationFilter::Initialize
    case StartEvent: {
      if (_Verbosity > 1) {
        os << "\nRegistration domain:\n\n";
        reg->_RegistrationDomain.Print(os, 2);
        if (reg->NumberOfImages() > 0) {
          os << "\nAttributes of downsampled images:\n";
          for (int n = 0; n < reg->NumberOfImages(); ++n) {
            os << "\n";
            reg->_Image[reg->_CurrentLevel][n].Attributes().Print(os, 2);
          }
          const int nterms = reg->_Energy.NumberOfTerms();
          for (int t = 0; t < nterms; ++t) {
            const ImageSimilarity * const sim = dynamic_cast<ImageSimilarity *>(reg->_Energy.Term(t));
            if (sim) {
              os << "\nAttributes of registered images (" << sim->Name() << "):\n\n";
              sim->Domain().Print(os, 2);
            }
          }
        }
        os << "\nInitial transformation:\n\n";
        HomogeneousTransformation *lin = dynamic_cast<HomogeneousTransformation *>(reg->_Transformation);
        if (lin) {
          const Matrix mat = lin->GetMatrix();
          Matrix pre (4, 4);
          Matrix post(4, 4);
          pre .Ident();
          post.Ident();
          pre (0, 3) = - reg->_TargetOffset._x;
          pre (1, 3) = - reg->_TargetOffset._y;
          pre (2, 3) = - reg->_TargetOffset._z;
          post(0, 3) = + reg->_SourceOffset._x;
          post(1, 3) = + reg->_SourceOffset._y;
          post(2, 3) = + reg->_SourceOffset._z;
          lin->PutMatrix(post * mat * pre);
          lin->Print(os, Indent(2, 2));
          lin->PutMatrix(mat);
        } else {
          reg->_Transformation->Print(os, Indent(2, 2));
        }
      }
      if (_Verbosity == 0) {
        if (_Color) os << xboldblack;
        os << "          Energy     ";
        string name;
        int    i;
        size_t w;
        for (i = 0; i < reg->_Energy.NumberOfTerms(); ++i) {
          if (reg->_Energy.Term(i)->Weight() != .0) {
            w = 17; // 17 instead of 16 to fit often used "Bending energy" header
            if (reg->_Energy.NumberOfActiveTerms() < 5 && (_AlwaysReportRawValue || reg->_Energy.Term(i)->DivideByInitialValue())) {
              w += 10;
            }
            name                   = reg->_Energy.Term(i)->Name();
            if (name.empty()) name = reg->_Energy.Term(i)->NameOfClass();
            if (name.length() > w-3) {
              name = name.substr(0, w-3), name += "...";
            }
            os << string(3 + (w - name.length())/2, ' ');
            os << name;
            os << string((w - name.length() + 1)/2, ' ');
          }
        }
        if (_Color) os << xreset;
      }
      _NumberOfGradientSteps = 0;
    } break;

    // Next gradient step
    case IterationEvent:
    case IterationStartEvent: {
      ++_NumberOfGradientSteps;
      os << "\n";
      if (debug_time || _Verbosity > 0) {
        os << "Iteration " << setw(3) << left << _NumberOfGradientSteps << right << " ";
      } else {
        os << " " << setw(3) << _NumberOfGradientSteps << " ";
      }
      if (debug_time) os << "\n";
    } break;

    // Start of line search given current gradient direction
    case LineSearchStartEvent: {

      reg->_Energy.Value(); // update cached value of individual terms if necessary

      if (_Verbosity > 0) {
        if (debug_time) os << "\n              ";
        os << "Line search with ";
        if (step->_Info && strcmp(step->_Info, "unbound") != 0) {
          os << step->_Info << " ";
        }
        os << "step length in [";
        os.unsetf(ios::floatfield);
        os << setprecision(2) << step->_MinLength << ", ";
        if (step->_Info && strcmp(step->_Info, "unbound") == 0) {
          os << "inf";
        } else {
          os << setprecision(2) << step->_MaxLength;
        }
        os << "]";
      }
      string name;
      int i, nterms = 0;
      bool divide_by_init = false;
      for (i = 0; i < reg->_Energy.NumberOfTerms(); ++i) {
        if (reg->_Energy.Term(i)->Weight() != .0) {
          ++nterms;
          if (reg->_Energy.Term(i)->DivideByInitialValue()) divide_by_init = true;
        }
      }
      if (_Verbosity > 0) {
        os << "\n\n  Energy = ";
        if (_Color) os << xboldblue;
        PrintNumber(os, step->_Current);
        if (_Color) os << xreset;
        for (i = 0; i < reg->_Energy.NumberOfTerms(); ++i) {
          if (reg->_Energy.IsActive(i)) {
            if (nterms > 1) {
              os << " = ";
              if (_AlwaysReportRawValue || reg->_Energy.Term(i)->DivideByInitialValue()) {
                PrintNormalizedNumber(os, reg->_Energy.Value(i));
                os << " (";
                PrintNumber(os, reg->_Energy.RawValue(i));
                os << ")";
              } else {
                if (divide_by_init) os << "          ";
                PrintNumber(os, reg->_Energy.Value(i));
                if (divide_by_init) os << " ";
              }
            }
            os << "    ";
            PrintWeight(os, reg->_Energy.Term(i)->Weight(), nterms);
            name                   = reg->_Energy.Term(i)->Name();
            if (name.empty()) name = reg->_Energy.Term(i)->NameOfClass();
            os << "% " << name << "\n";
            break;
          }
        }
        for (i = i + 1; i < reg->_Energy.NumberOfTerms(); ++i) {
          if (reg->_Energy.IsActive(i)) {
            os << "                            ";
            if (reg->_Energy.ExcludeConstraints() && reg->_Energy.IsConstraint(i)) {
              os << " ";
            } else {
              os << "+";
            }
            os << " ";
            if (_AlwaysReportRawValue || reg->_Energy.Term(i)->DivideByInitialValue()) {
              PrintNormalizedNumber(os, reg->_Energy.Value(i));
              os << " (";
              PrintNumber(os, reg->_Energy.RawValue(i));
              os << ")";
            } else {
              if (divide_by_init) os << "          ";
              PrintNumber(os, reg->_Energy.Value(i));
              if (divide_by_init) os << " ";
            }
            os << "    ";
            PrintWeight(os, reg->_Energy.Term(i)->Weight(), nterms);
            name                   = reg->_Energy.Term(i)->Name();
            if (name.empty()) name = reg->_Energy.Term(i)->DefaultName();
            os << "% " << name << "\n";
          }
        }
        if (_Color) os << xboldblack;
        os << "\n                 Energy        Step Length        Max. Delta\n\n";
        if (_Color) os << xreset;
      } else {
        if (_Color) os << xboldblue;
        PrintNumber(os, step->_Current);
        if (_Color) os << xreset;
        nterms = 0;
        if (reg->_Energy.ExcludeConstraints()) {
          for (i = 0; i < reg->_Energy.NumberOfTerms(); ++i) {
            if (reg->_Energy.IsActive(i) && !reg->_Energy.IsConstraint(i)) {
              if (reg->_Energy.NumberOfActiveTerms() < 5) {
                if (nterms++ == 0) {
                  os << " = ";
                } else {
                  os << " + ";
                }
                if (_AlwaysReportRawValue || reg->_Energy.Term(i)->DivideByInitialValue()) {
                  PrintNormalizedNumber(os, reg->_Energy.Value(i));
                  os << " (";
                  PrintNumber(os, reg->_Energy.RawValue(i));
                  os << ")";
                } else {
                  os << " ";
                  PrintNumber(os, reg->_Energy.Value(i));
                }
              } else {
                os << "    ";
                PrintNumber(os, reg->_Energy.RawValue(i));
              }
            }
          }
          for (i = 0; i < reg->_Energy.NumberOfTerms(); ++i) {
            if (reg->_Energy.IsActive(i) && reg->_Energy.IsConstraint(i)) {
              if (reg->_Energy.NumberOfActiveTerms() < 5) {
                if (nterms++ == 0) {
                  os << " = ";
                } else {
                  os << "   ";
                }
                if (_AlwaysReportRawValue || reg->_Energy.Term(i)->DivideByInitialValue()) {
                  PrintNormalizedNumber(os, reg->_Energy.Value(i));
                  os << " (";
                  PrintNumber(os, reg->_Energy.RawValue(i));
                  os << ")";
                } else {
                  os << " ";
                  PrintNumber(os, reg->_Energy.Value(i));
                }
              } else {
                os << "    ";
                PrintNumber(os, reg->_Energy.RawValue(i));
              }
            }
          }
        } else {
          for (i = 0; i < reg->_Energy.NumberOfTerms(); ++i) {
            if (reg->_Energy.IsActive(i)) {
              if (reg->_Energy.NumberOfActiveTerms() < 5) {
                if (nterms++ == 0) {
                  os << " = ";
                } else {
                  if (reg->_Energy.ExcludeConstraints() && reg->_Energy.IsConstraint(i)) {
                    os << "   ";
                  } else {
                    os << " + ";
                  }
                }
                if (_AlwaysReportRawValue || reg->_Energy.Term(i)->DivideByInitialValue()) {
                  PrintNormalizedNumber(os, reg->_Energy.Value(i));
                  os << " (";
                  PrintNumber(os, reg->_Energy.RawValue(i));
                  os << ")";
                } else {
                  os << " ";
                  PrintNumber(os, reg->_Energy.Value(i));
                }
              } else {
                os << "    ";
                PrintNumber(os, reg->_Energy.RawValue(i));
              }
            }
          }
        }
      }
      _NumberOfIterations = 0;
      _NumberOfSteps      = 0;
    } break;

    case LineSearchIterationStartEvent:
      if (_Verbosity > 0) {
        if (debug_time) os << "\nStep " << left << setw(3) << iter->Count() << "\n";
        else            os <<   "     " << setw(3) << iter->Count() << "   ";
      }
      break;

    case AcceptedStepEvent:
    case RejectedStepEvent: {
      if (_Verbosity > 0) {
        if (debug_time) os << "\n           ";
        if (_Color) os << (event == AcceptedStepEvent ? xgreen : xbrightred);
        PrintNumber(os, step->_Value ) << "  ";
        if (step->_Info && strcmp(step->_Info, "incremental") == 0) {
          PrintNumber(os, step->_Length) << "  ";
          PrintNumber(os, step->_Delta);
        } else {
          PrintNumber(os, step->_TotalLength + step->_Length) << "  ";
          PrintNumber(os, step->_TotalDelta  + step->_Delta);
        }
        if (_Color) os << xreset;
        else os << "    " << ((event == AcceptedStepEvent) ? "Accepted" : "Rejected");
        os << "\n";
      }
      if (event == AcceptedStepEvent) ++_NumberOfSteps;
      // Increment counter of actual line search iterations performed
      // Note that the Delta convergence criterium on the minimum maximum
      // DoF change can cause the line search to stop immediately
      ++_NumberOfIterations;
    } break;

    case LineSearchIterationEndEvent: {
      if (_Verbosity > 0 && iter->Count() == iter->Total()) {
        if (_Color) os << xboldblack;
        os << "\n              Maximum number of iterations exceeded\n";
        if (_Color) os << xreset;
      }
    } break;

    // End of line search
    case LineSearchEndEvent: {
      if (_Verbosity > 0) {
        // The minimum maximum DoF change convergence criterium kicked in immediately...
        if (_NumberOfIterations == 0) {
          os << "                              ";
          if (_Color) os << xboldred;
          PrintNumber(os, step->_Delta) << "\n";
        }
        // Report if line search was successful or no improvement
        if (_Color) os << xboldblack;
        if (_NumberOfSteps > 0) {
          if (_Verbosity > 0) {
            if (step->_Info && strcmp(step->_Info, "incremental") == 0) {
              os << "\n            Gradient scale = ";
              PrintNumber(os, step->_TotalLength) << " /";
              PrintNumber(os, step->_Unit) << "\n";
              os <<   "            Max. delta     = ";
              PrintNumber(os, step->_TotalDelta) << "\n";
            } else {
              os << "\n            Gradient scale = ";
              PrintNumber(os, step->_TotalLength / step->_Unit) << "\n";
            }
          }
        } else {
          if (step->_Delta < reg->_Optimizer->Delta()) {
            os << "\n         Converged, max. delta <= " << setprecision(1) << scientific << reg->_Optimizer->Delta() << "\n";
          } else {
            os << "\n         No further improvement within search range\n";
          }
        }
        if (_Color) os << xreset;
      }
    } break;

    // End of gradient step
    case IterationEndEvent: {
      if (_Verbosity > 0) {
        if (_Color) os << xboldblack;
        if (_NumberOfSteps > 0) {
          os << "            Energy slope   = ";
          PrintNumber(os, reg->_Optimizer->LastValuesSlope()) << "\n";
        }
        if (reg->_Optimizer->Converged()) {
          double epsilon = reg->_Optimizer->Epsilon();
          double value   = reg->_Optimizer->LastValues().back();
          if (epsilon < 0.) epsilon = abs(epsilon * value);
          if (abs(reg->_Optimizer->LastValuesSlope()) < epsilon) {
            os << "\n        Converged, energy slope < " << setprecision(1) << scientific << epsilon << "\n";
          }
        }
        os << "\n";
        if (_Color) os << xreset;
      }
    } break;

    // Energy function modified after convergence and optimization restarted
    case RestartEvent: {
      os << "\n";
      if (_Verbosity > 0) {
        os << "\nContinue with modified energy function\n";
      }
    } break;

    // End of optimization, but before GenericRegistrationFilter::Finalize
    case EndEvent: {
      HomogeneousTransformation *lin;
      lin = dynamic_cast<HomogeneousTransformation *>(reg->_Transformation);
      if (lin) {
        bool print_params;
        if (_Verbosity > 0) {
          os << "\nFinal transformation at level " << iter->Iter() << ":\n\n";
          print_params = true;
        } else {
          print_params = (iter->Iter() == 1);
          if (print_params) {
            os << "\n\n\nFinal " << ToPrettyString(reg->_CurrentModel) << ":\n\n";
          }
        }
        if (print_params) {
          const Matrix mat = lin->GetMatrix();
          Matrix pre (4, 4);
          Matrix post(4, 4);
          pre .Ident();
          post.Ident();
          pre (0, 3) = - reg->_TargetOffset._x;
          pre (1, 3) = - reg->_TargetOffset._y;
          pre (2, 3) = - reg->_TargetOffset._z;
          post(0, 3) = + reg->_SourceOffset._x;
          post(1, 3) = + reg->_SourceOffset._y;
          post(2, 3) = + reg->_SourceOffset._z;
          lin->PutMatrix(post * mat * pre);
          lin->Print(os, Indent(2, 2));
          lin->PutMatrix(mat);
        }
      }
    } break;

    // After GenericRegistrationFilter::Finalize which cleans up behind
    case FinishEvent: {
      os << "\n";
      if (_Verbosity > 0) {
        os << "Finished registration at level " << iter->Iter() << "\n";
      }
    } break;

    // Status message broadcasted by registration filter
    case StatusEvent: {
      if (_Color) os << xboldblack;
      os << msg;
      if (_Color) os << xreset;
    } break;

    // Log message broadcasted by registration filter
    case LogEvent: {
      if (_Verbosity > 0) os << msg;
    } break;

    default: break;
  }

  // Flush output buffer
  if (_FlushBuffer) os.flush();

  // Reset stream output format
  os.width(w);
  os.precision(p);
  os.flags(io_flags);
}


} // namespace mirtk
