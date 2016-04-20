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

#ifndef MIRTK_RegistrationFilter_H
#define MIRTK_RegistrationFilter_H

#include "mirtk/Observable.h"


namespace mirtk {


// Forward declaration of output type
class Transformation;


/**
 * Base class for registration filters
 *
 * \code
 * MyRegistrationFilter registration;
 * Transformation      *transformation = NULL;
 * registration.Input(...);
 * registration.Output(&transformation);
 * registration.Read(config_file);
 * registration.Run();
 * \endcode
 */
class RegistrationFilter : public Observable
{
  mirtkAbstractMacro(RegistrationFilter);

private:

  /// Pointer to output transformation
  Transformation **_Output;

private:

  /// Copy constructor
  /// \note Intentionally not implemented
  RegistrationFilter(const RegistrationFilter &);

  /// Assignment operator
  /// \note Intentionally not implemented
  void operator =(const RegistrationFilter &);

protected:

  /// Constructor
  RegistrationFilter();

public:

  /// Destructor
  virtual ~RegistrationFilter() = 0;

  /// Read registration parameters from input stream
  virtual bool Read(const char *, bool = false);

  /// Read registration parameters from input stream
  virtual bool Read(istream &, bool = false) = 0;

  /// Write registration parameters to file
  virtual void Write(const char *) const = 0;

  /// Runs the registration filter
  virtual void Run() = 0;

  /// Set pointer to output transformation
  void Output(Transformation **);

  /// Get output transformation
  Transformation *Output();

protected:

  /// Set current output transformation
  bool Output(Transformation *);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline RegistrationFilter::RegistrationFilter()
{
}

// -----------------------------------------------------------------------------
inline RegistrationFilter::~RegistrationFilter()
{
}

// -----------------------------------------------------------------------------
inline bool RegistrationFilter::Read(const char *fname, bool echo)
{
  ifstream from(fname);
  if (!from) return false;
  bool ok = this->Read(from, echo);
  from.close();
  return ok;
}

// -----------------------------------------------------------------------------
inline void RegistrationFilter::Output(Transformation **output)
{
  _Output = output;
}

// -----------------------------------------------------------------------------
inline bool RegistrationFilter::Output(Transformation *output)
{
  if (_Output) {
    (*_Output) = output;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
inline Transformation *RegistrationFilter::Output()
{
  return _Output ? *_Output : NULL;
}


} // namespace mirtk

#endif // MIRTK_RegistrationFilter_H
