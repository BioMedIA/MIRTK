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

#ifndef MIRTK_DataFiedelity_H
#define MIRTK_DataFiedelity_H

#include "mirtk/EnergyTerm.h"


namespace mirtk {


/**
 * Base class for energy terms measuring the amount of data fidelity
 *
 * Lower data fidelity corresponds to a better alignment of the registered
 * data. The optimizer should minimize the data fidelity.
 */
class DataFidelity : public EnergyTerm
{
  mirtkAbstractMacro(DataFidelity);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  DataFidelity(const char * = "", double = 1.0);

  /// Copy constructor
  DataFidelity(const DataFidelity &);

  /// Assignment operator
  DataFidelity &operator =(const DataFidelity &);

public:

  /// Destructor
  virtual ~DataFidelity();

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

};


} // namespace mirtk

#endif // MIRTK_EnergyTerm_H
