/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#ifndef MIRTK_Common_H
#define MIRTK_Common_H

// Configuration, preprocessor flags
#include "mirtk/Config.h"
#include "mirtk/CommonConfig.h"

// Standard containers
#include "mirtk/Exception.h"
#include "mirtk/Pair.h"
#include "mirtk/Array.h"
#include "mirtk/ArrayHeap.h"
#include "mirtk/List.h"
#include "mirtk/OrderedSet.h"
#include "mirtk/OrderedMap.h"
#include "mirtk/UnorderedSet.h"
#include "mirtk/UnorderedMap.h"
#include "mirtk/PriorityQueue.h"
#include "mirtk/Queue.h"
#include "mirtk/Stack.h"

// Standard functions/types
#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/String.h"
#include "mirtk/Path.h"
#include "mirtk/Version.h"
#include "mirtk/FastDelegate.h"
#include "mirtk/Status.h"
#include "mirtk/Algorithm.h"
#include "mirtk/Numeric.h"
#include "mirtk/Random.h"

// Base classes
#include "mirtk/Object.h"
#include "mirtk/ObjectFactory.h"

#include "mirtk/Observer.h"
#include "mirtk/Observable.h"
#include "mirtk/Configurable.h"

// I/O
#include "mirtk/Indent.h"
#include "mirtk/Stream.h"
#include "mirtk/Cfstream.h"

// Testing
#include "mirtk/TestProd.h"


#endif // MIRTK_Common_H
