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
#include <mirtkConfig.h>
#include <mirtkCommonConfig.h>

// Windows
#if WINDOWS
#  include <mirtkWindows.h>
#endif

// Standard containers
#include <mirtkPair.h>
#include <mirtkArray.h>
#include <mirtkArrayHeap.h>
#include <mirtkList.h>
#include <mirtkOrderedSet.h>
#include <mirtkOrderedMap.h>
#include <mirtkUnorderedSet.h>
#include <mirtkUnorderedMap.h>
#include <mirtkQueue.h>
#include <mirtkStack.h>

// Standard functions/types
#include <mirtkAssert.h>
#include <mirtkMemory.h>
#include <mirtkString.h>
#include <mirtkPath.h>
#include <mirtkMath.h>
#include <mirtkVersion.h>
#include <mirtkFastDelegate.h>
#include <mirtkStatus.h>
#include <mirtkAlgorithm.h>

// Base classes
#include <mirtkObject.h>
#include <mirtkObjectFactory.h>

#include <mirtkObserver.h>
#include <mirtkObservable.h>
#include <mirtkConfigurable.h>

// I/O
#include <mirtkIndent.h>
#include <mirtkCfstream.h>

// Testing
#include <mirtkTestProd.h>


#endif // MIRTK_Common_H
