// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/LSSVector.hh"

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

LSSVector::LSSVector() : m_useGPU(false)
{
}

//////////////////////////////////////////////////////////////////////////////

LSSVector::~LSSVector() {}

//////////////////////////////////////////////////////////////////////////////

void LSSVector::assembly()
  {
    beginAssembly();
    endAssembly();
  }

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework
} // namespace COOLFluiD










