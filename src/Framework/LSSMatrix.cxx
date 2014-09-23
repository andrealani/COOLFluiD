// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/LSSMatrix.hh"

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

LSSMatrix::LSSMatrix() : m_useGPU(false) 
{
}

//////////////////////////////////////////////////////////////////////////////

LSSMatrix::~LSSMatrix() {}

//////////////////////////////////////////////////////////////////////////////

void LSSMatrix::flushAssembly ()
{
  beginAssembly(FLUSH_ASSEMBLY);
  endAssembly(FLUSH_ASSEMBLY);
}

//////////////////////////////////////////////////////////////////////////////

void LSSMatrix::finalAssembly()
{
  beginAssembly(FINAL_ASSEMBLY);
  endAssembly(FINAL_ASSEMBLY);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework
} // namespace COOLFluiD
