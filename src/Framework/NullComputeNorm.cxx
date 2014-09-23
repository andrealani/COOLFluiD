// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullComputeNorm.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullComputeNorm, ComputeNorm, FrameworkLib, 1>
NullComputeNormProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void NullComputeNorm::GR_Combine (const CFreal & S1,
        const CFreal & S2,
        CFreal& D)
{
  D = S1 + S2;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NullComputeNorm::GR_GetLocalValue () const
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

RealVector NullComputeNorm::compute ()
{
  m_residuals = 0.;
  return m_residuals;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
