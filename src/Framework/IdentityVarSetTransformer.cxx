// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "IdentityVarSetTransformer.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<IdentityVarSetTransformer,
               VarSetTransformer,
               FrameworkLib,
               1>
identityVarSetTransformerProvider("Identity");

//////////////////////////////////////////////////////////////////////////////

IdentityVarSetTransformer::IdentityVarSetTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

IdentityVarSetTransformer::~IdentityVarSetTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
