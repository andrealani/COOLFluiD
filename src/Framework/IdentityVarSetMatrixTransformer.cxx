// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/IdentityVarSetMatrixTransformer.hh"
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

Environment::ObjectProvider<IdentityVarSetMatrixTransformer,
               VarSetMatrixTransformer,
               FrameworkLib,
               1>
identityVarSetMatrixTransformerProvider("Identity");

//////////////////////////////////////////////////////////////////////////////

IdentityVarSetMatrixTransformer::IdentityVarSetMatrixTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

IdentityVarSetMatrixTransformer::~IdentityVarSetMatrixTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
