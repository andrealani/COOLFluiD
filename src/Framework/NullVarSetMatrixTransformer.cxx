// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NullVarSetMatrixTransformer.hh"
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

Environment::ObjectProvider<NullVarSetMatrixTransformer,
               VarSetMatrixTransformer,
               FrameworkLib,
               1>
nullTransformerProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullVarSetMatrixTransformer::NullVarSetMatrixTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

NullVarSetMatrixTransformer::~NullVarSetMatrixTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

vector<State*>* NullVarSetMatrixTransformer::transform
(vector<State*> *const  states)
{
  CFLog(VERBOSE, "Calling NullVarSetMatrixTransformer::transform() => " <<
	"this is a Null VarSetMatrixTransformer " << "\n");

  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

State* NullVarSetMatrixTransformer::transform(State* const state)
{
  CFLog(VERBOSE, "Calling NullVarSetMatrixTransformer::transform() => " <<
	"this is a Null VarSetMatrixTransformer " << "\n");

  return CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

bool NullVarSetMatrixTransformer::getIsIdentityTransformation() const
{
  CFLog(VERBOSE, "Calling NullVarSetMatrixTransformer::getIsIdentityTransformation() => " <<
	"this is a Null VarSetMatrixTransformer " << "\n");

  return false;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
