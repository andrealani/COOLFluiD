// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullVarSetTransformer.hh"
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

Environment::ObjectProvider<NullVarSetTransformer,
               VarSetTransformer,
               FrameworkLib,
               1>
nullVarSetTransformerProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullVarSetTransformer::NullVarSetTransformer(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

NullVarSetTransformer::~NullVarSetTransformer()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSetTransformer::transform(const State& state, State& result)
{
  CFLog(VERBOSE, "Calling NullVarSetTransformer::transform() => " <<
	"this is a Null VarSetTransformer " << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSetTransformer::transformFromRef(const RealVector& pdata, State& result)
{
  CFLog(VERBOSE, "Calling NullVarSetTransformer::transformFromRef() => " <<
  "this is a Null VarSetTransformer " << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
