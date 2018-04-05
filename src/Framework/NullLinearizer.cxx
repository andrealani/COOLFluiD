// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Framework/NullLinearizer.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullLinearizer,
               JacobianLinearizer,
               FrameworkLib,
               1>
nullLinearizerScalarProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullLinearizer::NullLinearizer(Common::SafePtr<Framework::PhysicalModel> model) :
JacobianLinearizer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

NullLinearizer::~NullLinearizer()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullLinearizer::setPhysicalModel(PhysicalModelImpl* model)
{
}

//////////////////////////////////////////////////////////////////////////////

void NullLinearizer::linearize(const std::vector<State*>& statesInCell)
{
  CFLog(VERBOSE,"Calling NullSplitter::linearize() => this is a Null Linearizer."
  << "\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
