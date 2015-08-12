// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "EmptyConvergenceMethod/EmptyConvergenceMethod.hh"
#include "EmptyConvergenceMethod/EmptyIteratorData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace EmptyConvergenceMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<EmptyIteratorData>, 
		      EmptyIteratorData, 
		      EmptyConvergenceMethodLib> 
nullEmptyIteratorComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

//void EmptyIteratorData::defineConfigOptions(Config::OptionList& options)
//{
//}

//////////////////////////////////////////////////////////////////////////////

EmptyIteratorData::EmptyIteratorData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner)
{
  //  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

EmptyIteratorData::~EmptyIteratorData()
{
}

//////////////////////////////////////////////////////////////////////////////

void EmptyIteratorData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EmptyConvergenceMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

