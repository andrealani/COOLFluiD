// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "NullRemeshCondition.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullRemeshCondition,
                      SimpleMeshAdapterData,
                      SimpleGlobalMeshAdapterModule>
NullRemeshConditionProvider("NullRemeshCondition");

//////////////////////////////////////////////////////////////////////////////

NullRemeshCondition::NullRemeshCondition(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void NullRemeshCondition::execute()
{

  getMethodData().setNeedRemeshing(true);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
