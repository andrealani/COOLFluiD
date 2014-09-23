// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"
#include "SimpleGlobalMeshAdapter/InteractiveRemeshCondition.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InteractiveRemeshCondition,
                      SimpleMeshAdapterData,
                      SimpleGlobalMeshAdapterModule>
InteractiveRemeshConditionProvider("InteractiveRemeshCondition");

//////////////////////////////////////////////////////////////////////////////

void InteractiveRemeshCondition::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool, Config::DynamicOption<> >("Remesh","Need to remesh?");
}

//////////////////////////////////////////////////////////////////////////////

InteractiveRemeshCondition::InteractiveRemeshCondition(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
  addConfigOptionsTo(this);
  setParameter("Remesh",&m_remesh);
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveRemeshCondition::execute()
{
  CFAUTOTRACE;

  getMethodData().setNeedRemeshing(m_remesh);
  // set value back to false
  m_remesh = false;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
