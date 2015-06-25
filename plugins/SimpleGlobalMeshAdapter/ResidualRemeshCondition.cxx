// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "ResidualRemeshCondition.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ResidualRemeshCondition,
                      SimpleMeshAdapterData,
                      SimpleGlobalMeshAdapterModule>
residualRemeshConditionProvider("ResidualRemeshCondition");

//////////////////////////////////////////////////////////////////////////////

void ResidualRemeshCondition::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinValue","Minimum quality acceptable.");
}

//////////////////////////////////////////////////////////////////////////////

ResidualRemeshCondition::ResidualRemeshCondition(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
   addConfigOptionsTo(this);

   _minValue = -5.0;
   setParameter("MinValue",&_minValue);

}

//////////////////////////////////////////////////////////////////////////////

void ResidualRemeshCondition::execute()
{
  CFAUTOTRACE;

  const std::string otherNamespace = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> otherNsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(otherNamespace);
  const CFreal residual = SubSystemStatusStack::getInstance().
    getEntryByNamespace(otherNsp)->getResidual();
  
  if(residual < _minValue) getMethodData().setNeedRemeshing(true);
  else getMethodData().setNeedRemeshing(false);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
