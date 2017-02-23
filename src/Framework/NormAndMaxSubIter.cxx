// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NormAndMaxSubIter.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SimulationStatus.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NormAndMaxSubIter, StopCondition, FrameworkLib, 1>
NormAndMaxSubIterCondProvider("NormAndMaxSubIter");

//////////////////////////////////////////////////////////////////////////////

void NormAndMaxSubIter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Norm","Target value of the residual norm.");
   options.addConfigOption< CFuint >("MaxSubIter","Maximum number of sub-iterations to perform.");
   options.addConfigOption< std::string >("Namespace","Name of the namespace to check residual.");
}

//////////////////////////////////////////////////////////////////////////////

NormAndMaxSubIter::NormAndMaxSubIter(const std::string& name) :
  StopCondition(name)
{
  addConfigOptionsTo(this);

  m_norm = -8.0;
  setParameter("Norm",&m_norm);

  m_maxIter = std::numeric_limits<CFuint>::max();
  setParameter("MaxSubIter",&m_maxIter);

  setParameter("Namespace",&m_nspName);
}

//////////////////////////////////////////////////////////////////////////////

NormAndMaxSubIter::~NormAndMaxSubIter()
{
}

//////////////////////////////////////////////////////////////////////////////

bool NormAndMaxSubIter::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool NormAndMaxSubIter::isAchieved(const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;

  const std::string currentNspName = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getCurrentNamespace()->getName();
  
  if(m_nspName == currentNspName){
    const CFreal res  = SubSystemStatusStack::getActive()->getResidual();
    const CFuint subIter = SubSystemStatusStack::getActive()->getSubIter();

    if (subIter == 0) { return (subIter >= m_maxIter); }

    return ( (res < m_norm) || (subIter >= m_maxIter) );
  }

return false;

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
