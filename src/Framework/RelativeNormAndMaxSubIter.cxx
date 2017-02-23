// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/MathConsts.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/RelativeNormAndMaxSubIter.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/Framework.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RelativeNormAndMaxSubIter, StopCondition, FrameworkLib, 1>
relativeNormAndMaxSubIterCondProvider("RelativeNormAndMaxSubIter");

//////////////////////////////////////////////////////////////////////////////

void RelativeNormAndMaxSubIter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("RelativeNorm","Target value of the relative residual norm.");
   options.addConfigOption< CFuint >("MaxSubIter","Maximum number of sub-iterations to perform.");

   options.addConfigOption< std::string >("InterfaceName","Name of the Interface for the coupling residual.");
   options.addConfigOption< std::string >("TRSName","Name of the TRS for the coupling residual.");
   options.addConfigOption< std::string >("DataType","Type of the data contained in the coupling residual.");
   options.addConfigOption< std::string >("Namespace","Namespace of the data contained in the coupling residual.");
}

//////////////////////////////////////////////////////////////////////////////

RelativeNormAndMaxSubIter::RelativeNormAndMaxSubIter(const std::string& name) :
  StopCondition(name)
{
  addConfigOptionsTo(this);

  m_relNorm = -8.0;
  setParameter("RelativeNorm",&m_relNorm);

  m_firstNorm = MathTools::MathConsts::CFrealMax();

  m_maxIter = std::numeric_limits<CFuint>::max();
  setParameter("MaxSubIter",&m_maxIter);

///@todo here the name of the residual should be built from the TRS, Namespace (NamespaceSwitcher
  setParameter("InterfaceName",&m_interfaceName);

  setParameter("Namespace",&m_nspName);

  setParameter("TRSName",&m_trsName);

  setParameter("DataType",&m_dataType);

}

//////////////////////////////////////////////////////////////////////////////

RelativeNormAndMaxSubIter::~RelativeNormAndMaxSubIter()
{
}

//////////////////////////////////////////////////////////////////////////////

bool RelativeNormAndMaxSubIter::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool RelativeNormAndMaxSubIter::isAchieved(const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;
  const CFuint start = 0; // first residual comes after first iteration
  const CFuint first = 1; // first residual comes after first iteration
  
  const std::string subSysName = SubSystemStatusStack::getCurrentName();
  const std::string nspName = NamespaceSwitcher::getInstance(subSysName).getCurrentNamespace()->getName();
  const std::string residualName = 
    "COUPLING_" + m_interfaceName + "_" + m_trsName + "_" + nspName + "_" + subSysName +"_" + m_dataType + "_DATA";
  
  if(m_nspName == nspName){
    
  CFout << "Looking at residual : " << residualName << "\n";

  const CFreal res  = SimulationStatus::getInstance().getCouplingResidual(residualName);
  const CFuint subIter = SubSystemStatusStack::getActive()->getSubIter();

  if (subIter == start) { return (subIter >= m_maxIter); }

  if (subIter == first)
  {
    m_firstNorm = res;
  }

  const CFreal diff = res - m_firstNorm;

  return ( (diff < m_relNorm) || (subIter >= m_maxIter) );
}
return false;

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
