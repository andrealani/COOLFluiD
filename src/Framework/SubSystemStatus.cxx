// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/Factory.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void SubSystemStatus::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("TimeStepLayers","Number of Layers");
   options.addConfigOption< std::vector<CFreal> >("InnerDTRatio","Ratio of total timestep for each layer");
   options.addConfigOption< std::string >("ComputeDT","DT calculator.");
   options.addConfigOption< CFreal, Config::DynamicOption<> >("TimeStep","TimeStep value at start.");
}

//////////////////////////////////////////////////////////////////////////////

SubSystemStatus::SubSystemStatus(const std::string& name) :
  ConfigObject(name),
  m_var_registry ( new VarRegistry() ),
  m_iter(0),
  m_subIter(0),
  m_residual(0),
  m_monitored_var(0),
  m_subSystemName(),
  m_timeStep(-1.0),
  m_timeStepLayers(),
  m_innerDT(),
  m_innerDTRatio(),
  m_innerDTConf(0),
  m_previousTimeStep(0.),
  m_subIterationFlag(false),
  m_lastSubIter(false),
  m_stopSim(false)
{
  addConfigOptionsTo(this);

  m_max_time = 0.;

  m_timeStep = -1.0;
  setParameter("TimeStep",&m_timeStep);

  m_computeDTStr = "Null";
  setParameter("ComputeDT",&m_computeDTStr);

  m_timeStepLayers = 1;
  setParameter("TimeStepLayers",&m_timeStepLayers);

  setParameter("InnerDTRatio",&m_innerDTConf);

  m_residual = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

SubSystemStatus::~SubSystemStatus()
{
  deletePtr ( m_var_registry );
}

//////////////////////////////////////////////////////////////////////////////

void SubSystemStatus::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  // configure the DT computer
  m_computeDT.reset
    (FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeDT, m_computeDTStr)->
     create(m_computeDTStr));
  m_computeDT->setFactoryRegistry(getFactoryRegistry());
  
  configureNested ( m_computeDT.getPtr(), args );

  if (m_timeStepLayers > 1){
  m_innerDT.resize(m_timeStepLayers);
  m_innerDTRatio.resize(m_timeStepLayers);

  if(m_innerDT.size() != m_innerDTConf.size()) {
    CFout << "WARNING: Inner DeltaT not set correctly!!!" << "\n";
  }
  else {
      for (CFuint i = 0; i < m_timeStepLayers; ++i) {
        m_innerDTRatio[i] = m_innerDTConf[i];
        m_innerDT[i] = m_innerDTConf[i]*m_timeStep;
    }
  }
  }
}

//////////////////////////////////////////////////////////////////////////////

ConvergenceStatus SubSystemStatus::getConvergenceStatus()
{
  ConvergenceStatus status;

  status.iter    = getNbIter();
  status.subiter = getSubIter();
  status.res     = getResidual();
  status.time    = getCurrentTimeDim();

  return status;
}

//////////////////////////////////////////////////////////////////////////////

void SubSystemStatus::setResidual(const RealVector& residual)
{
  if(m_residual.size() != residual.size()) m_residual.resize(residual.size());

  m_residual = residual;
}

//////////////////////////////////////////////////////////////////////////////

CFreal SubSystemStatus::getResidual() const
{
  if (m_residual.size() == 0 || m_monitored_var >= m_residual.size())
  {
    return 0.0;
  }
  if (!m_global_res) 
  {
    return m_residual[m_monitored_var];
  }
  else 
  {
    CFreal globalResidual = 0.;
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

    for(CFuint iVar=0; iVar < nbEqs; ++iVar)
    {
      if(m_residual[iVar] != -MathTools::MathConsts::CFrealMax()) { globalResidual += m_residual[iVar];  }
    } 
    return globalResidual/nbEqs;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubSystemStatus::updateNbIter()
{
  ++m_iter;
}

//////////////////////////////////////////////////////////////////////////////

void SubSystemStatus::adimensionalizeTimeData()
{
  // adimensionalize m_timeStep
  m_timeStep /= (PhysicalModelStack::getActive()->getImplementor()->getRefTime());
}

//////////////////////////////////////////////////////////////////////////////

SubSystemStatusStack& SubSystemStatusStack::getInstance()
{
  static SubSystemStatusStack aSubSystemStatusStack;
  return aSubSystemStatusStack;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<SubSystemStatus> SubSystemStatusStack::getActive()
{
  cf_assert(getInstance().isEnabled());
  // returns the first on the stack
  return getInstance().top();
}

//////////////////////////////////////////////////////////////////////////////

SubSystemStatus* SubSystemStatusStack::createObject(const std::string& name)
{
  SubSystemStatus * ptr = new SubSystemStatus(name);
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

std::string
SubSystemStatusStack::getObjectName(const Common::SafePtr<Namespace>& nsp)
const
{
  return nsp->getSubSystemStatusName();
}

//////////////////////////////////////////////////////////////////////////////

CFreal SubSystemStatus::getCurrentTimeDim() const { return (PhysicalModelStack::getActive()->getImplementor()->getRefTime())*m_currentTime; }

void SubSystemStatus::setCurrentTimeDim(CFreal const p) { m_currentTime = p/(PhysicalModelStack::getActive()->getImplementor()->getRefTime()); }

CFreal SubSystemStatus::getMaxTimeDim() const { return m_max_time*(PhysicalModelStack::getActive()->getImplementor()->getRefTime()); }

void SubSystemStatus::setMaxTimeDim(CFreal time) { m_max_time = time/(PhysicalModelStack::getActive()->getImplementor()->getRefTime()); }

CFreal SubSystemStatus::getMaxDTDim() const { return m_maxDT*(PhysicalModelStack::getActive()->getImplementor()->getRefTime()); }

void SubSystemStatus::setMaxDTDim(CFreal p) { m_maxDT = p/(PhysicalModelStack::getActive()->getImplementor()->getRefTime()); }

CFreal SubSystemStatus::getDTDim() const { return m_timeStep*(PhysicalModelStack::getActive()->getImplementor()->getRefTime());  }

void SubSystemStatus::setDTDim(const CFreal DT)
{
  m_prevprevTimeStep = m_previousTimeStep;
  m_previousTimeStep = m_timeStep;
  m_timeStep = DT/(PhysicalModelStack::getActive()->getImplementor()->getRefTime());
}

//////////////////////////////////////////////////////////////////////////////

 void SubSystemStatus::setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr)
 {
   m_fr = fr;
 }

//////////////////////////////////////////////////////////////////////////////

 Common::SafePtr<Common::FactoryRegistry> SubSystemStatus::getFactoryRegistry() 
 {
#ifdef CF_HAVE_SINGLE_EXEC
   return m_fr;
#endif
 }

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

