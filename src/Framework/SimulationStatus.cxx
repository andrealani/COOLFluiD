// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.


#include "Framework/SimulationStatus.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

SimulationStatus::SimulationStatus() :
  m_iter(0),
  m_currentTime(0.),
  m_appendIter(false),
  m_restart(false),
  m_lastResidual(0.),
  m_residual(),
  m_couplingResiduals(0),
  m_couplingResidualNames(CFNULL),
  m_time(),
  m_subSystems(),
  m_lastOutputFiles(),
  m_lastOutputFileConst()
{
  m_couplingResidualNames = new Common::CFMap<std::string, CFuint>();
}

//////////////////////////////////////////////////////////////////////////////

SimulationStatus::~SimulationStatus()
{
}

//////////////////////////////////////////////////////////////////////////////

SimulationStatus& SimulationStatus::getInstance()
{
   static SimulationStatus singleton;
   return singleton;
}

//////////////////////////////////////////////////////////////////////////////

void SimulationStatus::resetResidual()
{
  fill(m_residual.begin(),m_residual.end(),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void SimulationStatus::resetTime()
{
  fill(m_time.begin(),m_time.end(),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void SimulationStatus::setSubSystems(const vector<std::string>& subSys)
{
  m_subSystems.clear();
  m_subSystems.resize(subSys.size());
  copy(subSys.begin(),subSys.end(),m_subSystems.begin());
  m_residual.resize(m_subSystems.size());
}

//////////////////////////////////////////////////////////////////////////////

void SimulationStatus::addCouplingResidualNames(const std::string residualName)
{
  const CFuint residualNameID = m_couplingResiduals.size();
  m_couplingResidualNames->insert(residualName,residualNameID);
  m_couplingResiduals.push_back(0.0);
}

//////////////////////////////////////////////////////////////////////////////

void SimulationStatus::setCouplingResidual(const CFreal residual, const std::string residualName)
{
  const CFuint residualID = m_couplingResidualNames->find(residualName);
  m_couplingResiduals[residualID] = residual;
}

//////////////////////////////////////////////////////////////////////////////

CFreal SimulationStatus::getCouplingResidual(const std::string residualName)
{
  const CFuint residualID = m_couplingResidualNames->find(residualName);
  return m_couplingResiduals[residualID];
}

//////////////////////////////////////////////////////////////////////////////

CFreal SimulationStatus::getSimulationTime(const std::string subSysName)
{
  CFAUTOTRACE;

  ///@todo

  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

