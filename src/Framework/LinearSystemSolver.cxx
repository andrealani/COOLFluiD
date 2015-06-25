// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/LinearSystemSolver.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/LSSData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFuint> >("MaskEquationIDs","IDs of the equations to solve with the current LSS");
}

//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::build_dynamic_functions()
{
  add_dynamic_function("solveSys",&LinearSystemSolver::solveSys);
}

//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::registActionListeners()
{
  Method::registActionListeners();
  
  // const std::string ssname = SubSystemStatusStack::getCurrentName();  
  //   event_handler->addListener(event_handler->key(ssname, "CF_ON_LSS_SOLVE"),
  //                              this,&LinearSystemSolver::solve);
}

//////////////////////////////////////////////////////////////////////////////

LinearSystemSolver::LinearSystemSolver(const std::string& name)
  : Method(name),
    m_lssData(CFNULL),
    m_nbSysEquations(),
    m_maskArray()
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);

  m_maskEquationIDs = vector<CFuint>();
  setParameter("MaskEquationIDs",&m_maskEquationIDs);
}

//////////////////////////////////////////////////////////////////////////////

LinearSystemSolver::~LinearSystemSolver()
{
}


//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::setMethodImpl()
{
  if (getMethodData().isNotNull()) {
    m_lssData = getMethodData().d_castTo<LSSData>();
  }
}

//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  const CFuint totalNbEqs = physModel->getNbEq();
  m_nbSysEquations = totalNbEqs;
  m_maskArray.resize(totalNbEqs);
  if (m_maskEquationIDs.size() == 0) {
    cf_assert(m_nbSysEquations == totalNbEqs);
    cf_assert(m_maskArray.size() == totalNbEqs);
    m_maskEquationIDs.resize(totalNbEqs);
    for (CFuint i = 0; i < totalNbEqs; ++i) {
      m_maskEquationIDs[i] = i;
    }
    m_maskArray = true;
    cf_assert(m_nbSysEquations == totalNbEqs);
  }
  else {
    m_maskArray = false;
    for (CFuint iEq = 0; iEq < m_maskEquationIDs.size(); ++iEq) {
      m_maskArray[m_maskEquationIDs[iEq]] = true;
    }
    m_nbSysEquations = m_maskEquationIDs.size();
  }

  cf_assert(m_maskEquationIDs.size() == m_nbSysEquations);
}

//////////////////////////////////////////////////////////////////////////////

void LinearSystemSolver::solveSys()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  solveSysImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

LSSIdxMapping& LinearSystemSolver::getLocalToGlobalMapping() 
{
  return m_lssData->getLocalToGlobalMapping(); 
}

//////////////////////////////////////////////////////////////////////////////

Framework::LSSIdxMapping& LinearSystemSolver::getLocalToLocallyUpdatableMapping() 
{
  return m_lssData->getLocalToLocallyUpdatableMapping();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
