// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SpaceMethodData.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void SpaceMethodData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("SolutionVar","Solution variable set.");
  options.addConfigOption< std::string >("UpdateVar","VarSet corresponding to Update variables.");
  options.addConfigOption< std::string >("DiffusiveVar","Diffusive variable set.");
  options.addConfigOption< bool >("FreezeSysMatrix","Freeze the system matrix in the sub-iterations?");
}

//////////////////////////////////////////////////////////////////////////////

SpaceMethodData::SpaceMethodData(Common::SafePtr<Method> owner)
  : MethodData(owner),
    _solutionVar(),
    _updateVar(),
    _diffusiveVar(), 
    _lss(),
    _numericalJacobian("Numerical"),
    _preconditionerData(),
    _isPerturb(false),
    m_isRestart(false),
    _iPerturbVar(0),
    _fillPreconditionerMatrix(false),
    _onlyPreprocessSolution(false),
    _computeJacobian(true),
    _sysMatFrozen(false)
{
  addConfigOptionsTo(this);
  
  _freezeSysMatEverIter = false;
  setParameter("FreezeSysMatrix",&_freezeSysMatEverIter);
  
  _updateVarStr = "Null";
  setParameter("UpdateVar",&_updateVarStr);
  
  _solutionVarStr = "Null";
  setParameter("SolutionVar",&_solutionVarStr);
  
  _diffusiveVarStr = "Null";
  setParameter("DiffusiveVar",&_diffusiveVarStr);
}

//////////////////////////////////////////////////////////////////////////////

SpaceMethodData::~SpaceMethodData()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethodData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "SpaceMethodData::configure() => before MethodData::configure()\n");
  
  MethodData::configure(args);
  
  CFLog(VERBOSE, "SpaceMethodData::configure() => after MethodData::configure()\n");
  
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  CFLog(VERBOSE, "SpaceMethodData::configure() => Configuring Update VarSet\n");
  std::string provName = (_updateVarStr != "Null") ? (physModel->getConvectiveName() + _updateVarStr) : "Null";
  CFLog(VERBOSE, "SpaceMethodData::configure() => VarSet is : " << provName << "\n");
  Common::SafePtr<ConvectiveVarSet::PROVIDER> varSetProv = getProvider<ConvectiveVarSet>(provName);
  CFLog(VERBOSE, "SpaceMethodData::configure() => before setting update VarSet ptr\n");
  _updateVar.reset(varSetProv->create(physModel->getImplementor()->getConvectiveTerm()));
  cf_assert(_updateVar.isNotNull());
  
  CFLogDebugMin("SpaceMethodData::configure() => Configuring Solution VarSet\n");
  provName = (_solutionVarStr != "Null") ? (physModel->getConvectiveName() + _solutionVarStr) : "Null";
  CFLogDebugMin("SpaceMethodData::configure() => VarSet is : " << provName << "\n");
  varSetProv = getProvider<ConvectiveVarSet>(provName);
  CFLog(VERBOSE, "SpaceMethodData::configure() => before setting solution VarSet ptr\n");
  _solutionVar.reset(varSetProv->create(physModel->getImplementor()->getConvectiveTerm()));
  cf_assert(_solutionVar.isNotNull());
  
  CFLogDebugMin("SpaceMethodData::configure() => Configuring Diffusive VarSet\n");
  provName = (physModel->getDiffusiveName() != "Null") ? (physModel->getDiffusiveName() + _diffusiveVarStr) : "Null";
  CFLogDebugMin("SpaceMethodData::configure() => VarSet is : " << provName << "\n");
  Common::SafePtr<DiffusiveVarSet::PROVIDER> diffVarSetProv = getProvider<DiffusiveVarSet>(provName);
  CFLog(VERBOSE, "SpaceMethodData::configure() => before setting diffusive VarSet ptr\n");
  _diffusiveVar.reset(diffVarSetProv->create(provName, physModel->getImplementor()));
  cf_assert(_diffusiveVar.isNotNull());
  configureNested ( _diffusiveVar.getPtr(), args );
  
  // configure the numerical jacobian
  configureNested(&_numericalJacobian, args);
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethodData::setup()
{
  MethodData::setup();

  RealVector refValues = PhysicalModelStack::getActive()->getImplementor()->getRefStateValues();
  _numericalJacobian.setRefValues(refValues);
}

//////////////////////////////////////////////////////////////////////////////

void SpaceMethodData::unsetup()
{
  MethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::CFL> SpaceMethodData::getCFL()
{
  return getConvergenceMethod()[0]->getCFL();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
