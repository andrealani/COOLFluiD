// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "PhysicalModelImpl.hh"
#include "Common/CFLog.hh"
#include "PhysicalPropertyLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelImpl::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("PropertyLibrary","Library computing physical properties.");
   options.addConfigOption< std::vector<CFreal> >("refValues","Reference values for variable scaling.");
   options.addConfigOption< CFreal >("refLength","Reference length for geometric scaling.");
   options.addConfigOption< bool >("Is2DHalf","Flag telling if the dimension is 2D and 1/2.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalModelImpl::PhysicalModelImpl(const std::string& name)
  : Common::OwnedObject(),
    ConfigObject(name),
    _jacobians(0),
    _isAdimensional(false),
    _refTime(),
    _refStateValuesConf(0),
    _refStateValues(),
    _refLength(),
    _eqSubSysDescriptor(),
    _physicalPropLib()
{
  addConfigOptionsTo(this);
  setParameter("refValues",&_refStateValuesConf);

  _refLength = 1.0;
  setParameter("refLength",&_refLength);

  _physicalPropLibStr = "Null";
  setParameter("PropertyLibrary",&_physicalPropLibStr);
  
  m_is2DHalf = false;
  setParameter("Is2DHalf",&m_is2DHalf);
}
    
//////////////////////////////////////////////////////////////////////////////

PhysicalModelImpl::~PhysicalModelImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelImpl::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  ConfigObject::configure(args);

  const CFuint nbEqs = getNbEquations();
  _refStateValues.resize(nbEqs);

  if(_refStateValues.size() != _refStateValuesConf.size()) {
    CFout << "WARNING: reference values not set !!!" << "\n";
    _refStateValues = 1.0;
  }
  else {
    for (CFuint i = 0; i < nbEqs; ++i) {
      _refStateValues[i] = _refStateValuesConf[i];
    }
  }

  _physicalPropLib = Environment::Factory<PhysicalPropertyLibrary>::getInstance().
    getProvider(_physicalPropLibStr)->create(_physicalPropLibStr);

  cf_assert(_physicalPropLib.isNotNull());

  configureNested ( _physicalPropLib.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModelImpl::setup()
{
  CFAUTOTRACE;

  CFLog(NOTICE,"-------------------------------------------------------------\n");
  CFLog(NOTICE,"Setting up Physical Model [" << getName() << "]\n" );
  // set the library
  _physicalPropLib->setup();
   // set the physical data
  computePhysicalData();
  // set the reference values
  setReferenceValues();
  // set reference value for time
  setReferenceTime();
  CFLog(NOTICE,"-------------------------------------------------------------\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

