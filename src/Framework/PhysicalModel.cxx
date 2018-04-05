// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

PhysicalModel::PhysicalModel(const std::string& name) :
  ConfigObject(name),
  _physicalModelImpl(),
  _setup(false),
  _dimension(0),
  _nbEquations(0),
  _nameObject(name),
  _nameImplementor(),
  _convectiveName(),
  _diffusiveName(),
  _sourceName()
{
}

//////////////////////////////////////////////////////////////////////////////

PhysicalModel::~PhysicalModel()
{
  unsetPhysicalModelImpl();
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModel::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModel::unsetPhysicalModelImpl()
{
  if(_setup) {
    _physicalModelImpl.release();
    _dimension      = 0;
    _nbEquations    = 0;
    _nameImplementor= "";
    _convectiveName = "";
    _diffusiveName  = "";
    _sourceName     = "";
    _setup = false;
    CFLog(VERBOSE,"PhysicalModel unset!" << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModel::setPhysicalModelImpl
(Common::SelfRegistPtr<PhysicalModelImpl> physicalModelImpl)
{
  if(_setup) {
    unsetPhysicalModelImpl();
  }

  // set the physical model implementor
  _physicalModelImpl = physicalModelImpl;
  cf_assert(_physicalModelImpl.isNotNull());

  _dimension         = _physicalModelImpl->getDimension();
  _nbEquations       = _physicalModelImpl->getNbEquations();
  _nameImplementor   = _physicalModelImpl->getName();
  _convectiveName    = _physicalModelImpl->getConvectiveName();
  _diffusiveName     = _physicalModelImpl->getDiffusiveName();
  //  _sourceName        = _physicalModelImpl->getSourceName();

  // set up the physical model implementor
  _setup = true;
    
  CFLog(VERBOSE,"Set PhysicalModel's name to: " << _nameImplementor << "\n");
  CFLog(VERBOSE,"Set PhysicalModel's convective name to: " << _convectiveName << "\n");
  CFLog(VERBOSE,"Set PhysicalModel's diffusive name to: "  << _diffusiveName << "\n");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalModel* PhysicalModelStack::createObject(const std::string& name)
{
  PhysicalModel * ptr = new PhysicalModel(name);
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

std::string
PhysicalModelStack::getObjectName(const Common::SafePtr<Namespace>& nsp) const
{
  return nsp->getPhysicalModelName();
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalModel::setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr)
 {
   m_fr = fr;
 }

//////////////////////////////////////////////////////////////////////////////

 Common::SafePtr<Common::FactoryRegistry> PhysicalModel::getFactoryRegistry() 
 {
#ifdef CF_HAVE_SINGLE_EXEC
  cf_assert(m_fr != CFNULL);
#endif
  return m_fr;
 }

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

