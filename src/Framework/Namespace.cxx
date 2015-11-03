// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Namespace.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void Namespace::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("PhysicalModelName","Name of PhysicalModel activated by this Namespace.");
   options.addConfigOption< std::string >("PhysicalModelType","Which PhysicalModel is activated by this Namespace.");
   options.addConfigOption< std::string >("MeshData","Which MeshData is activated by this Namespace.");
   options.addConfigOption< std::string >("SubSystemStatus","Which SusSystemStatus is activated by this Namespace.");
   options.addConfigOption< bool >("IsForCoupling","Flag indicating that the namespace has to be used for coupling (Default: false).");
}

//////////////////////////////////////////////////////////////////////////////

Namespace::Namespace(const std::string& name) :
  Config::ConfigObject(name),
  m_MeshDataName(),
  m_PhysicalModelName(),
  m_PhysicalModelType(),
  m_SubSystemStatusName()
{
   addConfigOptionsTo(this);
   
   setParameter("MeshData",&m_MeshDataName);
   setParameter("PhysicalModelName",&m_PhysicalModelName);
   setParameter("PhysicalModelType",&m_PhysicalModelType);
   setParameter("SubSystemStatus",&m_SubSystemStatusName); 
   
   m_isForCoupling = false;
   setParameter("IsForCoupling",&m_isForCoupling);
}

//////////////////////////////////////////////////////////////////////////////

Namespace::~Namespace()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string Namespace::getMeshDataName() const
{
  return m_MeshDataName;
}

//////////////////////////////////////////////////////////////////////////////

void Namespace::setMeshDataName(const std::string& theValue)
{
  m_MeshDataName = theValue;
}

//////////////////////////////////////////////////////////////////////////////

std::string Namespace::getPhysicalModelType() const
{
  return m_PhysicalModelType;
}

//////////////////////////////////////////////////////////////////////////////

void Namespace::setPhysicalModelType(const std::string& theValue)
{
  m_PhysicalModelType = theValue;
}

//////////////////////////////////////////////////////////////////////////////

std::string Namespace::getPhysicalModelName() const
{
  return m_PhysicalModelName;
}

//////////////////////////////////////////////////////////////////////////////

void Namespace::setPhysicalModelName(const std::string& theValue)
{
  m_PhysicalModelName = theValue;
}

//////////////////////////////////////////////////////////////////////////////

std::string Namespace::getSubSystemStatusName() const
{
  return m_SubSystemStatusName;
}

//////////////////////////////////////////////////////////////////////////////

void Namespace::setSubSystemStatusName(const std::string& theValue)
{
  m_SubSystemStatusName = theValue;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

