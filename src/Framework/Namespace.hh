// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Namespace_hh
#define COOLFluiD_Framework_Namespace_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class describes the names of the singletons that ar activated when
/// this Namespace is entered, thus specifying a sound context for the
/// methods to operate.
/// @author Tiago Quintino
class Framework_API Namespace : public Config::ConfigObject {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  Namespace(const std::string& name);

  /// Destructor
  ~Namespace();

  /// Sets the MeshData name associated to this Namespace
  /// @param theValue name of MeshData
  void setMeshDataName(const std::string& theValue);

  /// Gets the MeshData name associated to this Namespace
  /// @return name of MeshData
  std::string getMeshDataName() const;

  /// Sets the PhysicalModel name associated to this Namespace
  /// @param theValue name of PhysicalModel
  void setPhysicalModelName(const std::string& theValue);

  /// Gets the PhysicalModel name associated to this Namespace
  /// @return name of PhysicalModel
  std::string getPhysicalModelName() const;

  /// Sets the type of the PhysicalModel associated to this Namespace
  /// @param theValue type of the PhysicalModel
  void setPhysicalModelType(const std::string& theValue);

  /// Gets the type of the PhysicalModel associated to this Namespace
  /// @return name of PhysicalModel
  std::string getPhysicalModelType() const;

  /// Sets the SubSystemStatus name associated to this Namespace
  /// @param theValue name of SubSystemStatus
  void setSubSystemStatusName(const std::string& theValue);

  /// Gets the SubSystemStatus name associated to this Namespace
  /// @return name of SubSystemStatus
  std::string getSubSystemStatusName() const;
  
  /// @return the flag indicating if this namespace has to be used for coupling
  bool isForCoupling() const {return m_isForCoupling;}
  
private: // member data

  /// MeshData to be activated by this Namespace
  std::string m_MeshDataName;

  /// Name PhysicalModel to be activated by this Namespace
  std::string m_PhysicalModelName;

  /// Type of the PhysicalModel to be activated by this Namespace
  std::string m_PhysicalModelType;

  /// SubSystemStatus to be activated by this Namespace
  std::string m_SubSystemStatusName;
  
  /// flag indicating that the namespace has to be used for coupling purposes
  bool m_isForCoupling;
  
}; // end of class Namespace

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

// #include "Namespace.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Namespace_hh
