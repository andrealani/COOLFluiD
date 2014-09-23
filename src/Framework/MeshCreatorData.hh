// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MeshCreatorData_hh
#define COOLFluiD_Framework_MeshCreatorData_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Framework/MethodCommand.hh"
#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class SpaceMethod;

//////////////////////////////////////////////////////////////////////////////

/// This class collects some data shared by all MeshCreator's
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API MeshCreatorData : public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  MeshCreatorData(Common::SafePtr<Method> owner);

  /// Default destructor
  virtual ~MeshCreatorData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Sets the File name
  void setFileName(const std::string& theValue)
  {
    m_fileName = theValue;
  }

  /// Gets the File name
  boost::filesystem::path getFileName() const
  {
    return m_fileName;
  }

  /// Gets the scaling factor for the mesh
  CFreal getMeshScalingFactor() const
  {
    return m_meshScalingFactor;
  }

protected: // data

  /// The mesh scaling factor
  CFreal m_meshScalingFactor;

  /// The name of the file to open
  std::string m_fileName;

}; // end of class MeshCreatorData

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MeshCreatorData_hh
