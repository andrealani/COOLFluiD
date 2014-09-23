// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileWriter_CFmeshWriterData_hh
#define COOLFluiD_IO_CFmeshFileWriter_CFmeshWriterData_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Common/CFMap.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/OutputFormatterData.hh"
#include "CFmeshFileWriter/CFmeshFileWriter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// CFmeshFileWriterCom 's that compose the CFmeshFileWriter.
/// @see CFmeshWriterCom
/// @author Tiago Quintino
class CFmeshFileWriter_API CFmeshWriterData : public Framework::OutputFormatterData {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);


  /// Default constructor without arguments
  CFmeshWriterData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~CFmeshWriterData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Gets the Class name
  static std::string getClassName() {  return "CFmeshWriter";  }

  /// Sets the filename
  /// @param filepath pth to the file to be open
  void setFilename(const boost::filesystem::path& filepath)
  {
    m_filepath = filepath;
  }

  /// Gets the filename
  /// @return path to the file
  boost::filesystem::path getFilename() const
  {
    return m_filepath;
  }

  /// Gets the names and tags of the extra state variables to be written to file
  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > >
  getExtraVarSocketNamesAndTags()
  {
    return &_extraVarMap;
  }
  
  /// Gets the names and tags of the extra nodal variables to be written to file
  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > >
  getExtraNVarSocketNamesAndTags()
  {
    return &_extraNodalVarMap;
  }

  /// Gets the names and tags of the extra state variables to be written to file
  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > >
  getExtraSVarSocketNamesAndTags()
  {
    return &_extraStateVarMap;
  }

  /// Gets the flag for storing or not the past states
  bool storePastStates()
  {
    return _storePastStates;
  }


  /// Gets the flag for storing or not the past nodes
  bool storePastNodes()
  {
    return _storePastNodes;
  }

  /// Gets the flag for storing or not the inter states
  bool storeInterStates()
  {
    return _storeInterStates;
  }


  /// Gets the flag for storing or not the inter nodes
  bool storeInterNodes()
  {
    return _storeInterNodes;
  }

private:

  /// Filename to write solution to.
  boost::filesystem::path m_filepath;

  /// Flag to store past States/Nodes
  bool _storePastStates;
  bool _storePastNodes;
  /// Flag to store past Inter/Nodes
  bool _storeInterStates;
  bool _storeInterNodes;

  /// Name of the extra variables (for configuration)
  std::vector<std::string> _extraVarNames;
  std::vector<std::string> _extraNodalVarNames;
  std::vector<std::string> _extraStateVarNames;

  /// Tags of the extra variables (for configuration)
  std::vector<std::string> _extraVarTags;
  std::vector<std::string> _extraNodalVarTags;
  std::vector<std::string> _extraStateVarTags;

  /// Strides size of the extra variables (for configuration)
  std::vector<CFuint> _extraVarStrides;
  std::vector<CFuint> _extraNodalVarStrides;
  std::vector<CFuint> _extraStateVarStrides;

  /// Map containing the mapping between extra variable name and extra var file tag and stride
  Common::CFMap<std::string,std::pair<std::string,CFuint> > _extraVarMap;
  Common::CFMap<std::string,std::pair<std::string,CFuint> > _extraNodalVarMap;
  Common::CFMap<std::string,std::pair<std::string,CFuint> > _extraStateVarMap;

}; // end of class CFmeshWriterData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for CFmeshWriter
typedef Framework::MethodCommand<CFmeshWriterData> CFmeshWriterCom;

/// Definition of a command provider for CFmeshWriter
typedef Framework::MethodCommand<CFmeshWriterData>::PROVIDER CFmeshWriterComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileWriter_CFmeshWriterData_hh
