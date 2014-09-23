// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileReader_CFmeshReaderData_hh
#define COOLFluiD_IO_CFmeshFileReader_CFmeshReaderData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "Framework/MeshCreatorData.hh"
#include "CFmeshFileReader/CFmeshFileReaderAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// CFmeshFileReaderCom 's that compose the CFmeshFileReader.
/// @see CFmeshReaderCom
/// @author Tiago Quintino
/// @author Andrea Lani
class CFmeshFileReader_API CFmeshReaderData : public Framework::MeshCreatorData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  CFmeshReaderData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~CFmeshReaderData();

  /// Gets the file name to convert from
  /// @return the file name stored in m_convertFromFileName
  boost::filesystem::path getConvertFromFileName() const
  {
    return boost::filesystem::path(m_convertFromFileName);
  }

  /// Gets the file name to convert from
  /// @return the file name stored in m_convertFromFileName
  void setConvertFromFileName(std::string filename)
  {
    m_convertFromFileName = filename;
  }

  /// solution file with states from which each process can restart
  std::string getStateSolutionFile() const
  {
   return m_solutionFile;
  } 
 
  /// Gets the flag if the mesh has to be translated
  bool isTranslated()
  {
    return m_isTranslated;
  }

  /// Gets the flag if the mesh has to be translated
  bool isAnisotropicScaling()
  {
    return m_isAnisotropicScaling;
  }

  /// Sets the flag if the mesh has to be translated
  void setIsTranslated(const bool isTranslated)
  {
    m_isTranslated = isTranslated;
  }

  /// Gets the vector of translation of the mesh
  std::vector<CFreal> getTranslationVector()
  {
    return m_translationVector;
  }

  /// Gets the vector of transformation of the mesh
  std::vector<CFreal> getAnisotropicScalingVector()
  {
    return m_scalingVector;
  }

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Get initial values to set in the State
  std::vector<CFreal> getInitValues() const
  {
    return m_initValues;
  }

  /// Get values for scaling the states while reading the CFmesh
const std::vector<CFreal>& getStateScalingValues() const
  {
    return m_stateScalingValues;
  }


  /// Get initial values IDs to set in the State
  std::vector<CFuint> getInitValuesIDs() const
  {
    return m_initValuesIDs;
  }

  /// Get an array specifying which State variable need to be initialized
  /// with the given values
  std::vector<bool> getUseInitValues() const
  {
    return m_useInitValues;
  }

  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > > getExtraNVarSocketNamesAndTags()
  {
    return &m_extraNodalVarMap;
  }

  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > > getExtraSVarSocketNamesAndTags()
  {
    return &m_extraStateVarMap;
  }
  
  Common::SafePtr<Common::CFMap<std::string, std::pair<std::string,CFuint> > > getExtraVarSocketNamesAndTags()
  {
    return &m_extraVarMap;
  }

  /// Gets the flag for storing or not the past states
  bool readPastStates()
  {
    return m_readPastStates;
  }

  /// Gets the flag for storing or not the past nodes
  bool readPastNodes()
  {
    return m_readPastNodes;
  }

  /// Gets the flag for storing or not the intermediate states
  bool readInterStates()
  {
    return m_readInterStates;
  }

  /// Gets the flag for storing or not the intermediate nodes
  bool readInterNodes()
  {
    return m_readInterNodes;
  }
 
  /// array driving a switch in node coordinates: e.g. [0,2,1] means to read [X,Y,Z] as [X,Z,Y]   
  std::vector<CFuint> getNodeSwitchIDs() const
  {
    return m_nodeSwitchIDs;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "CFmeshReader";
  }

private: // data

  /// the filename of the format to convert from
  std::string m_convertFromFileName;

  /// the name of the type of polynomial interpolation
  std::string m_polynomialTypeName;

  /// should the mesh be translated
  bool m_isTranslated;

  /// should an anisotropic scaling be applied on the mesh
  bool m_isAnisotropicScaling;

  /// Translation vector
  std::vector<CFreal> m_translationVector;

  /// Scaling vector (for non isotropic scaling)
  std::vector<CFreal> m_scalingVector;

  /// array of flags telling if to use the given initial
  /// values to initialize each variable in the State's
  std::vector<bool> m_useInitValues;

  /// array of initial values to set in the State's
  std::vector<CFreal> m_initValues;

  /// array of initial values IDs to set in the State's
  std::vector<CFuint> m_initValuesIDs;

  /// Flag to store past States/Nodes
  bool m_readPastStates;
  bool m_readPastNodes;
  
  /// Flag to store intermediate States/Nodes
  bool m_readInterStates;
  bool m_readInterNodes;

  /// Name of the extra variables (for configuration)
  std::vector<std::string> m_extraNodalVarNames;
  std::vector<std::string> m_extraStateVarNames;
  std::vector<std::string> m_extraVarNames;

  /// Tags of the extra variables (for configuration)
  std::vector<std::string> m_extraNodalVarTags;
  std::vector<std::string> m_extraStateVarTags;
  std::vector<std::string> m_extraVarTags;

  /// Strides size of the extra variables (for configuration)
  std::vector<CFuint> m_extraNodalVarStrides;
  std::vector<CFuint> m_extraStateVarStrides;
  std::vector<CFuint> m_extraVarStrides;

  std::vector<CFreal>  m_stateScalingValues;

  /// Map containing the mapping between extra variable name and extra var file tag
  Common::CFMap<std::string,std::pair<std::string,CFuint> > m_extraNodalVarMap;
  Common::CFMap<std::string,std::pair<std::string,CFuint> > m_extraStateVarMap;
  Common::CFMap<std::string,std::pair<std::string,CFuint> > m_extraVarMap;

  /// name of solution file with only states from which each process can restart
  std::string m_solutionFile;

  /// array driving a switch in node coordinates: e.g. [0,2,1] means to read [X,Y,Z] as [X,Z,Y]
  std::vector<CFuint> m_nodeSwitchIDs; 

}; // end of class CFmeshReaderData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for CFmeshReader
typedef Framework::MethodCommand<CFmeshReaderData> CFmeshReaderCom;

/// Definition of a command provider for CFmeshReader
typedef Framework::MethodCommand<CFmeshReaderData>::PROVIDER CFmeshReaderComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileReader_CFmeshReaderData_hh
