// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CFmeshFileReader_CFmeshReader_hh
#define COOLFluiD_IO_CFmeshFileReader_CFmeshReader_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MeshCreator.hh"
#include "Framework/MeshFormatConverter.hh"

#include "CFmeshFileReader/CFmeshReaderData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CFmeshReader.
/// @author Tiago Quintino
/// @author Andrea Lani
class CFmeshFileReader_API CFmeshReader : public Framework::MeshCreator {
public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit CFmeshReader(const std::string& name);

  /// Destructor
  ~CFmeshReader();
  
  /// Tell whether the mesh has been generated
  virtual bool isMeshGenerated() const {return !m_onlyConversion;}
  
  /// Configures the method, by allocating the it's dynamic members.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Modify the filename of the Mesh
  /// the new file is assumed to be a CFmesh file
  virtual void modifyFileNameForRestart(const std::string filename);

  /// Modify the filename of the Mesh
  virtual void modifyFileName(const std::string filename);

protected: // abstract interface implementations

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Generates the Mesh in the MeshData and the connectivity.
  /// @see MeshCreator::generateMeshData()
  virtual void generateMeshDataImpl();

  /// Process the CFMeshData to do for example: renumbering, FVM-FEM mesh conversion.
  /// @see MeshCreator::processMeshDataImpl()
  virtual void processMeshDataImpl();

  /// Builds the Mesh from the CFMeshData.
  /// @see MeshCreator::buildMeshData()
  virtual void buildMeshDataImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

private:  // helper functions

  /// Convert the mesh format to the CFmesh one
  void convertFormat();
  
  /// Helper function that actually does the job for converting
  void convert(Common::SelfRegistPtr<Framework::MeshFormatConverter> converter);
  
private: // data

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The Setup command to use
  Common::SelfRegistPtr<CFmeshReaderCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<CFmeshReaderCom> m_unSetup;

  /// string for configuring the readCFmesh command
  std::string m_readCFmeshStr;

  /// The ReadCFmesh command to use
  Common::SelfRegistPtr<CFmeshReaderCom> m_readCFmesh;

  /// The data to share between CFmeshReader commands
  Common::SharedPtr<CFmeshReaderData> m_data;

  /// string for configuring the converter
  std::string m_converterStr;

  /// option to convert back to the original format
  /// usefull only for debugging
  bool m_convertBack;
  
  /// option to choose to only convert the mesh without loading it into memory
  bool m_onlyConversion;
  
  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for delayed configuration of an object (MeshFormatConverter)
  Config::ConfigArgs m_stored_args;

}; // class CFmeshReader

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CFmeshFileReader_CFmeshReader_hh
