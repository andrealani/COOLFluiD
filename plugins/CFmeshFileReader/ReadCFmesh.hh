// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_ReadCFmesh_hh
#define COOLFluiD_CFmeshFileReader_ReadCFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CFmeshFileReader.hh"

#include "CFmeshFileReader/ReadBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class is a MethodCommand that reads the data in the file.
/// This implementation loads in each processor the totality of the mesh
/// before partitioning, therefore cannot be used in testcases with big meshes.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Dries Kimpe
class CFmeshFileReader_API ReadCFmesh : public ReadBase {

public: // functions

  /// Constructor.
  explicit ReadCFmesh(const std::string& name);

  /// Destructor.
  virtual ~ReadCFmesh();

  /// Configures the command.
  void configure ( Config::ConfigArgs& args );

  /// Execute Processing actions
  void execute();

private: // data

  /// reader that will actually read the mesh
  Framework::CFmeshFileReader <Framework::CFmeshReaderSource> m_reader;

}; // class ReadCFmesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_ReadCFmesh_hh

