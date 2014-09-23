// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_ParReadCFmesh_hh
#define COOLFluiD_CFmeshFileReader_ParReadCFmesh_hh

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileReader/ReadBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

/// This class is a MethodCommand that reads the data in the file.
/// This implementation uses properly the MPI bindings to only load
/// in each processor a part of the mesh.
/// @author Andrea Lani
/// @author Tiago Quintino
/// @author Dries Kimpe
template <typename READER>
class CFmeshFileReader_API ParReadCFmesh : public ReadBase {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ParReadCFmesh(const std::string& name);

  /// Destructor.
  virtual ~ParReadCFmesh();

  /// Configures the command.
  void configure ( Config::ConfigArgs& args );

  /// Execute Processing actions
  void execute();

private: // data

  /// CFmesh file reader
  READER m_reader;

  /// user option to renumber the states
  bool m_renumber;

}; // class ParReadCFmesh

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "CFmeshFileReader/ParReadCFmesh.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_ParReadCFmesh_hh
