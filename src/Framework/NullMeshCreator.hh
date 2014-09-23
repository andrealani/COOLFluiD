// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullMeshCreator_hh
#define COOLFluiD_Framework_NullMeshCreator_hh

//////////////////////////////////////////////////////////////////////////////

#include "MeshCreator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullMeshCreator.
/// @author Tiago Quintino
class Framework_API NullMeshCreator : public MeshCreator {
public:

  /// Constructor.
  explicit NullMeshCreator(const std::string& name);

  /// Destructor
  ~NullMeshCreator();

  /// Modify the filename of the Mesh
  /// the new file is assumed to be a CFmesh file
  void modifyFileNameForRestart(const std::string filename);

  /// Modify the filename of the Mesh
  void modifyFileName(const std::string filename);

  /// Checks if this object is a Null object.
  /// Since this is NullMeshCreator
  /// @return true
  virtual bool isNull() const { return true; }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Generates the Mesh in the MeshData and the connectivity.
  /// @note since this is a Null Method it doesn't do anything but
  /// issue a warning.
  /// @see MeshCreator::generateMeshData()
  virtual void generateMeshDataImpl();

  /// Builds the Mesh from the CFMeshData.
  /// @note since this is a Null Method it doesn't do anything but
  /// issue a warning.
  /// @see MeshCreator::buildMeshData()
  virtual void buildMeshDataImpl();

  /// Process the CFMeshData to do for example: renumbering, FVM-FEM mesh conversion.
  /// @note since this is a Null Method it doesn't do anything but
  /// issue a warning.
  /// @see MeshCreator::processMeshDataImpl()
  virtual void processMeshDataImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

}; // class NullMeshCreator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullMeshCreator_hh
