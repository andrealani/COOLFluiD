// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_TecplotWriter_WriteSolutionHighOrder_hh
#define COOLFluiD_IO_TecplotWriter_WriteSolutionHighOrder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FileWriter.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

#include "TecplotWriter/TecWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action
/// to write the MeshData solution to a Tecplot format file
/// for visualization. Each high-order cell is written to a different zone,
/// (is in fact a miniature unstructured mesh)
/// This writer is suited for methods like discontinuous Galerkin, spectral volume,
/// spectral difference and related methods (the shape functions should be implemented though!!!)
/// @author Kris Van den Abeele
class TecplotWriter_API WriteSolutionHighOrder : public TecWriterCom,
                      public Framework::FileWriter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit WriteSolutionHighOrder(const std::string& name);

  /// Destructor.
  ~WriteSolutionHighOrder()
  {
  }

    /// Set up private data
  void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

  /// Gets the file extension to append to the file name
  const std::string getWriterFileExtension() const
  {
    return std::string(".plt");
  }

protected:

  /// Write the Tecplot file in Binmary format
  /// @throw Common::FilesystemException
  void writeToBinaryFile();

  /// Write the to the given file stream the MeshData.
  /// @throw Common::FilesystemException
  void writeToFileStream(std::ofstream& fout);

  /// Write the boundary surface data
  void writeBoundarySurface();

  /// Get the name of the writer
  const std::string getWriterName() const;

  /// returns the Tecplot cell shape
  static std::string getTecplotCellShape(CFGeoShape::Type shape,CFuint geoOrder);

  /// returns the mapped coordinates of the output points in a cell of given shape and order
  static std::vector< RealVector > getOutputPntsMappedCoords(CFGeoShape::Type shape,CFuint solOrder);

  /// returns the local subcell node connectivity in a cell of given shape and order
  static std::vector< std::vector< CFuint > > getOutputCellNodeConn(CFGeoShape::Type shape,CFuint solOrder);
  
private:

  /// File format to write in (ASCII or Binary)
  std::string _fileFormatStr;

  /// Builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >  m_stdTrsGeoBuilder;

}; // class WriteSolutionHighOrder

//helper functions

  /// Compute the midpoints of all 6 edges of a tetrahedron
  std::vector<RealVector> compute_midpoints(std::vector<RealVector> tetrahedron);
  // Subdivide a tetrahedron into 8 sub-tetrahedra
  std::vector<std::vector<RealVector> > subdivide_tetrahedron(std::vector<RealVector> tetrahedron, std::vector<RealVector> midpoints);
//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_WriteSolutionHighOrder_hh

