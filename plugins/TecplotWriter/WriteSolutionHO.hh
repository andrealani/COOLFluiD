// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_TecplotWriter_WriteSolutionHO_hh
#define COOLFluiD_IO_TecplotWriter_WriteSolutionHO_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FileWriter.hh"
#include "Framework/DataSocketSink.hh"

#include "TecplotWriter/TecWriterData.hh"
#include "TecplotWriter/MapGeoEntToTecplot.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class State; }

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents writes the MeshData to a Tecplot format file
/// for visualization.
/// This particular implementation subdivides the high-order elements
/// into linear elements
/// @author Tiago Quintino
class TecplotWriter_API WriteSolutionHO : public TecWriterCom,
                        public Framework::FileWriter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit WriteSolutionHO(const std::string& name);

  /// Virtual destructor
  virtual ~WriteSolutionHO();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup the object
  virtual void setup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

  /// Gets the file extension to append to the file name
  const std::string getWriterFileExtension() const
  {
    return std::string(".plt");
  }

protected: // functions

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

  /// Writes the Tecplot file header
  void write_tecplot_header(std::ofstream& fout);

protected: // data

  /// socket for State's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

private: // data

  /// File format to write the tecplot file in ASCII or binary form
  std::string m_fileFormatStr;
  /// dimensionality of the domain
  CFuint m_dimension;
  /// number of equations
  CFuint m_nbEqs;
  /// reference lenght
  CFreal m_refLenght;
  /// map geometric entities to Teplot
  MapGeoEntToTecplot m_mapgeoent;

}; // class WriteSolutionHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_WriteSolutionHO_hh

