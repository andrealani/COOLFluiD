// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_ParaViewWriter_WriteSolutionHighOrderHighOrder_hh
#define COOLFluiD_IO_ParaViewWriter_WriteSolutionHighOrder_hh

//////////////////////////////////////////////////////////////////////////////

#include "ParaWriterData.hh"
#include "Framework/FileWriter.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action
 * to write the MeshData solution to a ParaView format file
 * for visualization. Each high-order cell is written to a different piece,
 * (is in fact a miniature unstructured mesh)
 * This writer is suited for methods like discontinuous Galerkin, spectral volume,
 * spectral difference and related methods (the shape functions should be implemented though!!!)
 *
 * @author Kris Van den Abeele
 */
class WriteSolutionHighOrder : public ParaWriterCom,
                      public Framework::FileWriter {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit WriteSolutionHighOrder(const std::string& name);

  /**
   * Destructor.
   */
  ~WriteSolutionHighOrder()
  {
  }

  /**
    * Set up private data
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Gets the file extension to append to the file name
   */
  const std::string getWriterFileExtension() const
  {
    return std::string(".vtu");
  }

protected:

  /**
   * Write the ParaView file in Binmary format
   * @throw Common::FilesystemException
   */
  void writeToBinaryFile();

  /**
   * Write the to the given file stream the MeshData.
   * @throw Common::FilesystemException
   */
  void writeToFileStream(std::ofstream& fout);

  /**
   * Write the boundary surface data
   */
  void writeBoundarySurface();

  /**
   * Get the name of the writer
   */
  const std::string getWriterName() const;


  /**
   * function to check the endianness
   */
  bool isLittleEndian()
  {
    short int word = 0x0001;
    char *byte = (char *) &word;
    return byte[0];
  }

protected:

  // File format to write in (ASCII or Binary)
  std::string m_fileFormatStr;

}; // class WriteSolutionHighOrder

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_ParaViewWriter_WriteSolutionHighOrder_hh

