// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_ParaViewWriter_WriteSolution_hh
#define COOLFluiD_IO_ParaViewWriter_WriteSolution_hh

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
 * for visualization.
 *
 * @author Kris Van den Abeele
 */
class WriteSolution : public ParaWriterCom,
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
  explicit WriteSolution(const std::string& name);

  /**
   * Destructor.
   */
  ~WriteSolution()
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

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// socket for State Proxy
  Framework::DataSocketSink<
                            Framework::ProxyDofIterator<RealVector>*> socket_nstatesProxy;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
private:

  /// File format to write in (ASCII or Binary)
  std::string m_fileFormatStr;

}; // class WriteSolution

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_ParaViewWriter_WriteSolution_hh

