// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_ParaViewWriter_ParaWriter_hh
#define COOLFluiD_IO_ParaViewWriter_ParaWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/OutputFormatter.hh"
#include "Common/NotImplementedException.hh"
#include "ParaWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a ParaView writer.
 *
 * @author Kris Van den Abeele
 */
class ParaWriter : public Framework::OutputFormatter {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit ParaWriter(const std::string& name);

  /**
   * Destructor
   */
  ~ParaWriter();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the extension of the files of this format
   * @return std::string with extension ".vtu"
   */
  std::string getFormatExtension() const;

  /**
   * Set the data to be written.
   * To make the writer free the internal data structure,
   * call this with a null pointer.
   */
  void bindDataImpl ()
  {
    // nothing to do here
  }

protected: // abstract interface implementations

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Opens the file for writing.
   * @see OutputFormatter::open()
   */
  virtual void openImpl();

  /**
   * Writes the solution in the Domain in the specified format.
   * @see OutputFormatter::write()
   */
  virtual void writeImpl();

  /**
   * Closes the file
   * @see OutputFormatter::close()
   */
  virtual void closeImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /**
   * Sets the data of the method.
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<ParaWriterCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<ParaWriterCom> m_unSetup;

  ///The writeSolution command to use
  Common::SelfRegistPtr<ParaWriterCom> m_writeSolution;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The writeSolution for configuration
  std::string m_writeSolutionStr;

  ///The data to share between ParaWriterCom commands
  Common::SharedPtr<ParaWriterData> m_data;

}; // end ParaWriter

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_ParaViewWriter_ParaWriter_hh
