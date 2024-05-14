// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_IO_CGNSWriter_CGWriter_hh
#define COOLFluiD_IO_CGNSWriter_CGWriter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/OutputFormatter.hh"
#include "Common/NotImplementedException.hh"
#include "CGNSWriter/CGWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CGNSWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CGNS writer.
/// @author Rayan Dhib


class CGNSWriter_API CGWriter : public Framework::OutputFormatter {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit CGWriter(const std::string& name);

  /// Destructor
  ~CGWriter();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the extension of the files of this format
  /// @return std::string with extension ".plt"
  std::string getFormatExtension() const;
  
protected: // abstract interface implementations

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Opens the file for writing.
  /// @see OutputFormatter::open()
  virtual void openImpl();

  /// Writes the solution in the Domain in the specified format.
  /// @see OutputFormatter::write()
  virtual void writeImpl();

  /// Closes the file
  /// @see OutputFormatter::close()
  virtual void closeImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<CGWriterCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<CGWriterCom> m_unSetup;

  ///The writeSolution command to use
  Common::SelfRegistPtr<CGWriterCom> m_writeSolution;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The writeSolution for configuration
  std::string m_writeSolutionStr;

  ///The data to share between CGWriterCom commands
  Common::SharedPtr<CGWriterData> m_data;

}; // end CGWriter

//////////////////////////////////////////////////////////////////////////////

    } // namespace CGNSWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_CGNSWriter_CGWriter_hh
