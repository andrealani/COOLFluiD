// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_OutputFormatter_hh
#define COOLFluiD_Framework_OutputFormatter_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Framework/Method.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/MeshDataInputSource.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class SubSystemStatusStack;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a OutputFormatter.
/// This is an abstract class.
/// @author Tiago Quintino
class Framework_API OutputFormatter : public Method,
                    public Common::DynamicFunctionCaller<OutputFormatter> {

public: // typedefs

  typedef Environment::ConcreteProvider<OutputFormatter,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  OutputFormatter(const std::string& name);

  /// Default destructor.
  virtual ~OutputFormatter();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Opens the file for writing.
  /// @post pushs and pops the Namespace to which this Method belongs
  void open();

  /// Writes the solution in the Domain in the specified format.
  /// @post pushs and pops the Namespace to which this Method belongs
  void write();

  /// Closes the file
  /// @post pushs and pops the Namespace to which this Method belongs
  void close();

  /// Should the output formatter write now to file?
	/// @param force_write overrides the rate of saving
	/// @return if it is OK to save now
  bool isSaveNow( const bool force_write );

  /// Returns the extension of the files of this format
  /// @return std::string with extension
  virtual std::string getFormatExtension() const = 0;

  /// Sets the output file name
  void setOutputFileName(const std::string filename);
  
  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<OutputFormatter>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName() { return "OutputFormatter"; }

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

protected: // abstract interface implementations

  /// Opens the file for writing.
  /// This is the abstract function that the concrete methods must implement.
  virtual void openImpl() = 0;

  /// Writes the solution in the Domain in the specified format.
  /// This is the abstract function that the concrete methods must implement.
  virtual void writeImpl() = 0;

  /// Closes the file
  /// This is the abstract function that the concrete methods must implement.
  virtual void closeImpl() = 0;

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

protected: // helper methods

  /// Computes the m_fullOutputName and sets it
  virtual void computeFullOutputName();
  
protected: // member data

  /// Name of Solution File where to write
  boost::filesystem::path m_filename;

  /// Filename of the output file
  boost::filesystem::path m_fullOutputName;

  /// configuration string for solution file
  std::string m_filenameStr;

  /// Rate to save intermediate solution to solution file.
  CFuint  m_saveRate;

  /// If to save after the final iteration
  bool  m_savefinal;

  /// Save each iteration to different Tecplot file (with suffix m_iter#).
  bool  m_appendIter;

  /// Save each iteration to different Tecplot file (with suffix m_iter#).
  bool  m_appendTime;
  
  /// Append processor rank to file name
  bool  m_appendRank;
  
}; // class OutputFormatter

//////////////////////////////////////////////////////////////////////////////

  } // namespace OutputFormatter
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(OutputFormatter) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_OutputFormatter_hh
