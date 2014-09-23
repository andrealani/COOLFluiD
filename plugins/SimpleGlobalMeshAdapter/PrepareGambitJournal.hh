// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_PrepareGambitJournal_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_PrepareGambitJournal_hh

//////////////////////////////////////////////////////////////////////////////

#include "SimpleMeshAdapterData.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to generate a new mesh
   *
   * @author Thomas Wuilbaut
   *
   */
class PrepareGambitJournal : public SimpleMeshAdapterCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit PrepareGambitJournal(const std::string& name);

  /**
   * Destructor.
   */
  ~PrepareGambitJournal()
  {
  }

  /**
   * Configures this command
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Sets up this command
   */
  virtual void setup();

  /**
   * Execute Processing actions
   */
  void execute();

private:

  //Read old file and write new one
  void readWriteFile();

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);

private:

  //Name of the mesh file
  std::string _filename;

  //Name of the old journal file
  std::string _journalFile;

  //Name of the new journal file
  std::string _newJournalFile;

  /// storage for the current values of the vars
  RealVector _variables;

  /// storage for the current values of the gambit parameters
  RealVector _result;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the vars
  std::vector<std::string> _vars;

  /// a vector of string to hold the gambit parameters
  std::vector<std::string> _gambitParams;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // class PrepareGambitJournal

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_PrepareGambitJournal_hh

