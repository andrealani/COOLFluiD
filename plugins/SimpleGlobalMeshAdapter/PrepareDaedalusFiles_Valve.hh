// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_PrepareDaedalusFiles_Valve_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_PrepareDaedalusFiles_Valve_hh

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
class PrepareDaedalusFiles_Valve : public SimpleMeshAdapterCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit PrepareDaedalusFiles_Valve(const std::string& name);

  /**
   * Destructor.
   */
  ~PrepareDaedalusFiles_Valve()
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

private:

  //Name of the mesh file
  std::string _script;

  //Name of the script start file
  std::string _scriptFile;

}; // class PrepareDaedalusFiles_Valve

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_PrepareDaedalusFiles_Valve_hh

