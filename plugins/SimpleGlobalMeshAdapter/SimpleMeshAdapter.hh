// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleMeshAdapter_hh
#define COOLFluiD_Numerics_SimpleGlobalMeshAdapter_SimpleMeshAdapter_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/MeshAdapterMethod.hh"
#include "SimpleMeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a MeshAdapterMethod that implements
 * a simple global mesh adaptation
 *
 * @author Thomas Wuilbaut
 *
 */
class SimpleMeshAdapter : public Framework::MeshAdapterMethod {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   *
   * @param name missing documentation
   */
  explicit SimpleMeshAdapter(const std::string& name);

  /**
   * Default destructor
   */
  ~SimpleMeshAdapter();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set the collaborating MeshCreator and OutputFormatter non-root
   */
  virtual void setNonRootMethods();

protected: // abstract interface implementations

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Adapts of the mesh
   * @see MeshAdapterMethod::adaptMesh()
   */
  virtual void adaptMeshImpl();

  /**
   * Remeshes
   * @see MeshAdapterMethod::remesh()
   */
  virtual void remeshImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /**
   * UnSets the data of the method.
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _unSetup;

  ///The Prepare commands to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _prepare;

  ///The Prepare commands to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _runMeshGenerator;

  ///The Prepare commands to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _readNewMesh;

  ///The Prepare commands to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _interpolateSol;

  ///The Prepare commands to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _writeInterpolatedMesh;

  ///The RemeshCondition commands to use
  Common::SelfRegistPtr<SimpleMeshAdapterCom> _remeshCondition;

  ///The Setup string for configuration
  std::string _setupStr;

  ///The prepare type for configuration
  std::string _prepareStr;

  ///The MeshAdapter command strings for configuration
  std::string _runMeshGeneratorStr;
  std::string _readNewMeshStr;
  std::string _interpolateSolStr;
  std::string _writeInterpolatedMeshStr;
  std::string _remeshConditionStr;
  std::string _unSetupStr;

  ///The data to share between SimpleGlobalMeshAdapterMethod commands
  Common::SharedPtr<SimpleMeshAdapterData> _data;

}; // class SimpleMeshAdapter

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SimpleGlobalMeshAdapter_hh
