#ifndef COOLFluiD_Numerics_RemeshMeandros_RMeshMe_hh
#define COOLFluiD_Numerics_RemeshMeandros_RMeshMe_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/MeshAdapterMethod.hh"
#include "RMeshMeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

/// @author Jurek Majevski
class RMeshMe : public Framework::MeshAdapterMethod
{
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit RMeshMe(const std::string& name);

  /// Default destructor
  ~RMeshMe();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args arguments from where to read the configuration
  virtual void configure ( Config::ConfigArgs& args );

protected: // abstract interface implementations

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Adapts of the mesh
  /// @see MeshAdapterMethod::adaptMesh()
  virtual void adaptMeshImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<RMeshMeCom> _setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<RMeshMeCom> _unSetup;

  ///The command used for writing control space for meandros
  Common::SelfRegistPtr<RMeshMeCom> _writeControlSpcCom;

  ///The command used for calling meandros
  Common::SelfRegistPtr<RMeshMeCom> _meandrosCallCom;

  ///The command used used for interpolating solution to the new grid
  Common::SelfRegistPtr<RMeshMeCom> _interpolCom;

  ///The Setup string for configuration
  std::string _setupStr;

  ///The UnSetup string for configuration
  std::string _unSetupStr;

  ///The _writeControlSpc string for configuration
  std::string _writeControlSpcStr;

  ///The _writeMeandrosCall string for configuration
  std::string _meandrosCallStr;

  ///The _writeInterpol string for configuration
  std::string _interpolStr;


  ///The data to share between RemeshMeandrosMethod commands
  Common::SharedPtr<RMeshMeData> _data;

}; // class RMeshMe

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RemeshMeandros_hh
