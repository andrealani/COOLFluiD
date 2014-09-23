#ifndef COOLFluiD_Numerics_RemeshMeandros_RMeshMeData_hh
#define COOLFluiD_Numerics_RemeshMeandros_RMeshMeData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/MeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a Data Object that is accessed by the different
  /// RemeshMeandrosCom 's that compose the RemeshMeandros.
  /// @see RemeshMeandrosCom
  /// @author Jurek Majewski
class RMeshMeData : public Framework::MeshAdapterData {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  RMeshMeData(Common::SafePtr<Framework::Method> owner);

  /// Default destructor
  ~RMeshMeData();

  /// Configure the data from the supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "RMeshMe";
  }



  /// @return if the hessian should be inverted before smoothing
  const std::string& getMeandrosDir() const   	{ return _meandrosDirName; }



private:

  std::string  _meandrosDirName;



}; // end of class RMeshMeData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for RemeshMeandros
typedef Framework::MethodCommand<RMeshMeData> RMeshMeCom;  // TODO ::

/// Definition of a command provider for RemeshMeandros
typedef Framework::MethodCommand<RMeshMeData>::PROVIDER RMeshMeComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RemeshMeandros_RMeshMeData_hh
