#ifndef COOLFluiD_Numerics_MeshFEMMove_FEMMoveData_hh
#define COOLFluiD_Numerics_MeshFEMMove_FEMMoveData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "MathTools/RealVector.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/MeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * MeshFEMMoveCom 's that compose the MeshFEMMove.
   *
   * @see MeshFEMMoveCom
   *
   * @author Thomas Wuilbaut
   */
class FEMMoveData : public Framework::MeshAdapterData {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  FEMMoveData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~FEMMoveData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "FEMMoveData";
  }

  /**
   * Sets the collaborator convergence methods into the data
   */
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod>& convergenceMtd)
  {
    _convergenceMtd = convergenceMtd;
  }

  /**
   * Gets the convergence method
   */
  Framework::MultiMethodHandle<Framework::ConvergenceMethod>& getConvergenceMethod()
  {
    return _convergenceMtd;
  }

  /**
   * Gets the convergence method
   */
  std::string getOtherNamespace()
  {
    return _otherNamespace;
  }

private: //member data

  /// ConvergenceMethod used to discretize the domain
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> _convergenceMtd;

  /// The namespace of the subsystem
  std::string _otherNamespace;

}; // end of class FEMMoveData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for MeshFEMMove
typedef Framework::MethodCommand<FEMMoveData> FEMMoveCom;

/// Definition of a command provider for MeshFEMMove
typedef Framework::MethodCommand<FEMMoveData>::PROVIDER FEMMoveComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_FEMMoveData_hh
