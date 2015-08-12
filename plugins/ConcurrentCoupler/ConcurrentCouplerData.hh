#ifndef COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CouplerData.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * ConcurrentCouplerCom commands that compose @see ConcurrentCoupler
   *
   * @author Andrea Lani
   * @author Thomas Wuilbaut
   */
class ConcurrentCouplerData : public Framework::CouplerData {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Default constructor without arguments
   */
  ConcurrentCouplerData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~ConcurrentCouplerData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Sets up the FiniteElementData
   */
  void setup();

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> >
  getStdTrsGeoBuilder()
  {
    return &_stdTrsGeoBuilder;
  }

  /**
   * @return the GeometricEntity builder for faces in FVMCC
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
  getFaceTrsGeoBuilder()
  {
    return &_faceTrsGeoBuilder;
  }
  
  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ConcurrentCouplerData";
  }

  /**
   * Sets the SpaceMethod which this Coupler uses
   */
  void setSpaceMethod(Framework::MultiMethodHandle<Framework::SpaceMethod> spaceMtd)
  {
    _spaceMethod= spaceMtd;
  }
  
  /**
   * Gets the SpaceMethod which this DataProcessing uses
   * @return pointer to the SpaceMethod
   */
  Framework::MultiMethodHandle<Framework::SpaceMethod> getSpaceMethod() const
  {
    return _spaceMethod;
  }
  
  /**
   * Get the transformer from send (source) to recv (target) variables
   */
  Common::SafePtr<Framework::VarSetTransformer> getSendToRecvVecTrans() const
  {
    cf_assert(_sendToRecvVecTrans.isNotNull());
    return _sendToRecvVecTrans.getPtr();
  }
  
private:
  
  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;
  
  /// builder for face TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;
  
  /// handle to the space method
  Framework::MultiMethodHandle<Framework::SpaceMethod> _spaceMethod;
  
  /// vector transformer from send (source) to recv (target) variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _sendToRecvVecTrans;
  
  /// Name of the vector transformer from send (source) to recv (target) variables
  std::string _sendToRecvVecTransStr;
  
}; // end of class ConcurrentCouplerData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Coupler
typedef Framework::MethodCommand<ConcurrentCouplerData> ConcurrentCouplerCom;

/// Definition of a command provider for Coupler
typedef Framework::MethodCommand<ConcurrentCouplerData>::PROVIDER ConcurrentCouplerComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh
