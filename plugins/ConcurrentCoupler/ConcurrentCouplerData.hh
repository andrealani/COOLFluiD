#ifndef COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CouplerData.hh"

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
  
  /// @return the DataStorage corresponding to the given namespace
  Common::SafePtr<Framework::DataStorage> getDataStorage(const std::string& nspName);
  
  /// @return true if the current rank has to be involved in the data transfer
  bool isActiveRank(const std::vector<int>& flags) const
  {
    const int nspRank = Common::PE::GetPE().GetRank(getNamespace());
    cf_assert(nspRank < flags.size());
    return (flags[nspRank] == 1); 
  }
  
private:
  
  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;
  
  /// builder for face TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;
  
  /// handle to the space method
  Framework::MultiMethodHandle<Framework::SpaceMethod> _spaceMethod;
  
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
