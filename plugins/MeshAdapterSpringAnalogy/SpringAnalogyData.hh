#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_SpringAnalogyData_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_SpringAnalogyData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/MeshAdapterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

  namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Data Object that is accessed by the different
   * MeshSpringAnalogyCom 's that compose the MeshAdapterSpringAnalogy.
   *
   * @see SpringAnalogyCom
   *
   * @author Thomas Wuilbaut
   */
class SpringAnalogyData : public Framework::MeshAdapterData {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  SpringAnalogyData(Common::SafePtr<Framework::Method> owner);

  /**
   * Default destructor
   */
  ~SpringAnalogyData();

  /**
   * Configure the data from the supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "SpringAnalogy";
  }

  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> >
  getGeoWithNodesBuilder()
  {
    return &_geoWithNodesBuilder;
  }

  CFuint getTotalNbSteps()
  {
    return _totalNbSteps;
  }

  CFuint getCurrentStep()
  {
    return _currentStep;
  }

  void resetCurrentStep()
  {
    _currentStep = 1;
  }

  void updateCurrentStep()
  {
    _currentStep++;
  }

private:

  // builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _geoWithNodesBuilder;

  CFuint _totalNbSteps;

  CFuint _currentStep;

}; // end of class SpringAnalogyData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for MeshAdapterSpringAnalogy
typedef Framework::MethodCommand<SpringAnalogyData> SpringAnalogyCom;

/// Definition of a command provider for MeshAdapterSpringAnalogy
typedef Framework::MethodCommand<SpringAnalogyData>::PROVIDER SpringAnalogyComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_SpringAnalogyData_hh
