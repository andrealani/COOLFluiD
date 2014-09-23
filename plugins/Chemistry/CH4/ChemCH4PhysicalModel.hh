#ifndef COOLFluiD_Physics_Chemistry_CH4_ChemCH4PhysicalModel_hh
#define COOLFluiD_Physics_Chemistry_CH4_ChemCH4PhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalModelImpl.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a ChemCH4PhysicalModel.
 *
 * @author Tiago Quintino
 *
 */
class ChemCH4PhysicalModel : public Framework::PhysicalModelImpl {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  ChemCH4PhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~ChemCH4PhysicalModel();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("ChemCH4" + Common::StringOps::to_str(getDimension()) + "D");
  }

  /**
   * Get the convective name
   */
  std::string getConvectiveName() const
  {
    return getTypeName();
  }

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const
  {
    std::string result;
    result = getTypeName() + "Diffusive";
    return result;
  }

  /**
   * Get the Linear Source name
   */
  std::string getSourceName() const
  {
    return "Null";
  }

  /**
   * Set up
   */
  void setup();

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the convective term
   */
  virtual Common::SafePtr<Framework::BaseTerm> getConvectiveTerm() const
  {
    return CFNULL;
  }

  /**
   * Get the diffusive term
   */
  virtual Common::SafePtr<Framework::BaseTerm> getDiffusiveTerm() const
  {
    return CFNULL;
  }

  /**
   * Get the source term
   */
  virtual Common::SafePtr<Framework::BaseTerm> getSourceTerm() const
  {
    return CFNULL;
  }

  /**
   * Set the maximum nb of states data
   */
  void setMaxNbStatesData(const CFuint maxNbStatesData)
  {
    throw Common::NotImplementedException (FromHere(),"StructMechPhysicalModel::setMaxNbStatesData()");
  }

  /**
   * This struct groups all the data that is typical of this PhysicalModel
   */
  struct PhysicalData {

  }; // end PhysicalData

  /**
   * Get the Physical data
   */
  PhysicalData* getPhysicalData()
  {
    return &_physicalData;
  }

  /**
   * Get the R constant for the gas
   */
  CFreal getR() const
  {
    return _R;
  }

  /**
   * Get the Temperature of the gas
   */
  CFreal getTemperature() const
  {
    return _Temperature;
  }

  /**
   * Get the Pressure of the gas
   */
  CFreal getPressure() const
  {
    return _Pressure;
  }

protected:

  /// physical quantities to be computed in a linearized state
  PhysicalData _physicalData;

  /// constant temperature of the gas
  CFreal _Temperature;

  /// constant R of the gas
  CFreal _R;

  /// constant pressure of the gas
  CFreal _Pressure;

}; // end of class ChemCH4PhysicalModel

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Chemistry_CH4_ChemCH4PhysicalModel_hh
