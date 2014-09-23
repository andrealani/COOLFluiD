#ifndef COOLFluiD_Physics_Maxwell_MaxwellModel_hh
#define COOLFluiD_Physics_Maxwell_MaxwellModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionPM.hh"
#include "Maxwell/ConvMaxwellTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for Maxwell
 * model.
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 *
 */
template <int DIM>
class MaxwellModel : public Framework::ConvectionPM<ConvMaxwellTerm> {
  
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  MaxwellModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MaxwellModel();
   /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("Maxwell" +  Common::StringOps::to_str(DIM) + "D");
  }

  /**
   * Get the convective name
   */
   std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const;

  /**
   * Get the source name
   */
  std::string getSourceName() const;

  /**
   * @return the space dimension of the SubSystem
   */
  virtual CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  virtual CFuint getNbEquations() const;

 private:

  /**
   * Set the reference values
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceValues();

  /**
   * Set the reference value for time
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceTime();
  
private:
  
}; // end of class EulerNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MaxwellModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_MaxwellModel_hh
