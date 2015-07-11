#ifndef COOLFluiD_Physics_Poisson_PoissonModel_hh
#define COOLFluiD_Physics_Poisson_PoissonModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionDiffusionReactionPM.hh"
#include "Poisson/PoissonReacTerm.hh"
#include "Poisson/PoissonDiffTerm.hh"
#include "Poisson/PoissonConvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a Poisson
 * model.
 *
 * @author Alejandro Alvarez Laguna
 *
 */
template <int DIM>
class PoissonModel : public Framework::ConvectionDiffusionReactionPM
<PoissonConvTerm, PoissonDiffTerm, PoissonReacTerm> {
  
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  PoissonModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~PoissonModel();
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
    return std::string("Poisson" +  Common::StringOps::to_str(DIM) + "D");
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
// DATA
  
}; // end of class PoissonModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "PoissonModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Poisson_PoissonModel_hh
