#ifndef COOLFluiD_Physics_NEQ_EulerNEQ_hh
#define COOLFluiD_Physics_NEQ_EulerNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/EulerTerm.hh"
#include "Framework/ConvectionReactionPM.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NEQ/NEQReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a Euler LTE model with
 * demixing (variable elemental fractions).&
 *
 * @author Andrea Lani
 *
 */

template <int DIM>
class EulerNEQ : public Framework::ConvectionReactionPM
<Framework::MultiScalarTerm<NavierStokes::EulerTerm>, NEQReactionTerm> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  EulerNEQ(const std::string& name);

  /**
   * Default destructor
   */
  ~EulerNEQ();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("Euler" + Common::StringOps::to_str(DIM) + "DNEQ");
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
    return "Null";
  }

  /**
   * Get the source name
   */
  std::string getSourceName() const
  {
    return getTypeName();
  }

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

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

  /// number of species
  CFuint _nbSpecies;

  /// number of Euler equations
  CFuint _nbEulerEqs;

  /// number of equations for vibrational energies
  CFuint _nbVibEnergyEqs;

  /// number of electron temperature
  CFuint _nbTe;

}; // end of class EulerNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "EulerNEQ.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_EulerNEQ_hh
