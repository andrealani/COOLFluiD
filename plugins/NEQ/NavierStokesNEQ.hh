#ifndef COOLFluiD_Physics_NEQ_NavierStokesNEQ_hh
#define COOLFluiD_Physics_NEQ_NavierStokesNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokes Chemical NEQ
 * model.
 *
 * @author Andrea Lani
 * @author Janos Molnar
 *
 */
template <int DIM, typename CTERM, typename DTERM, typename RTERM>
class NavierStokesNEQ : public Framework::ConvectionDiffusionReactionPM<CTERM, DTERM, RTERM> {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  NavierStokesNEQ(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesNEQ();
   /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  virtual std::string getTypeName() const
  {
    return std::string("NavierStokes" +  Common::StringOps::to_str(DIM) + "D" + "NEQ");
  }

  /**
   * Get the convective name
   */
  virtual std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  virtual std::string getDiffusiveName() const;

  /**
   * Get the source name
   */
  virtual std::string getSourceName() const;

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;
  
  /**
   * @return the number of equations of the SubSystem
   */
  virtual CFuint getNbEquations() const;
  
protected:

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

  /**
   * Compute the reference chemistry time
   */
  virtual void computeTauChem();
  
protected:
  
  /// number of species
  CFuint _nbSpecies;

  /// number of Euler equations
  CFuint _nbEulerEqs;

  /// number of equations for vibrational energies
  CFuint _nbVibEnergyEqs;

  /// number of electron temperature
  CFuint _nbTe;

  /// free stream composition
  std::vector<CFreal> _yInf;

}; // end of class NavierStokesNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesNEQ.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_NavierStokesNEQ_hh
