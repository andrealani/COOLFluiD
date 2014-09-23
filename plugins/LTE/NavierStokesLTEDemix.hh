#ifndef COOLFluiD_Physics_LTE_NavierStokesLTEDemix_hh
#define COOLFluiD_Physics_LTE_NavierStokesLTEDemix_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTerm.hh"
//#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionDiffusionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokes LTE model
 * with elemental demixing.
 *
 * @author Andrea Lani
 *
 */
template <int DIM, class BASE>
class NavierStokesLTEDemix :
	public Framework::ConvectionDiffusionPM
	<Framework::MultiScalarTerm<BASE>,
	 NavierStokes::NSTerm> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  NavierStokesLTEDemix(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesLTEDemix();

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the convective name
   */
   std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const;

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

  /// number of elements
  CFuint _nelem;

}; // end of class NavierStokesLTEDemix

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesLTEDemix.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_NavierStokesLTEDemix_hh
