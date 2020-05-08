#ifndef COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaPhysicalModel_hh
#define COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTurbTerm.hh"
#include "KOmegaReactionTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesKLogOmegaPhysicalModel.
 *
 * @author Ray Vandenhoeck
 *
 */

template <int DIM>
class NavierStokesKLogOmegaPhysicalModel :
	public Framework::ConvectionDiffusionReactionPM
	<Framework::MultiScalarTerm<NavierStokes::EulerTerm>,
	 NavierStokes::NSTurbTerm, KOmegaReactionTerm<Framework::BaseTerm> > {
public:

   /**
    * Constructor without arguments
    */
  NavierStokesKLogOmegaPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesKLogOmegaPhysicalModel();

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
    return std::string("NavierStokes" +  Common::StringOps::to_str(DIM) + "DKLogOmega");
  }

  /**
   * Get the convective name
   */
  virtual  std::string getConvectiveName() const;

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const
  {
    return getTypeName();
  }

  /**
   * Get the source name
   */
  virtual std::string getSourceName() const
  {
    return getTypeName();
  }

  /**
   * @return the space dimension of the SubSystem
   */
  virtual CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  virtual CFuint getNbEquations() const;

  /**
   * Set the reference values
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceValues();

}; // end of class NavierStokesKLogOmegaPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

#include "NavierStokesKLogOmegaPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_KOmega_NavierStokesKLogOmegaPhysicalModel_hh
