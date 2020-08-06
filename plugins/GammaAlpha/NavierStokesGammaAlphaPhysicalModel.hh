#ifndef COOLFluiD_Physics_GammaAlpha_NavierStokesGammaAlphaPhysicalModel_hh
#define COOLFluiD_Physics_GammaAlpha_NavierStokesGammaAlphaPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTurbTerm.hh"
#include "GammaAlphaReactionTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionDiffusionReactionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesGammaAlphaPhysicalModel.
 *
 * @author Ray Vandenhoeck
 *
 */

template <int DIM>
class NavierStokesGammaAlphaPhysicalModel :
	public Framework::ConvectionDiffusionReactionPM
	<Framework::MultiScalarTerm<NavierStokes::EulerTerm>,
	 NavierStokes::NSTurbTerm, GammaAlphaReactionTerm> {
public:

   /**
    * Constructor without arguments
    */
  NavierStokesGammaAlphaPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesGammaAlphaPhysicalModel();

   /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure (Config::ConfigArgs& args);

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("NavierStokes" +  Common::StringOps::to_str(DIM) + "DGammaAlpha");
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

}; // end of class NavierStokesGammaAlphaPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace GammaAlpha

  } // end of namespace Physics

} // end of namespace COOLFluiD

#include "NavierStokesGammaAlphaPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_GammaAlpha_NavierStokesGammaAlphaPhysicalModel_hh
