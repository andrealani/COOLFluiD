#ifndef COOLFluiD_Physics_LSvki_ClarkPhysicalModel_hh
#define COOLFluiD_Physics_LSvki_ClarkPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTurbTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/ConvectionDiffusionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LESvki {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesPhysicalModel.
 *
 * @author Andrea Lani
 *
 */

template <int DIM>
 class ClarkPhysicalModel :
	public Framework::ConvectionDiffusionPM
	<NavierStokes::EulerTerm, NavierStokes::NSTurbTerm> {
 public:

   /**
    * Constructor without arguments
    */
  ClarkPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ClarkPhysicalModel();

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
    return std::string("Clark" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class SmagorinskyPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

#include "ClarkPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesPhysicalModel_hh
