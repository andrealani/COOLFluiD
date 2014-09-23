#ifndef COOLFluiD_Physics_NavierStokes_NavierStokesPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokesPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NSTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/ConvectionDiffusionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a NavierStokesPhysicalModel.
 *
 * @author Andrea Lani
 *
 */

template <int DIM>
 class NavierStokesPhysicalModel :
	public Framework::ConvectionDiffusionPM
	<EulerTerm, NSTerm> {
 public:

   /**
    * Constructor without arguments
    */
  NavierStokesPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesPhysicalModel();

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
    return std::string("NavierStokes" + Common::StringOps::to_str(DIM) + "D");
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
   
   /// Set the reference values for the compressible case
  void setReferenceValuesCompressible();
   
   /// Set the reference values for the incompressible case
  void setReferenceValuesIncompressible();
   
   /**
   * Set the reference value for time
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceTime();

}; // end of class NavierStokesPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

#include "NavierStokesPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesPhysicalModel_hh
