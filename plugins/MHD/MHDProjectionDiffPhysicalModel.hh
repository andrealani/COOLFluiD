#ifndef COOLFluiD_Physics_MHD_MHDProjectionDiffPhysicalModel_hh
#define COOLFluiD_Physics_MHD_MHDProjectionDiffPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHDProjectionTerm.hh"
#include "MHD/MHDProjectionDiffTerm.hh"
#include "Framework/ConvectionDiffusionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDProjectionDiffPhysicalModel.
 *
 * @author Andrea Lani
 *
 */

template <int DIM>
 class MHDProjectionDiffPhysicalModel :
	public Framework::ConvectionDiffusionPM
	<MHDProjectionTerm, MHDProjectionDiffTerm> {
 public:

   /**
    * Constructor without arguments
    */
  MHDProjectionDiffPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MHDProjectionDiffPhysicalModel();

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
    return std::string("MHD" + Common::StringOps::to_str(DIM) + "DProjectionDiff");
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

}; // end of class MHDProjectionDiffPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

#include "MHDProjectionDiffPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionDiffPhysicalModel_hh
