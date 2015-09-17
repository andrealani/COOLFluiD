#ifndef COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHDModel_hh
#define COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHDModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionDiffusionReactionPM.hh"
#include "Framework/MultiScalarTerm.hh"
#include "MultiFluidMHD/ReacMFMHDTerm.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "MultiFluidMHD/DiffMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a multi-fluid MHD 
 * model.
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 *
 */
template <int DIM>
class MultiFluidMHDModel : public Framework::ConvectionDiffusionReactionPM
<Framework::MultiScalarTerm<EulerMFMHDTerm>, DiffMFMHDTerm, ReacMFMHDTerm> {
  
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  MultiFluidMHDModel(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MultiFluidMHDModel();
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
    return (this->is2DHalf()) ? 
      std::string("MultiFluidMHD" +  Common::StringOps::to_str(DIM) + "DHalf") : 
      std::string("MultiFluidMHD" +  Common::StringOps::to_str(DIM) + "D");
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
  
  /**
   * @return the number of species
   */
  virtual CFuint getNbSpecies() const {return _nbSpecies;}  
  

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
  
}; // end of class MultiFluidMHDModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "MultiFluidMHDModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHDModel_hh
