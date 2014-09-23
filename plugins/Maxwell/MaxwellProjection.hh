#ifndef COOLFluiD_Physics_Maxwell_MaxwellProjection_hh
#define COOLFluiD_Physics_Maxwell_MaxwellProjection_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionPM.hh"
#include "MaxwellProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MaxwellProjection.
 *
 * @author Alejandro Alvarez
 * @author Andrea Lani
 *
 */

template <int DIM>
class MaxwellProjection : public Framework::ConvectionPM<MaxwellProjectionTerm> {

public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  MaxwellProjection(const std::string& name);

  /**
   * Default destructor
   */
  ~MaxwellProjection();

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;
  
   /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );  

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("Maxwell" + Common::StringOps::to_str(DIM) + "DProjection");
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
  virtual void setReferenceTime()
  {
    /// @warning the reference time set below is only okay if the reference length is equal
    /// to one, or if the diffusion coefficient is equal to zero. To make it general, the diffusion
    /// coefficient should be rescaled also.
    _refTime = getRefLength();
  }
}; // end of class MaxwellProjection

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

#include "MaxwellProjection.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_MaxwellProjection_hh
