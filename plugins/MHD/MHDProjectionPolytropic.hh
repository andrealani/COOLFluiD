#ifndef COOLFluiD_Physics_MHD_MHDProjectionPolytropic_hh
#define COOLFluiD_Physics_MHD_MHDProjectionPolytropic_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionPM.hh"
#include "MHDProjectionPolytropicTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a MHDProjectionPolytropic.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */

template <int DIM>
class MHDProjectionPolytropic : public Framework::ConvectionPM<MHDProjectionPolytropicTerm> {

public:

  /**
   * Constructor without arguments
   */
  MHDProjectionPolytropic(const std::string& name);

  /**
   * Default destructor
   */
  ~MHDProjectionPolytropic();

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

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
    return std::string("MHD" + Common::StringOps::to_str(DIM) + "DProjectionPolytropic");
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
  virtual void setReferenceValues()
  {
  }

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
}; // end of class MHDProjectionPolytropic

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

#include "MHDProjectionPolytropic.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MHD_MHDProjectionPolytropic_hh
