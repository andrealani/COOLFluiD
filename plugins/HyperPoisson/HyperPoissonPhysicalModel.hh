#ifndef COOLFluiD_Physics_HyperPoisson_HyperPoissonPhysicalModel_hh
#define COOLFluiD_Physics_HyperPoisson_HyperPoissonPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "HyperPTerm.hh"
#include "Framework/ConvectionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a Hyperpolized Poisson Physical Model.
 *
 * @author Rayan Dhib
 * @author Andrea Lani
 */
template <int DIM>
class HyperPoissonPhysicalModel : public Framework::ConvectionPM<HyperPTerm> {
public:

  /**
   * Constructor without arguments
   */
  HyperPoissonPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~HyperPoissonPhysicalModel();

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

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
    return std::string("HyperPoisson" + Common::StringOps::to_str(DIM) + "D");
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
  virtual void setReferenceTime();

}; // end of class HyperPoissonPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

#include "HyperPoissonPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_HyperPoisson_HyperPoissonPhysicalModel_hh
