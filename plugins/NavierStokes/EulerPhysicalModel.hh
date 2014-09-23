#ifndef COOLFluiD_Physics_NavierStokes_EulerPhysicalModel_hh
#define COOLFluiD_Physics_NavierStokes_EulerPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "EulerTerm.hh"
#include "Framework/ConvectionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a EulerPhysicalModel.
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
template <int DIM>
class EulerPhysicalModel : public Framework::ConvectionPM<EulerTerm> {
public:

  /**
   * Constructor without arguments
   */
  EulerPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~EulerPhysicalModel();

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
    return std::string("Euler" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class EulerPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

#include "EulerPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_EulerPhysicalModel_hh
