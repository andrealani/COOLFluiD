#ifndef COOLFluiD_Physics_LinEuler_LinEulerPhysicalModel_hh
#define COOLFluiD_Physics_LinEuler_LinEulerPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinEulerTerm.hh"
#include "Framework/ConvectionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a LinEulerPhysicalModel.
 *
 *
 * @author Lilla Edit Koloszar
 */

template <int DIM>
class LinEulerPhysicalModel : public Framework::ConvectionPM<LinEulerTerm> {
public:

  /**
   * Constructor without arguments
   */
  LinEulerPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~LinEulerPhysicalModel();

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
    return std::string("LinEuler" + Common::StringOps::to_str(DIM) + "D");
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

}; // end of class LinEulerPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

#include "LinEulerPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_LinEulerPhysicalModel_hh
