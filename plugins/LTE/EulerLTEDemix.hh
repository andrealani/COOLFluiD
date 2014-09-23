#ifndef COOLFluiD_Physics_LTE_EulerLTEDemix_hh
#define COOLFluiD_Physics_LTE_EulerLTEDemix_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarTerm.hh"
#include "Framework/ConvectionPM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a Euler LTE model with
 * demixing (variable elemental fractions).&
 *
 * @author Andrea Lani
 *
 */

template <int DIM, class BASE>
class EulerLTEDemix : public Framework::ConvectionPM
<Framework::MultiScalarTerm<BASE> > {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  EulerLTEDemix(const std::string& name);

  /**
   * Default destructor
   */
  ~EulerLTEDemix();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("Euler" + Common::StringOps::to_str(DIM) + "DLTEDemix");
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
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * @return the number of elemental fractions
   */
  CFuint getNbElemFractions() const
  {
    return _nelem;
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
  virtual void setReferenceValues();

  /**
   * Set the reference value for time
   * This function is virtual to allow to use a different
   * adimensionalizaton, if needed
   */
  virtual void setReferenceTime();

private:

  /// number of elements
  CFuint _nelem;

}; // end of class EulerLTEDemix

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "EulerLTEDemix.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LTE_EulerLTEDemix_hh
