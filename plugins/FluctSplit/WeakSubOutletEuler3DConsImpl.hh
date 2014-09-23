#ifndef COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler3DConsImpl_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler3DConsImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC3DImpl.hh"
#include "NavierStokes/Euler3DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a weak sub inlet bc for Euler3D
 *
 * @author Andrea Lani
 * @author Tiago Quintino
 *
 */
class WeakSubOutletEuler3DConsImpl : public WeakBC3DImpl {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSubOutletEuler3DConsImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSubOutletEuler3DConsImpl();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );


 protected:

  /**
   * Set the additional flux and the jacobian of the fluxes
   */
  void computeFluxAndJacob(std::vector<Framework::State*>& states,
         RealVector& flux,
         RealMatrix& fluxJacob);

private:

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler3DCons> m_varSet;

  /// temporary storage for the jacobian of the gstate
  RealMatrix m_dUoutDu;

  /// storage for the identity matrix (only the diagonal is stored)
  RealVector m_identity;

  /// total temperature
  CFreal m_pressure;

}; // end of class WeakSubOutletEuler3DConsImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSubOutletEuler3DConsImpl_hh
