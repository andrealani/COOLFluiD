#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEulerMassFlowT_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEulerMassFlowT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/SubInletEulerFunc.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class EulerTerm;
    }
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;		
    template <typename BASE> class MultiScalarTerm;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions 
   * given for mass flow and temperature
   *
   * @author Andrea Lani
   */
class SubInletEulerMassFlowT : public SubInletEulerFunc {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  SubInletEulerMassFlowT(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubInletEulerMassFlowT();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face);

private: // data

  /// pointer o the Euler term
  Common::SafePtr<Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> > m_term;
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// flag telling if the state has partial densities
  bool m_stateHasPartialDensities;
  
  /// molar masses
  RealVector m_masses;
  
  /// inlet values for y_i
  std::vector<CFreal> m_yi;
  
}; // end of class SubInletEulerMassFlowT

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEulerMassFlowT_hh
