#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVTLTE_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVTLTE_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/SubInletEulerFunc.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with the initial conditions 
   * given for mass flow and temperature
   *
   * @author Radek Honzatko
   *
   */
class SubInletEuler2DPuvtUVTLTE : public SubInletEulerFunc {
public:

  /**
   * Constructor
   */
  SubInletEuler2DPuvtUVTLTE(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SubInletEuler2DPuvtUVTLTE();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
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

  /// physical model var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_varSet;
  
  /// physico-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// use the old (buggy) implementation
  bool m_useOld;
  
}; // end of class SubInletEuler2DPuvtUVTLTE

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVTLTE_hh
