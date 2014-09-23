#ifndef COOLFluiD_Numerics_FiniteVolume_FarField2DYiPuvt_hh
#define COOLFluiD_Numerics_FiniteVolume_FarField2DYiPuvt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  

  namespace Physics {
    namespace NavierStokes {
      template <class BASEVS> class MultiScalarVarSet;
      class Euler2DVarSet;
    }
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic outlet command
 *
 * @author Khalil Bensassi
 *
 */
class FarField2DYiPuvt : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FarField2DYiPuvt(const std::string& name);

  /**
   * Default destructor
   */
   ~FarField2DYiPuvt();
  //virtual ~SubOutletEuler2DYiPuvt();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  //virtual void setup();
   void setup();

  /**
   * Apply boundary condition on the given face
   */
  //virtual void setGhostState(Framework::GeometricEntity *const face);
   void setGhostState(Framework::GeometricEntity *const face);

 private: // data
  


 /// physical model var set
  Common::SafePtr<Physics::NavierStokes::MultiScalarVarSet<Physics::NavierStokes::Euler2DVarSet> > _varSet;
  

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;


  // static pressure or pressure variation
  CFreal m_pressure;

  /// number of dimensions
  //CFuint m_dim;

 /// physical model data
  RealVector _dataInnerState;

 // molar masses
  RealVector _mm;
  RealVector _Yi;

  /// Yi mass fraction
  //std::vector<CFreal> _Yi;  

 /// Vibrational temperatures
  std::vector<CFreal>  _Tv;

  /// inlet values to be fixed [V T Tv]
  /// the partial pressures will be extrapolated from inside 
  std::vector<CFreal> m_vTTv;

  /// mass species fractions in the internal cell
  RealVector m_ysIn;
  
  /// physical data array
  RealVector m_pData;
}; // end of class FarFieldEuler2DYiPuvt

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FarFieldEuler2DYiPuvt_hh
