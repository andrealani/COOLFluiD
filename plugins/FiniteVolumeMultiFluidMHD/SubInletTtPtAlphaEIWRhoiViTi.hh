#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletTtPtAlphaEIWRhoiViTi_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletTtPtAlphaEIWRhoiViTi_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MultiFluidMHD {
      class DiffMFMHD2DVarSet;
    }

    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell2DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

    /**
     * This class represents a Subsonic Inlet imposing the Total Pressure and Total Temperature
     * Maxwell Equations: Perfectly Conducting Wall Condition
     *
     * @author Alejandro Alvarez
     *
     */
class SubInletTtPtAlphaEIWRhoiViTi : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletTtPtAlphaEIWRhoiViTi(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletTtPtAlphaEIWRhoiViTi();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

 private: // data

  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell2DProjectionVarSet> > _varSet;
  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

  /// total temperature
  CFreal     _tTotal;

  /// total pressure
  CFreal     _pTotal;

  /// alpha
  CFreal     _alpha;

}; // end of class SubInletTtPtAlphaEIWRhoiViTi

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletTtPtAlphaEIWRhoiViTi_hh
