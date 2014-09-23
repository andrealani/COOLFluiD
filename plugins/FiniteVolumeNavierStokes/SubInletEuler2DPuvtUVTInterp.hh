#ifndef COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVTInterp_hh
#define COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVTInterp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OldLookupTable.hh"
#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a subsonic inlet command with
   * the initial conditions given for u(y),v(y),T(y) from file
   *
   * @author Andrea Lani
   *
   */
class SubInletEuler2DPuvtUVTInterp : public FVMCC_BC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SubInletEuler2DPuvtUVTInterp(const std::string& name);

  /**
   * Default destructor
   */
  ~SubInletEuler2DPuvtUVTInterp();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

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
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  /// look up table for u(y)
  Common::LookUpTable<CFreal,CFreal> _lookupTableU;

  /// look up table for v(y)
  Common::LookUpTable<CFreal,CFreal> _lookupTableV;

  /// look up table for T(y)
  Common::LookUpTable<CFreal,CFreal> _lookupTableT;

  /// input file name
  std::string _datafile;

}; // end of class SubInletEuler2DPuvtUVTInterp

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SubInletEuler2DPuvtUVTInterp_hh
