#ifndef COOLFluiD_Numerics_FiniteVolume_OpenBoundaryEasy3Fin_hh
#define COOLFluiD_Numerics_FiniteVolume_OpenBoundaryEasy3Fin_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MultiFluidMHD {
      class DiffMFMHD3DVarSet;
    }

    namespace MultiFluidMHD {
      template <class BASE> class MultiFluidMHDVarSet;
    }
    namespace Maxwell {
      class Maxwell3DProjectionVarSet;
    }
  }

 
  namespace Numerics {
    
    namespace FiniteVolume {
          
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a supersonic outlet command in 3D
   * for projection scheme
   * 
   * @author Michaela Brchnelova 
   *
   */
class OpenBoundaryEasy3Fin : public FVMCC_BC {
public: 

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  OpenBoundaryEasy3Fin(const std::string& name);
  
  /**
   * Default destructor
   */
  ~OpenBoundaryEasy3Fin();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object by complementing the 
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

private: //data

  /// physical model var set
  Common::SafePtr<Physics::MultiFluidMHD::MultiFluidMHDVarSet<Physics::Maxwell::Maxwell3DProjectionVarSet> > _updateVarSet;

  //Common::SafePtr<Physics::MHD::MHD3DProjection_TF_exp3VarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;
  CFreal _mptome;


}; // end of class OpenBoundaryEasy3Fin

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_OpenBoundaryEasy3Fin_hh
