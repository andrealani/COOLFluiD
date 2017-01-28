#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_BC_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_BC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/FVMCC_BaseBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
/**
 * This class represents the base class for any BC for cell center FV
 *
 * @author Andrea Lani
 *
 */
class FVMCC_BC : public Framework::FVMCC_BaseBC<CellCenterFVMCom> {

public: // functions

  /**
   * Constructor
   */
  FVMCC_BC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FVMCC_BC();

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face) = 0;
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Transfer the gradients data
   */
  virtual void transferGradientsData() {}
  
};
      
////////////////////////////////////////////////////////////////////////////// 


 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_BC_hh
