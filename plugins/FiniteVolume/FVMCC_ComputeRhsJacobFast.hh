#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFast_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFast_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Numerics {

    namespace FiniteVolume {
  
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_ComputeRhsJacobFast : public FVMCC_ComputeRhsJacob {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobFast(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_ComputeRhsJacobFast();

protected:
  
  /// Restore the backed up left states
  virtual void restoreState(CFuint iCell)
  {
    _polyRec->getExtrapolatedPhysicaData()[iCell] = _polyRec->getBackupPhysicaData()[iCell];
  }
  
  /**
   * Compute the physical data in the states
   */
  virtual void computePhysicalData()
  {
    computeAndBackUpStatesData();
  }
  
  /// Compute the left perturbed states data
  void computePerturbedStatesData(CFuint iVar, CFuint iCell)
  {
    _reconstrVar->computePerturbedPhysicalData(*_polyRec->getExtrapolatedValues()[iCell], 
					       _polyRec->getBackupPhysicaData()[iCell],
					       _polyRec->getExtrapolatedPhysicaData()[iCell], iVar);
  }
  
  /// Compute convective and diffusive fluxes
  virtual void computeConvDiffFluxes(CFuint iVar, CFuint iCell);
  
}; // class FVMCC_ComputeRhsJacobFast

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobFast_hh
