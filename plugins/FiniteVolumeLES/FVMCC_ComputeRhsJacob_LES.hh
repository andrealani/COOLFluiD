#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_LES_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_LES_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_ComputeRhsJacob.hh"
#include "LES/LESVarSet.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      class FVMCC_BC; 
      class FVMCC_ComputeRhsJacob;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS and Jacobian using
 * standard cell center FVM schemes for LES.
 * This adds volume to the diffusive flux calculation.
 *
 * @author Kris Van den Abeele
 * @author Willem Deconinck
 *
 */
class FVMCC_ComputeRhsJacob_LES : public FVMCC_ComputeRhsJacob {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacob_LES(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacob_LES();

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
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected:

  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();

  /**
   * Compute convective and diffusive fluxes
   */
  virtual void computeConvDiffFluxes(CFuint iVar, CFuint iCell);

protected:

  /// pointer to the LES varset
  Common::SafePtr<LES::LESVarSet> m_lesVarSet;

}; // class FVMCC_ComputeRhsJacob_LES

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacob_LESNumerics_hh
