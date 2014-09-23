#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobDiag_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobDiag_hh

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
class FVMCC_ComputeRhsJacobDiag : public FVMCC_ComputeRhsJacob {
public:
  
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRhsJacobDiag(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRhsJacobDiag();
  
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
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();

private:
  
  /**
   * Compute the jacobian contribution of the current (internal) face
   */
  virtual void computeBothJacobTerms();
  
  /**
   * Compute only the jacobian contribution to one of the two states 
   * (this is used in parallel computing)
   */
  virtual void computeJacobTerm(const CFuint idx);
  
  /**
   * Compute the jacobian contribution of the current (boundary) face
   */
  virtual void computeBoundaryJacobianTerm();
  
  /**
   * Compute the source term and its analytical jacobian
   * @param ist ID of the source term
   */
  virtual void computeSourceTermAnalytJacob(CFuint ist);
  
  /**
   * Compute the source term and its numerical jacobian
   * @param ist ID of the source term
   */
  virtual void computeSourceTermNumJacob(CFuint ist);	
  
protected:
  
  /// storage of diagonal block matrices 
  Framework::DataSocketSource<CFreal> socket_diagMatrices;
  
  /// storage of the local updatable IDs or -1 (ghost) for all local states
  Framework::DataSocketSource<CFint> socket_upLocalIDsAll;
            
  /// matrix iterator
  RealMatrix _matIter0;
  
  /// matrix iterator
  RealMatrix _matIter1;
    
}; // class FVMCC_ComputeRhsJacobDiag

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobDiag_hh
