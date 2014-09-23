#ifndef COOLFluiD_MutationUsage_ComputeTransportCoefs_hh
#define COOLFluiD_MutationUsage_ComputeTransportCoefs_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeWithMutation.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace MutationUsage {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the transport coefficients for LTE and Demixing 
 * mixture
 *
 * @author Janos Molnar
 *
 */
class ComputeTransportCoefs : public ComputeWithMutation {
public:
  
  /**
   * Constructor
   */
  ComputeTransportCoefs
  (Common::SafePtr<Physics::Mutation::MutationLibrary> ptr) : 
    ComputeWithMutation(ptr)
  {
  }
  
  /**
   * Default destructor
   */
  ~ComputeTransportCoefs() 
  {
  }
  
  /**
   * Compute the required quantities
   */
  virtual void compute();
    
}; // end of class ComputeTransportCoefs

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MutationUsage_ComputeTransportCoefs_hh
