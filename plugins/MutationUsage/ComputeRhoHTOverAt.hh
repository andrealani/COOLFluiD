#ifndef COOLFluiD_MutationUsage_ComputeRhoHTOverAt_hh
#define COOLFluiD_MutationUsage_ComputeRhoHTOverAt_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeWithMutation.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace MutationUsage {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the LTE composition for a chemically reactive 
 * mixture
 *
 * @author Andrea Lani
 *
 */
class ComputeRhoHTOverAt : public ComputeWithMutation {
public:
  
  /**
   * Constructor
   */
  ComputeRhoHTOverAt
  (Common::SafePtr<Physics::Mutation::MutationLibrary> ptr) : 
    ComputeWithMutation(ptr)
  {
  }
  
  /**
   * Default destructor
   */
  ~ComputeRhoHTOverAt() 
  {
  }
  
  /**
   * Compute the required quantities
   */
  virtual void compute();
    
}; // end of class ComputeRhoHTOverAt

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MutationUsage_ComputeRhoHTOverAt_hh
