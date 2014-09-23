#ifndef COOLFluiD_MutationUsage_ComputeLTEComposition_hh
#define COOLFluiD_MutationUsage_ComputeLTEComposition_hh

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
class ComputeLTEComposition : public ComputeWithMutation {
public:
  
  /**
   * Constructor
   */
  ComputeLTEComposition
  (Common::SafePtr<Physics::Mutation::MutationLibrary> ptr) : 
    ComputeWithMutation(ptr)
  {
  }
  
  /**
   * Default destructor
   */
  ~ComputeLTEComposition() 
  {
  }
  
  /**
   * Compute the required quantities
   */
  virtual void compute();
    
}; // end of class ComputeLTEComposition

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MutationUsage_ComputeLTEComposition_hh
