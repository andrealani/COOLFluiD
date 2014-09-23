#ifndef COOLFluiD_MutationUsage_ComputeWithMutation_hh
#define COOLFluiD_MutationUsage_ComputeWithMutation_hh

//////////////////////////////////////////////////////////////////////////////

#include "MutationI/MutationLibrary.hh"
#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace MutationUsage {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a base class to use Mutation for the
 * computation of thermodynamic quantities or trasport properties
 *
 * @author Andrea Lani
 *
 */
class ComputeWithMutation : public Common::OwnedObject {
public:
  
  typedef Environment::ConcreteProvider<ComputeWithMutation,1> PROVIDER;
  typedef Common::SafePtr<Physics::Mutation::MutationLibrary> ARG1;
    
  /**
   * Constructor
   */
  ComputeWithMutation
  (Common::SafePtr<Physics::Mutation::MutationLibrary> ptr) : 
    Common::OwnedObject(),
    _ptr(ptr)
  {
  }
  
  /**
   * Default destructor
   */
  virtual ~ComputeWithMutation() 
  {
  }
  
  /**
   * Compute the required quantities
   */
  virtual void compute() = 0;
  
  /**
   * Get the class name
   */
  static std::string getClassName() 
  {
    return "ComputeWithMutation";
  }
  
protected:
  
  /**
   * Get the library
   */
  Common::SafePtr<Physics::Mutation::MutationLibrary> getLibrary() const
  {
    return _ptr;
  }
  
private:
  
  /// pointer to the library
   Common::SafePtr<Physics::Mutation::MutationLibrary> _ptr;
  
}; // end of class ComputeWithMutation

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MutationUsage_ComputeWithMutation_hh
