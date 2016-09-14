#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolator_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NodalStatesExtrapolator.hh"
#include "Common/ConnectivityTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object
 *
 * @author Andrea Lani
 *
 */
template <typename DATA> 
class DistanceBasedExtrapolator : public Framework::NodalStatesExtrapolator<DATA> {
public:
    
  /**
   * Constructor
   */
  DistanceBasedExtrapolator(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolator();

  /**
   * Set up private data needed by the computation
   */
  virtual void setup(); 

  /**
   * Unsetup private data needed by the computation
   */
  virtual void unsetup();
  
  /**
   * Extrapolate the solution in all mesh nodes
   */
  virtual void extrapolateInAllNodes();

  /**
   * Extrapolate the solution in the given nodes
   */
  virtual void extrapolateInNodes(const std::vector<Framework::Node*>& nodes);
  
protected:
  
  /**
   * Apply the nodal extrapolation
   */
  virtual void extrapolate();
  
  /**
   * Apply the nodal extrapolation in the inner part
   */
  void applyInner();
  
  /**
   * Update the weights
   */
  void updateWeights(CFuint nodeID, const Framework::Node& currNode);
   
protected:
  
  /// inverse of the weights
  /// it is worth storing this table in CFfloat (saves 50% of the memory
  /// than storing in CFdouble) because this is a potentially very
  /// big storage, depending on the grid complexity
  std::vector<CFreal*> _weights;
  
  // array of weights
  std::vector<CFreal> _weightsStorage;
  
}; // end of class DistanceBasedExtrapolator

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DistanceBasedExtrapolator.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolator_hh
