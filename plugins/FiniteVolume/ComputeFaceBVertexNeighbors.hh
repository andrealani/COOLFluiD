#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeFaceBVertexNeighbors_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeFaceBVertexNeighbors_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeStencil.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a functor that creates the set of face+vertex
 * neighbors (distance-1 neighbors) of each State
 *
 * @author Andrea Lani
 */
class ComputeFaceBVertexNeighbors : public ComputeStencil {
public:

  /**
   * Constructor
   */
  ComputeFaceBVertexNeighbors(const std::string& name);

  /**
   * Destructor
   */
  ~ComputeFaceBVertexNeighbors();

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() ();

 private: // helper functions

  /**
   * Create a list with all the boundary nodes
   */
  void createBoundaryNodesList(std::vector<CFuint>& bNodes) const;

}; // end of class ComputeFaceBVertexNeighbors

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeFaceBVertexNeighbors_hh
