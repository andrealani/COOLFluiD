#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeFaceVertexNeighborsPlusGhost_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeFaceVertexNeighborsPlusGhost_hh

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
class ComputeFaceVertexNeighborsPlusGhost : public ComputeStencil {
public:
  
  //   static std::pair<CFuint,CFuint> GLOBAL_STATE_ID;
  
  /**
   * Constructor
   */
  ComputeFaceVertexNeighborsPlusGhost(const std::string& name);

  /**
   * Destructor
   */
  ~ComputeFaceVertexNeighborsPlusGhost();

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() ();

}; // end of class ComputeFaceVertexNeighborsPlusGhost

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeFaceVertexNeighborsPlusGhost_hh
