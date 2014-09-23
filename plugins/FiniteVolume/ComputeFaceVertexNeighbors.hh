#ifndef COOLFluiD_Numerics_FiniteVolume_ComputeFaceVertexNeighbors_hh
#define COOLFluiD_Numerics_FiniteVolume_ComputeFaceVertexNeighbors_hh

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
class ComputeFaceVertexNeighbors : public ComputeStencil {
public:

  /**
   * Constructor
   */
  ComputeFaceVertexNeighbors(const std::string& name);

  /**
   * Destructor
   */
  ~ComputeFaceVertexNeighbors();

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   */
  void operator() ();

}; // end of class ComputeFaceVertexNeighbors

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_ComputeFaceVertexNeighbors_hh
