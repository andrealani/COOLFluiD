#ifndef COOLFluiD_ShapeFunctions_SetQuadLagrangeP2SpectralFDP1StateCoord_hh
#define COOLFluiD_ShapeFunctions_SetQuadLagrangeP2SpectralFDP1StateCoord_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SetElementStateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }
  namespace Framework { class State; }

  namespace SpectralFD {


//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a functor and offers an abstract interface
 * for setting the corresponding space coordinates (Framework::Node) in
 * the State's in a quadrangle with P2 geometrical and P1 solution
 * interpolation.
 *
 * @author Kris Van den Abeele
 */
class SetQuadLagrangeP2SpectralFDP1StateCoord : public Framework::SetElementStateCoord {

public:

  /**
   * Constructor
   */
  SetQuadLagrangeP2SpectralFDP1StateCoord() : Framework::SetElementStateCoord()
  {
  }

  /**
   * Destructor
   */
  ~SetQuadLagrangeP2SpectralFDP1StateCoord()
  {
  }

  /**
   * Overloading of the operator () to make this class act as a
   * functor
   * @param nodes   list of the nodes in the current element
   * @param states  list of the states in the current element
   */
  void operator() (const std::vector<Framework::Node*>& nodes,
                   std::vector<Framework::State*>& states);

  /**
   * Function allowing to update the StateCoord
   * @param nodes   list of the nodes in the current element
   * @param states  list of the states in the current element
   */
  void update(const std::vector<Framework::Node*>& nodes,
                            std::vector<Framework::State*>& states);

}; // end of class SetQuadLagrangeP2SpectralFDP1StateCoord

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_SetQuadLagrangeP2SpectralFDP1StateCoord_hh
