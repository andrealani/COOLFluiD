#ifndef COOLFluiD_ShapeFunctions_SetTriagLagrangeP1FluxReconstructionP5StateCoord_hh
#define COOLFluiD_ShapeFunctions_SetTriagLagrangeP1FluxReconstructionP5StateCoord_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SetElementStateCoord.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }
  namespace Framework { class State; }

  namespace FluxReconstructionMethod {


//////////////////////////////////////////////////////////////////////////////

/**
 * This class is a functor and offers an abstract interface
 * for setting the corresponding space coordinates (Framework::Node) in
 * the State's in a triangle with P1 geometrical and P5 solution
 * interpolation.
 *
 * @warning No coordinates are set, as the states are not really related to a specific points in space
 *
 * @author Rayan Dhib
 */
class SetTriagLagrangeP1FluxReconstructionP5StateCoord : public Framework::SetElementStateCoord {

public:

  /**
   * Constructor
   */
  SetTriagLagrangeP1FluxReconstructionP5StateCoord() : Framework::SetElementStateCoord()
  {
  }

  /**
   * Destructor
   */
  ~SetTriagLagrangeP1FluxReconstructionP5StateCoord()
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

}; // end of class SetTriagLagrangeP1FluxReconstructionP5StateCoord

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_SetTriagLagrangeP1FluxReconstructionP5StateCoord_hh
