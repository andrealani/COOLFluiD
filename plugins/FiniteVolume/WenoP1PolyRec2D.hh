#ifndef COOLFluiD_Numerics_FiniteVolume_WenoP1PolyRec2D_hh
#define COOLFluiD_Numerics_FiniteVolume_WenoP1PolyRec2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "LeastSquareP1PolyRec2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a least square polynomial reconstructor in 2D for FVM
 *
 * @author Andrea Lani
 *
 */
class WenoP1PolyRec2D : public LeastSquareP1PolyRec2D {
public:

  /**
   * Constructor
   */
  WenoP1PolyRec2D(const std::string& name);

  /**
   * Default destructor
   */
  ~WenoP1PolyRec2D();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /**
   * Compute the gradients
   */
  void computeGradients();

  /**
   * Set up the private data
   */
  void setup();

  /**
   * Update the weights when nodes are moving
   */
  void updateWeights();

protected:
  
  /// socket for weights
  Framework::DataSocketSink<CFreal> socket_dr;
  
}; // end of class WenoP1PolyRec2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_WenoP1PolyRec2D_hh
