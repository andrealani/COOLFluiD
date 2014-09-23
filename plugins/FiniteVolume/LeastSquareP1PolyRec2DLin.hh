#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DLin_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DLin_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_PolyRecLin.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a least square polynomial reconstructor in 2D for FVM
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 */
class LeastSquareP1PolyRec2DLin : public FVMCC_PolyRecLin {
public:

  /**
   * Constructor
   */
  LeastSquareP1PolyRec2DLin(const std::string& name);

  /**
   * Default destructor
   */
  ~LeastSquareP1PolyRec2DLin();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

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

  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face);

  /**
   * Extrapolate the solution in the face quadrature points
   */
  void extrapolateImpl(Framework::GeometricEntity* const face,
		       CFuint iVar, CFuint leftOrRight);

private:

  /// socket for stencil
  Framework::DataSocketSink<
                            std::vector<Framework::State*> > socket_stencil;

  /// socket for weights
  Framework::DataSocketSink<CFreal> socket_weights;

  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;

  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;

  /// socket for Linearized State's
  Framework::DataSocketSink<Framework::State*> socket_linearizedStates;

  /// socket for Linearized Ghost State's
  Framework::DataSocketSink<Framework::State*> socket_linearizedGhostStates;

  RealVector  _l11;

  RealVector  _l12;

  RealVector  _l22;

  RealVector  _lf1;

  RealVector  _lf2;

}; // end of class LeastSquareP1PolyRec2DLin

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2DLin_hh
