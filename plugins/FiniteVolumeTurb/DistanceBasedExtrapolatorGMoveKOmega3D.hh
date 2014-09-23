#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveKOmega3D_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveKOmega3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/DistanceBasedExtrapolatorGMove.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesTurbVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used
 * in combination with BCs that moves the ghost states - extended to 3D
 *
 * @author Thomas Wuilbaut
 * @author Milan Zaloudek
 *
 */
class DistanceBasedExtrapolatorGMoveKOmega3D : public DistanceBasedExtrapolatorGMove {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveKOmega3D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveKOmega3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Extrapolate the solution in all mesh nodes
   */
  virtual void extrapolateInAllNodes();

  /**
   * Extrapolate the solution in the given nodes
   */
  virtual void extrapolateInNodes(const std::vector<Framework::Node*>& nodes);

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = DistanceBasedExtrapolatorGMove::needsSockets();

    result.push_back(&socket_wallDistance);

    return result;
  }

 private:

  /// socket for the wall distance values storage
  Framework::DataSocketSink<CFreal> socket_wallDistance;

  /// physical model convective variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesTurb3DVarSet> _diffVarTurb;

}; // end of class DistanceBasedExtrapolatorGMoveKOmega3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveKOmega3D_hh
