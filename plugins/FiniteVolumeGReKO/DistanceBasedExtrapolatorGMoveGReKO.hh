#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveGReKO_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveGReKO_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/DistanceBasedExtrapolatorGMove.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used
 * in combination with BCs that moves the ghost states (nodes) like
 * @see NoSlipWallAdiabaticTurb2D
 *
 * @author Khalil Bensassi
 *
 */
template <class CVARSET, class DVARSET>
class DistanceBasedExtrapolatorGMoveGReKO : public DistanceBasedExtrapolatorGMove {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveGReKO(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveGReKO();

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
  Common::SafePtr<CVARSET> _varSetTurb;

  /// physical model diffusive variable set
  Common::SafePtr<DVARSET> _diffVarTurb;

}; // end of class DistanceBasedExtrapolatorGMoveGReKO

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "DistanceBasedExtrapolatorGMoveGReKO.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveGReKO_hh
