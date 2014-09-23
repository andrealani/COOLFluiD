#ifndef COOLFluiD_Numerics_FiniteVolume_VolumeBasedExtrapolator_hh
#define COOLFluiD_Numerics_FiniteVolume_VolumeBasedExtrapolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NodalStatesExtrapolator.hh"
#include "Common/ConnectivityTable.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object
 *
 * @author Thomas Wuilbaut
 *
 */
 class VolumeBasedExtrapolator : public Framework::NodalStatesExtrapolator<CellCenterFVMData> {
public:

  /**
   * Constructor
   */
  VolumeBasedExtrapolator(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~VolumeBasedExtrapolator();

  /**
   * Set up private data needed by the computation
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
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
      Framework::NodalStatesExtrapolator<CellCenterFVMData>::needsSockets();

    result.push_back(&socket_volumes);

    return result;
  }

protected: // helper function

  /**
   * Add the boundary neighbors to the list of the state neighbors
   * for each node
   */
  virtual void addBoundaryNeighbors
  (std::vector< std::vector<CFuint> >& bFacesPerNode,
   Common::CFMap<CFuint, CFint>& mapFaceTrsID);

protected:

  /// socket for cell volumes
  Framework::DataSocketSink< CFreal> socket_volumes;

  /// inverse of the weights
  /// it is worth storing this table in CFfloat (saves 50% of the memory
  /// than storing in CFdouble) because this is a potentially very
  /// big storage, depending on the grid complexity
  //   std::vector<float*> _weights;
  std::vector<CFreal*> _weights;

  // array of weights
  // std::vector<float> _weightsStorage;
  std::vector<CFreal> _weightsStorage;

  ///vector with the ID
  std::vector<CFuint> _ghost2InnerStates;

}; // end of class VolumeBasedExtrapolator

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_VolumeBasedExtrapolator_hh
