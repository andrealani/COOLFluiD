#ifndef COOLFluiD_Numerics_FluctSplit_TwoLayerWeakBC2D_hh
#define COOLFluiD_Numerics_FluctSplit_TwoLayerWeakBC2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler2D
 *
 * @author Thomas Wuilbaut
 *
 */
class TwoLayerWeakBC2D : public WeakBC {

public:

  /**
   * Constructor.
   */
  TwoLayerWeakBC2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~TwoLayerWeakBC2D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Set the state vector in the ghost State's
   */
  virtual void setGhostState(const Framework::State& state,
                             Framework::State& gstate) = 0;

  /**
   * Compute the face normal in case of a moving boundary
   */
  void setMovingFaceNormal(const CFuint faceID);

protected:

  /// socket for the interRhs data
  Framework::DataSocketSink< CFreal> socket_interRhs;

  /// socket for the pastNormals data
  Framework::DataSocketSink< InwardNormalsData*> socket_pastNormals;

  /// socket for the interNormals data
  Framework::DataSocketSink< InwardNormalsData*> socket_interNormals;

  /// socket for the Nodes data
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for the past nodes data
  Framework::DataSocketSink< Framework::Node*> socket_pastNodes;

  /// socket for the past states data
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// socket for the intermidiate states data
  Framework::DataSocketSink< Framework::State*> socket_interStates;

  /// socket for the isUpdated data
  Framework::DataSocketSink< bool> socket_isUpdated;

  /// distribution coefficient
  CFreal                         m_alpha;

  /// SubSystem TimeStep
  CFreal                         m_dt;

  /// Index of the current layer
  CFuint                       m_layer;

  /// Speed of the face
  RealVector m_wallSpeed;

}; // end of class TwoLayerWeakBC2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_TwoLayerWeakBC2D_hh
