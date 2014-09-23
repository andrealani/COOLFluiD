#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2DImpl_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2DImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakSlipWall2DImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an implicit weak moving slip wall bc (E. V. Weide) for Euler2D
 *
 * @author Thomas Wuilbaut
 *
 */
class UnsteadyWeakSlipWallEuler2DImpl : public WeakSlipWall2DImpl {
public:

  /**
   * Constructor.
   */
  UnsteadyWeakSlipWallEuler2DImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyWeakSlipWallEuler2DImpl();

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
  virtual void executeOnTrs();

  /**
   * Compute the normal flux and the corresponding jacobian
   */
  virtual void computeNormalFluxAndJacob(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux,
				 RealMatrix& fluxJacob);

  /**
   * Compute the face normal
   */
  void setFaceNormal(const CFuint faceID,
		     RealVector& normal);

protected:

  // the socket to the data handle of the state's
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  // the socket to the data handle of the past normals
  Framework::DataSocketSink<InwardNormalsData*> socket_pastNormals;

  // the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // the socket to the data handle of the pastNode's
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;

  /// Euler var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  /// Speed of the face
  RealVector _wallSpeed;

}; // end of class UnsteadyWeakSlipWallEuler2DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2DImpl_hh
