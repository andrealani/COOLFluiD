#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2D_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/WeakSlipWallEuler2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
    
/**
 * This class represents a strong slip wall bc for Euler2D, indipendent from
 * the (update) variables that are used to store the solution
 *
 * @author Thomas Wuilbaut
 * @author Andrea Lani
 */
class UnsteadyWeakSlipWallEuler2D : public WeakSlipWallEuler2D {
public:

  /**
   * tructor.
   */
  UnsteadyWeakSlipWallEuler2D(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~UnsteadyWeakSlipWallEuler2D();

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
   * Compute the normal flux
   */
  virtual void computeNormalFlux(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux);
  
  /**
   * Compute the face normal
   */
  virtual void setFaceNormal(const CFuint faceID, RealVector& normal);
  
private:

  /// socket for Past InwardNormals's
  Framework::DataSocketSink<InwardNormalsData*> socket_pastNormals;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;
  
  /// socket for Past Node's
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;

  /// socket for past State's
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// SubSystem TimeStep
  CFreal                         _dt;

  /// Speed of the face
  RealVector _wallSpeed;
  
}; // end of class UnsteadyWeakSlipWallEuler2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2D_hh
