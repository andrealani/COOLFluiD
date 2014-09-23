#ifndef COOLFluiD_Numerics_FiniteVolume_WMUSCLPolyRec_hh
#define COOLFluiD_Numerics_FiniteVolume_WMUSCLPolyRec_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_PolyRec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
/**
 * This class implements an upwind biased polynomial reconstructor for FVM, 
 * only suitable for structured  fully quadrilateral or hexahedral meshes
 *
 * @author Andrea Lani
 */
class WMUSCLPolyRec : public FVMCC_PolyRec {
public:

  /**
   * Constructor
   */
  WMUSCLPolyRec(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~WMUSCLPolyRec();

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
  virtual void computeGradients();

  /**
   * Set up the private data
   */
  virtual void setup();

  /**
   * Update the weights when nodes are moving
   */
  virtual void updateWeights();

protected:

  /**
   * Extrapolate the solution in the face quadrature points
   */
  virtual void extrapolateImpl(Framework::GeometricEntity* const face);

  /**
   * Extrapolate the solution in the face quadrature points
   */
  virtual void extrapolateImpl(Framework::GeometricEntity* const face,
			       CFuint iVar, CFuint leftOrRight);
  
  /**
   * Extrapolate the solution on a inner face
   */
  virtual void extrapolateInner(Framework::GeometricEntity* const face);
  
  /**
   * Extrapolate the solution on a boundary face
   */
  virtual void extrapolateBoundary(Framework::GeometricEntity* const face);
  
  /**
   * Compute the mid node on the given face
   */
  void computeMidFaceNode(Framework::GeometricEntity* const face);
  
  /**
   * Compute the linear reconstruction for the left state
   */
  void computeReconstruction(Framework::GeometricEntity* const face, 
			     CFuint leftOrRight, CFuint iVar);
    
protected:
  
  /// socket for stencil
  Framework::DataSocketSink<
                            std::vector<Framework::State*> > socket_stencil;
  
  // update states array
  std::vector<Framework::State*> _upStates; 
  
  // back up left state
  RealVector _lStateBkp;
  
  // back up right state
  RealVector _rStateBkp;
  
  // limiter value
  RealVector _vLimiter;

  // temporary solution variations (U_i - U_i-1)
  RealVector _drUlUll;
  
  // temporary solution variations (U_i+1 - U_i)
  RealVector _drUrUl; 
  
  // temporary solution variations (U_i+2 - U_i+1)
  RealVector _drUrrUr;

  // temporary limiter argument
  RealVector _r1;
  
  // temporary limiter argument
  RealVector _r2;
  
  // temporary face normal
  RealVector _faceNormal;
    
  //temporary face mid node (lying on the normal to the face and passing via the inner node)
  RealVector _midNode;
  
  // distance between the inner cell center and the boundary face
  CFreal _drXiXw;
  
  // distance between x_0 and x_1, x_0 and x_mid_face
  RealMatrix _dx;
  
  // reconstruction stencil pattern
  MathTools::CFMat<CFuint> _stencilPattern;
  
}; // end of class LeastSquareP1PolyRec2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_WMUSCLPolyRec_hh
