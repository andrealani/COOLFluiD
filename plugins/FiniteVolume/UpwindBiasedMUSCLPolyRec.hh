#ifndef COOLFluiD_Numerics_FiniteVolume_UpwindBiasedMUSCLPolyRec_hh
#define COOLFluiD_Numerics_FiniteVolume_UpwindBiasedMUSCLPolyRec_hh

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
class UpwindBiasedMUSCLPolyRec : public FVMCC_PolyRec {
public:

  /**
   * Constructor
   */
  UpwindBiasedMUSCLPolyRec(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~UpwindBiasedMUSCLPolyRec();

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
   * Compute the limiters for all the cells
   */
  virtual void computeLimiters();

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
  
protected:
  
  /// constant for 1./6.
  const CFreal _sixth;
  
  /// constant for 1./3.
  const CFreal _third;
  
  /// constant 
  const CFreal _k1;
    
  /// socket for stencil
  Framework::DataSocketSink<
                            std::vector<Framework::State*> > socket_stencil;
  
  // update states array
  std::vector<Framework::State*> _upStates; 
  
  // back up left state
  RealVector _lStateBkp;
  
  // back up right state
  RealVector _rStateBkp;
    
}; // end of class LeastSquareP1PolyRec2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UpwindBiasedMUSCLPolyRec_hh
