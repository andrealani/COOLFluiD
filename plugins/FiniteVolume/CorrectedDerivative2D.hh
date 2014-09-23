#ifndef COOLFluiD_Numerics_FiniteVolume_CorrectedDerivative2D_hh
#define COOLFluiD_Numerics_FiniteVolume_CorrectedDerivative2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "DerivativeComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes derivatives using a diamond volume in 2D
 *
 * @author Andrea Lani
 *
 */
class CorrectedDerivative2D : public DerivativeComputer {
public:

  /**
   * Constructor
   */
  CorrectedDerivative2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CorrectedDerivative2D();

  /**
   * Set up the member data
   */
  virtual void setup();

  /*
   * Compute the gradients
   */
  virtual void computeGradients(const RealMatrix& values,
			std::vector<RealVector*>& gradients);
  
  /**
   * Compute the control volume around the current face
   */
  virtual void computeControlVolume(std::vector<RealVector*>& states, 
				    Framework::GeometricEntity *const geo);
  
  /**
   * Compute the average values corresponding to the given values
   */
  virtual void computeAverageValues(Framework::GeometricEntity *const geo,
				    const std::vector<RealVector*>& values,
				    RealVector& avValues);
  
 
  /**
   * Get the maximum number of vertices in the control volume
   */
  virtual CFuint getMaxNbVerticesInControlVolume() const
  {
    return 100;
  }

  /**
   * Get the current number of vertices in the control volume
   */
  virtual CFuint getNbVerticesInControlVolume(Framework::GeometricEntity *const geo) const
  {
    return _nbDofsInCV;
  }

  /**
   * Get the jacobian of the gradients
   */
  virtual Common::SafePtr<std::vector<RealVector> > getGradientsJacob();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:
  
  /**
   * Set the control volume states starting form the given data
   */
  virtual void setControlVolumeStates(std::vector<RealVector*>& states);
  
  /*
   * Compute the inner gradients
   */
  virtual void computeInnerGradients(const RealMatrix& values,
				     std::vector<RealVector*>& gradients);

  /*
   * Compute the boundary gradients
   */
  virtual void computeBoundaryGradients(const RealMatrix& values,
					std::vector<RealVector*>& gradients);

  /*
   * Correct the given gradient
   */
  void correctGradient(const CFreal& grad, const RealVector& gradLR, RealVector& result) const
  {
    result = grad*_eRLdotNN + _oEminusNN*gradLR;
  }
 
  /*
   * Compute 1 - n x n
   */
  virtual void computeOEminusNN(const RealVector& n);
  
protected: // data

  /// degrees of freedom in the control volume
  CFuint _nbDofsInCV;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of states
  Framework::DataSocketSink< Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of nodes
  Framework::DataSocketSink< Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// socket for stencil
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// storage for uX
  Framework::DataSocketSink<CFreal> socket_uX;

  /// storage for uY
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// pointer to the cell TRS
  Common::SafePtr<Framework::TopologicalRegionSet> _cellsTRS;

  /// distance between the states RL
  CFreal _dr;

  /// direction vector for RL
  RealVector _eRL;

  /// (e dot n) multiplied normal
  RealVector _eRLdotNN;

  /// (1 - n x n)
  RealMatrix _oEminusNN;

  /// temporary gradient
  RealVector _gradLR;

  /// temporary array for RHS in LS algorithm
  RealVector _lf1;

  /// temporary array for RHS in LS algorithm
  RealVector _lf2;

  /// cell builder
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _cellBuilder;
  
  /// array storing the neighbor states
  std::vector<Framework::State*> _stencil;

}; // end of class CorrectedDerivative2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CorrectedDerivative2D_hh
