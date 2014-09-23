#ifndef COOLFluiD_Numerics_FiniteVolume_CorrectedDerivativeGG2D_hh
#define COOLFluiD_Numerics_FiniteVolume_CorrectedDerivativeGG2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeometricEntityPool.hh"
#include "Framework/CellTrsGeoBuilder.hh"

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
class CorrectedDerivativeGG2D : public DerivativeComputer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CorrectedDerivativeGG2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CorrectedDerivativeGG2D();

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
  Common::SafePtr<std::vector<RealVector> > getGradientsJacob();

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

  /*
   * Compute the gradients with Green Gauss
   */
  void computeGradientsGG(Common::SafePtr<Common::Table<CFuint> > faceNodes,
			  std::vector<CFuint> cellFaceIDs,
			  Framework::State *const state,
			  const RealMatrix& values,
			  CFuint counter,
			  RealMatrix& ux);

protected: // data

  /// degrees of freedom in the control volume
  CFuint _nbDofsInCV;

  /// storage of face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of states
  Framework::DataSocketSink< Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of nodes
  Framework::DataSocketSink< Framework::Node* , Framework::GLOBAL > socket_nodes;

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

  /// gradients XY of all left variables
  RealMatrix _uxL;

  /// gradients XY of all right variables
  RealMatrix _uxR;

  /// temporary face normal
  RealVector _fnormal;

  /// cell builder
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> _cellBuilder;

  /// array storing the neighbor states
  std::vector<Framework::State*> _stencil;

  /// matrix storing the IDs of the nodal states in the left cell, face by face
  Common::SafePtr<Common::Table<CFuint> >  _faceNodesL;

  /// matrix storing the IDs of the nodal states in the right cell, face by face
  Common::SafePtr<Common::Table<CFuint> > _faceNodesR;

  /// current face nodes IDs in the extrapolated values storage
  std::vector<CFuint> _currFaceNodeIDs;

  /// face IDs beloning to the left cell
  std::vector<CFuint> _cellFaceIDsL;

  /// face IDs beloning to the right cell
  std::vector<CFuint> _cellFaceIDsR;

  /// number of nodes in left cell
  CFuint _nbNodesInCellL;

  /// number of nodes in right cell
  CFuint _nbNodesInCellR;

  /// flag telling if to use weights in the LS algo
  bool _useWeights;

}; // end of class CorrectedDerivativeGG2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CorrectedDerivativeGG2D_hh
