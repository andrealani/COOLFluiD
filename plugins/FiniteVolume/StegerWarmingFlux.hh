#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the StegerWarming flux
 *
 * @author Andrea Lani
 *
 */
class StegerWarmingFlux : public FVMCC_FluxSplitter {
public:

  /**
   * Constructor
   */
  StegerWarmingFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~StegerWarmingFlux();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
  /**
   * Compute the left flux jacobian
   */
  virtual void computeLeftJacobian();
  
  /**
   * Compute the right flux jacobian
   */
  virtual void computeRightJacobian();
  
protected:
  
  /**
   * Compute numerically the variable transformation matrix
   */
  void computeTransformMatrix(Framework::State* currState);
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _solutionVarSet;
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _updateVarSet;
  
  /// variable set transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> _updateToSolutionVarTrans;
  
  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
  /// Numerical jacobian calculator
  std::auto_ptr<Framework::NumericalJacobian> _numJacob;
  
  /// matrix of right eigenvectors
  RealMatrix   _jacobPlus;

  /// matrix of left eigenvectors
  RealMatrix   _jacobMin;
  
  /// matrix of right eigenvectors
  RealMatrix   _jacobPlusTrasf;
  
  /// matrix of left eigenvectors
  RealMatrix   _jacobMinTrasf;
  
  /// dummy matrix of eigenvectors
  RealMatrix   _jacobDummy;
  
  /// vector of eigenvalues
  RealVector   _eValues;

  /// vector of right state eigenvalues
  RealVector   _rightEvalues;

  /// vector of left state eigenvalues
  RealVector   _leftEvalues;
  
  /// temporary state
  RealVector   _tState;
  
  /// temporary flux
  RealVector   _tFlux;
  
  /// temporary unit normal
  RealVector   _tempUnitNormal;
  
}; // end of class StegerWarmingFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingFlux_hh
