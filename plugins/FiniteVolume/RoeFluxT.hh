#ifndef COOLFluiD_Numerics_FiniteVolume_RoeFluxT_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeFluxT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux corresponding to the Euler
 * physical model 2D (in conservative variables)
 *
 * @author Andrea Lani
 *
 */
template <int N>
class RoeFluxT : public FVMCC_FluxSplitter {
public:

  /**
   * Constructor
   */
  RoeFluxT(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~RoeFluxT();
  
  /** 
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure
   */
  virtual void configure ( Config::ConfigArgs& args );
  
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
  
protected: // helper functions
  
  /**
   * Compute numerically the variable transformation matrix
   */
  void computeTransformMatrix(Framework::State* currState);
    
  /**
   * Set the abs of the eigen values
   */
  virtual void setAbsEigenValues();
  
  /**
   * Perform ad-hoc linearization
   */
  virtual void linearize();
  
  /**
   * Compute the artificial diffusion reduction coefficient
   */
  virtual CFreal getReductionCoeff() {
    return _currentDiffRedCoeff;
  }
  
private:
  
  /// Coefficient to reduce the diffusive part
  CFreal _currentDiffRedCoeff;
  
protected:
  
  /// array storing the sum of the right and left flux
  MathTools::CFVec<CFreal,N> _sumFlux;
  
  /// vector of eigenvalues
  MathTools::CFVec<CFreal,N> _absEvalues;
  
  /// abs of the jacobian matrix
  MathTools::CFMat<CFreal,N,N> _absJacob;
  
  /// matrix of right eigenvectors
  RealMatrix   _rightEv;

  /// matrix of left eigenvectors
  RealMatrix   _leftEv;
  
  /// vector of eigenvalues
  RealVector   _eValues;

  /// vector of right state eigenvalues
  RealVector   _rightEvalues;

  /// vector of left state eigenvalues
  RealVector   _leftEvalues;
  
  /// temporary state
  RealVector   _tState;
  
  /// right jacobian matrix
  RealMatrix _jRight;
  
  /// left jacobian matrix
  RealMatrix   _jLeft;
  
  /// jacobian matrix
  RealMatrix   _jacob;
  
  /// jacobian matrix
  RealMatrix   _jacobLeftTransf;
  
  /// jacobian matrix
  RealMatrix   _jacobRightTransf;
  
  /// dummy jacobian matrix
  RealMatrix   _jacobDummy;

  /// temporary unit normal
  RealVector   _tempUnitNormal;

  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;

}; // end of class RoeFluxT

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "RoeFluxT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeFluxT_hh
