#ifndef COOLFluiD_Numerics_FiniteVolume_RoeFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeFlux_hh

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
class RoeFlux : public FVMCC_FluxSplitter {
public:

  /**
   * Constructor
   */
  RoeFlux(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~RoeFlux();
  
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
  RealVector   _sumFlux;

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

  /// vector of eigenvalues
  RealVector   _absEvalues;

  /// temporary state
  RealVector   _tState;

  /// abs of the jacobian matrix
  RealMatrix   _absJacob;

  /// right jacobian matrix
  RealMatrix   _jRight;
  
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

}; // end of class RoeFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeFlux_hh
