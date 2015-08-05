#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_FluxSplitter_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_FluxSplitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/FluxSplitter.hh"
#include "Framework/VectorialFunction.hh"

#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a generic FVMCC flux splitter
 *
 * @author Andrea Lani
 *
 */
class FVMCC_FluxSplitter : public Framework::FluxSplitter<CellCenterFVMData> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FVMCC_FluxSplitter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FVMCC_FluxSplitter();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);
  
protected:
  
  /**
   * Compute the mesh speed
   */
  virtual void computeMeshSpeed() 
  {
  }
  
  /**
   * Integrate the flux on the current face
   */
  void integrateFluxOnly(RealVector& result)
  {
    compute(result);
    result *= _faceArea;
  }
  
  /**
   * Integrate the flux and the flux jacobian on the current face
   */
  void integrateFluxAndJacob(RealVector& result)
  {
    integrateFluxOnly(result);
    computeLeftJacobian();
    computeRightJacobian(); 
    _lFluxJacobian *= _faceArea;
    _rFluxJacobian *= _faceArea;
  }
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result) = 0;
  
  /**
   * Compute the left flux jacobian
   */
  virtual void computeLeftJacobian();
  
  /**
   * Compute the right flux jacobian
   */
  virtual void computeRightJacobian();
  
  /**
   * Evaluate the dissipation control function
   */
  virtual void evaluateDissipationControlFunction();
  
  /**
   * Get dissipation control value
   */
  CFreal getDissipationControlCoeff()  
  {
    if (_dissipationControlDef != "") {
      evaluateDissipationControlFunction();
    }
    return _dissipationControlCoeff;
  }

protected:
  
  /// compute flux derivative with respect to pressure
  virtual void dFdP(CFuint side, CFuint iVar, CFreal* row) {}
  
  /// compute flux derivative with respect to density
  virtual void dFdRho(CFuint side, CFuint iVar, CFreal* row) {}
  
  /// compute flux derivative with respect to temperature
  virtual void dFdT(CFuint side, CFuint iVar, CFreal* row) {}
  
  /// compute flux derivative with respect to velocity u
  virtual void dFdU(CFuint side, CFuint iVar, CFreal* row) {}
  
  /// compute flux derivative with respect to velocity v
  virtual void dFdV(CFuint side, CFuint iVar, CFreal* row) {}
  
  /// compute flux derivative with respect to velocity w
  virtual void dFdW(CFuint side, CFuint iVar, CFreal* row) {}
    
  /// compute flux derivative with respect to k
  virtual void dFdK(CFuint side, CFuint iVar, CFreal* row) {}

  /// compute flux derivative with respect to rhou
  virtual void dFdRhoU(CFuint side, CFuint iVar, CFreal* row) {}

      /// compute flux derivative with respect to rhov
  virtual void dFdRhoV(CFuint side, CFuint iVar, CFreal* row) {}

      /// compute flux derivative with respect to rhow
  virtual void dFdRhoW(CFuint side, CFuint iVar, CFreal* row) {}

      /// compute flux derivative with respect to rhoE
  virtual void dFdRhoE(CFuint side, CFuint iVar, CFreal* row) {}
  
private:
  
  /// dissipation control coefficient (must be inaccessible by subclasses)
  CFreal _dissipationControlCoeff;
  
protected:
  
  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// storage of the cell volumes
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// storage of the nodal states
  Framework::DataSocketSink<RealVector> socket_nstates;
  
  /// storage of the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of the nodes
  Framework::DataSocketSink<Framework::Node*,  Framework::GLOBAL> socket_nodes;
  
  /// typedef for pointers to member functions for computing jacobians
  typedef void (FVMCC_FluxSplitter::*FluxDerivative)(CFuint, CFuint, CFreal*);
  
  /// array of objects computing flux derivative  
  std::vector<FluxDerivative> m_fluxDerivative;
  
  /// face area
  CFreal _faceArea;
  
  /// storage of the residual of the first iteration
  CFreal _firstResidual;
  
  /// storage of the previous residual
  CFreal _lastResidual;
  
  /// storage of the maximum residual reached until now
  CFreal _maxResidual;
  
  /// temporary jacobian matrix
  RealMatrix _tmpJacobMatrix;
  
  /// input array (i,r,ri,rl,rmax,cfl) for the dissipation control function
  RealVector _dissipationControlInput;
  
  /// the VectorialFunction to use
  Framework::VectorialFunction _dissipationControlFunction;
  
  /// a single string holding the variables for the function definition
  std::vector<std::string> _dissipationControlVars;
  
  /// a single string holding the function definition
  std::string _dissipationControlDef;
  
}; // end of class FVMCC_FluxSplitter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_FluxSplitter_hh
