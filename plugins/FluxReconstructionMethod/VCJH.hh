#ifndef COOLFluiD_FluxReconstructionMethod_VCJH_hh
#define COOLFluiD_FluxReconstructionMethod_VCJH_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {
    
    class FluxReconstructionElementData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Vincent-Castonguay-Jameson-Huyn Correction Function (also called Energy Stable FR schemes)
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 * @author Rayan Dhib
 */
class VCJH : public BaseCorrectionFunction {

public:  // methods

  /// Defines the Config Options of this class
  static void defineConfigOptions(Config::OptionList& options);
    
  /// Configures this Method.
  virtual void configure ( Config::ConfigArgs& args );
    
  /// Constructor
  VCJH(const std::string& name);

  /// Destructor
  ~VCJH();
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "VCJH";
  }
  
  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();
    
  /**
   * Compute the VCJH correction function of an instance of FluxReconstructionElementData. 
   * corrfcts contains for each solution point the contrbution of each flux point (being a realvector with the dimensionality of the test case)
   */
  void computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > >& corrcts);
  
  /**
   * Compute the a specified polynomial order VCJH correction function of an instance of FluxReconstructionElementData, used for LLAV. 
   * corrfcts contains for each solution point the contrbution of each flux point (being a realvector with the dimensionality of the test case)
   */
  void computeCorrectionFunction(const CFPolyOrder::Type solOrder, const CFreal factor, Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrcts);
    
  /**
   * Compute the divergence of the VCJH correction function of an instance of FluxReconstructionElementData. 
   * corrfcts contains for each solution point the contrbution of each flux point
   */
  void computeDivCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrcts);
    
private : // helper functions
  /// Compute the value of the VCJH 1D correction function of order p and cfactor at the 1D coordinate ksi
  CFreal computeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor);

  /// Compute the value of the derivative of the VCJH 1D correction function of order p and cfactor at the 1D coordinate ksi
  CFreal computeDerivativeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor);

  /// Compute the rhs of the system of equations to be solved for Triangles correction function (for triag)
  RealVector computeIntRHSTriag(const CFPolyOrder::Type solOrder, CFuint f, CFuint j);

  /// Compute the rhs of the system of equations to be solved for Triangles correction function (for Tetra)
  RealVector computeIntRHSTetra(const CFPolyOrder::Type solOrder, CFuint f, CFuint j);

  /// Compute the orthonormal dubiner basis in 2D (for triag)
  RealVector computeDubiner2D(const CFPolyOrder::Type solOrder,CFreal ksi, CFreal eta);

  /// Compute the orthonormal dubiner basis in 3D (for Tetra)
  RealVector computeDubiner3D(const CFPolyOrder::Type solOrder,CFreal ksi, CFreal eta, CFreal zta);

  /// Compute the orthonormal normalized jacobi poly (for simplex)
  CFreal ComputeJacobi(CFuint N,CFuint alpha, CFuint beta,CFreal x);

  /// Compute the sigmas for building the div of the correction function on triangles (for triag)
  RealVector computeSigmasTriag(const CFPolyOrder::Type solOrder,CFuint f, CFuint j, CFreal ksi, CFreal eta, CFreal cfactor);

  /// Compute the sigmas for building the div of the correction function on triangles (for Tetra)
  RealVector computeSigmasTetra(const CFPolyOrder::Type solOrder,CFuint f, CFuint j, CFreal ksi, CFreal eta, CFreal zta, CFreal cfactor);

  /// Compute the divergence of the correction functions for a solution point (for triag)
  RealMatrix computeTriagDivCorrFct(const CFPolyOrder::Type solOrder, CFreal cfactor, CFreal ksi, CFreal eta);

  /// Compute the divergence of the correction functions for a solution point (for Tetra)
  RealMatrix computeTetraDivCorrFct(const CFPolyOrder::Type solOrder, CFreal cfactor, CFreal ksi, CFreal eta, CFreal zta);

  /// Compute the left hand side matrix of the system of equations to be solved (for triag)
  RealMatrix computeLhsTriag(const CFPolyOrder::Type solOrder, CFreal ksi, CFreal eta);

  /// Compute the left hand side matrix of the system of equations to be solved (for Tetra)
  RealMatrix computeLhsTetra(const CFPolyOrder::Type solOrder);

  /// Computes the factorial
  CFreal factorial(CFreal n);

  /// Round a floating-point value to the nearest integer
  CFuint round(double x);

  /// Compute inverse matrix
  void InvertMatrix(RealMatrix A, RealMatrix& AI);

  /// Swap rows
  void SwapRows(RealMatrix& A, CFuint row1, CFuint row2);

private : // private data
  /// Value of the C factor for VCJH 1D correction function
  CFreal  m_cfactor;

}; // class VCJH

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_VCJH_hh
