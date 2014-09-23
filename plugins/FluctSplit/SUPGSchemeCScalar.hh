#ifndef COOLFluiD_Numerics_FluctSplit_SUPGSchemeCScalar_hh
#define COOLFluiD_Numerics_FluctSplit_SUPGSchemeCScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterScalar.hh"
#include "FluctSplit/FluctSplitScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the SUPG scheme for RDS space discretization
/// based on the CRD approach
/// The implementation was taken as-is from THOR code.
/// @author Tiago Quintino
class FluctSplitScalar_API SUPGSchemeCScalar : public RDS_SplitterScalar {
public:

  /// Default constructor.
  SUPGSchemeCScalar(const std::string& name);

  /// Default destructor
  ~SUPGSchemeCScalar();

  /// Set up
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  /// Distribute the residual
  virtual void distributePart(std::vector<RealVector>& residual);

  /// Compute all the contributions for the Picard jacobian
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private: // data

  /// C1 constant of the scheme
  CFreal m_C1;

  /// C2 constant of the scheme
  CFreal m_C2;

  /// average cell size
  CFreal m_h;

  /// dimension
  CFuint  m_dim;

  /// inverse of the dimension
  CFreal m_invdim;

  /// inverse of the cell area
  CFreal m_invCellArea;

  /// tau parameter
  RealVector m_tau;

  /// lambda parameter
  RealVector m_lambda;

  /// average speed norm parameter
  std::vector<RealVector> m_avgSpeed;

  /// norm of the gradient of the variables
  RealVector m_normGrad;

  /// kappa hat parameter
  RealVector m_kappa_hat;

  /// gradient of each variable
  std::vector<RealVector> m_grad;

  /// normal of eact state
  std::vector<RealVector> m_stateNormal;

  /// p parameter
  std::vector<RealVector> m_p;

}; // end of class SUPGSchemeCScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SUPGSchemeCScalar_hh
