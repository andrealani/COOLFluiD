#ifndef COOLFluiD_Numerics_FluctSplit_BSchemeCSysScalar_hh
#define COOLFluiD_Numerics_FluctSplit_BSchemeCSysScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterSysScalar.hh"
#include "FluctSplit/BSchemeBase.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the B scheme (system+scalar) for RDS space discretization
 * based on the CRD approach
 *
 * @author Andrea Lani
 *
 */
class BSchemeCSysScalar : public BSchemeBase<RDS_SplitterSysScalar> {
public:
	
  /**
   * Default constructor.
   */
  BSchemeCSysScalar(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~BSchemeCSysScalar();

  /**
   * Set up
   */
  virtual void setup();
  
  /**
   * Distribute the residual
   */
  virtual void distribute(std::vector<RealVector>& residual);
  
  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

protected: // functions

  /// Compute the blending coefficients
  virtual void computeBlendingCoeff();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation();

private:

  RealVector _sumKplusU;

  RealMatrix _sumKplus;

  RealMatrix _invK;

  RealVector _uInflow;
  
  RealVector _uDiff;
  
  RealVector _phiLDA;
  
  std::vector<RealVector> _phiN;
  
  RealVector m_sumKplusU;
  
  RealVector m_sumKplus;
  
  RealVector m_invCoeff;
  
  RealVector m_uInflow;
 
  RealVector m_uTempSc;
    
}; // end of class BSchemeCSysScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BSchemeCSysScalar_hh
