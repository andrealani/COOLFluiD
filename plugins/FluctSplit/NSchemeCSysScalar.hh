#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeCSysScalar_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeCSysScalar_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_SplitterSysScalar.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the N scheme for RDS space discretization
 * based on the CRD approach
 *
 * @author Andrea Lani
 *
 */
class NSchemeCSysScalar : public RDS_SplitterSysScalar {
public:
  
  /**
   * Default constructor.
   */
  NSchemeCSysScalar(const std::string& name);

  /**
   * Default destructor
   */
  ~NSchemeCSysScalar();

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
  
private:

  RealVector  _sumKplusU;

  RealMatrix  _sumKplus;

  RealMatrix  _invK;
  
  RealMatrix _tempMat;
  
  RealVector m_sumKplusU;
  
  RealVector m_sumKplus;
  
}; // end of class NSchemeCSysScalar

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeCSysScalar_hh
