#ifndef COOLFluiD_Numerics_FluctSplit_RDHLLSchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_RDHLLSchemeCSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the HLL scheme for RDS space discretization
 * based on the CRD approach
 *
 * @author Andrea Lani
 * @author Jesus Garicano Mena
 *
 */
class RDHLLSchemeCSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  RDHLLSchemeCSys(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RDHLLSchemeCSys();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Distribute the residual
   */
  virtual void distribute(std::vector<RealVector>& residual);
  /**
   * Distribute the residual
   */
  virtual void distributePart(std::vector<RealVector>& residual);
  
  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);
  
protected:

  RealVector  _sumKplusU;

  RealMatrix  _sumKplus;

  RealMatrix  _invK;

  RealVector  _uInflow;

  RealVector  _uDiff;

  RealVector  _temp;
  
  RealVector  _tempBkp;
  
  RealMatrix _tempMat;

  RealMatrix _tmp;
  
  RealVector _sumKU;
  
  RealVector _pData;
  
  RealVector n_IJ;
  RealVector n_IK;
  
  /// Euler flux function along x at node I
  RealVector F_I;
  /// Euler flux function along y at node I
  RealVector G_I;
  /// Id
  RealVector F_J;
  /// Id
  RealVector G_J;
  /// Id
  RealVector F_K;
  /// Id
  RealVector G_K;

  /// Numerical FV flux across section IJ of the convex hull of node I
  RealVector Fn_IJ;
  ///  Idem for section IK
  RealVector Fn_IK;

}; // end of class RDHLLSchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RDHLLSchemeCSys_hh
