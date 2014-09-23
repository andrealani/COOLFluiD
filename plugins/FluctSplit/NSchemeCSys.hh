#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeCSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_SplitterSys.hh"

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
class NSchemeCSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  NSchemeCSys(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NSchemeCSys();

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

}; // end of class NSchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeCSys_hh
