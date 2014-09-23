#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the N scheme for RDS space discretization
 *
 * @author Andrea Lani
 *
 */
class NSchemeSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  NSchemeSys(const std::string& name);
  
  /**
   * Default destructor
   */
  ~NSchemeSys();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);

  /**
   * Distribute part of the residual
   */
  void distributePart(std::vector<RealVector>& residual);
    
  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

private:

  RealVector _sumKminU;

  RealMatrix _sumKmin;

  RealMatrix _invK;

  RealVector _uInflow;

  RealVector _uDiff;

  RealVector _temp;

  RealMatrix _tempMat;

  RealMatrix _sumKplus;

  RealMatrix _betaLDA;

}; // end of class NSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeSys_hh
