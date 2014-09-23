#ifndef COOLFluiD_Numerics_FluctSplit_LDASchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_LDASchemeSys_hh

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
 *
 *
 */
class LDASchemeSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  LDASchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~LDASchemeSys();

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

  RealVector _phi;

  RealVector _uTemp;

  RealVector _temp;

  RealMatrix _sumKplus;

  RealMatrix _k;

  RealMatrix _invK;

  RealMatrix _beta;

}; // end of class LDASchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeSys_hh
