#ifndef COOLFluiD_Numerics_FluctSplit_PSISchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_PSISchemeCSys_hh

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
   */
class PSISchemeCSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  PSISchemeCSys(const std::string& name);

  /**
   * Default destructor
   */
  ~PSISchemeCSys();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Distribute the residual
   */
  virtual void distribute(std::vector<RealVector>& residual);

  /**
   * Distribute part of the residual
   */
  virtual void distributePart(std::vector<RealVector>& residual);

private:

  RealMatrix _sumKplus;

  RealMatrix _invK;

  RealVector _sumKplusU;

  RealVector _uInflow;

  RealVector _uDiff;

  RealVector _sumBeta;

  RealVector _invCoeff;

  RealVector _temp;

}; // end of class PSISchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PSISchemeCSys_hh
