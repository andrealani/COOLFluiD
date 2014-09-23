#ifndef COOLFluiD_Numerics_FluctSplit_PSISchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_PSISchemeSys_hh

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
class PSISchemeSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  PSISchemeSys(const std::string& name);

  /**
   * Default destructor
   */
  ~PSISchemeSys();

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

private:

  RealMatrix _k;

  RealMatrix _sumKmin;

  RealMatrix _invK;

  RealVector _sumKminU;

  RealVector _uTemp;

  RealVector _uDiff;

  RealVector _phi;

  RealVector _sumBeta;

  RealVector _invCoeff;

  RealVector _temp;

}; // end of class PSISchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PSISchemeSys_hh
