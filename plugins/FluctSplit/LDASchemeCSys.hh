#ifndef COOLFluiD_Numerics_FluctSplit_LDASchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_LDASchemeCSys_hh

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
class LDASchemeCSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  LDASchemeCSys(const std::string& name);

  /**
   * Default destructor
   */
  ~LDASchemeCSys();

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
private:

  RealMatrix _sumKplus;

  RealMatrix _invK;

  RealVector _uTemp;

  RealMatrix _k;

  RealMatrix _beta;

}; // end of class LDASchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeCSys_hh
