#ifndef COOLFluiD_Numerics_FluctSplit_LDASchemeCSysCUDA_hh
#define COOLFluiD_Numerics_FluctSplit_LDASchemeCSysCUDA_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/RDS_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the N scheme for RDS space discretization
   *
   * @author Andrea Lani
   */
class LDASchemeCSysCUDA : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  LDASchemeCSysCUDA(const std::string& name);

  /**
   * Default destructor
   */
  ~LDASchemeCSysCUDA();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Unset up
   */
  virtual void unsetup();

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
  
  double* dev_a;
  double* dev_b;
  double* dev_c;
  double* dev_d;
  
}; // end of class LDASchemeCSysCUDA

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LDASchemeCSysCUDA_hh
