#ifndef COOLFluiD_Numerics_FluctSplit_RusanovSchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_RusanovSchemeCSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Rusanov scheme for RDS space discretization
 * based on the CRD approach
 *
 * @author Tiago Quintino
 */
class RusanovSchemeCSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  RusanovSchemeCSys(const std::string& name);
  
  /**
   * Default destructor
   */
  ~RusanovSchemeCSys();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Distribute the residual
   */
  virtual void distribute(std::vector<RealVector>& residual);
  
private: // data

   // diffusion coefficient
   CFreal m_alpha;

   // maximum eigenvalue
   CFreal m_maxEV;

   // minimum eigenvalue
   CFreal m_minEV;

   RealVector m_sumUmin;

}; // end of class RusanovSchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RusanovSchemeCSys_hh
