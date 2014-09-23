#ifndef COOLFluiD_Numerics_FluctSplit_NSchemeSysT_hh
#define COOLFluiD_Numerics_FluctSplit_NSchemeSysT_hh

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
template <CFuint N>
class NSchemeSysT : public RDS_SplitterSys {
public:
  typedef NSchemeSysT<N> SELF;
  
  /**
   * Default constructor.
   */
  NSchemeSysT(const std::string& name);

  /**
   * Default destructor
   */
  ~NSchemeSysT();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Distribute the residual
   */
  void distribute(std::vector<RealVector>& residual);

private:
  
  MathTools::CFVec<CFreal,N> _sumKminU;
  
  MathTools::CFVec<CFreal,N> _uInflow;

  MathTools::CFVec<CFreal,N> _uDiff;

  MathTools::CFVec<CFreal,N> _temp;

  MathTools::CFMat<CFreal,N,N> _tempMat;

  MathTools::CFMat<CFreal,N,N> _sumKplus;
  
  MathTools::CFMat<CFreal,N,N> _betaLDA;
  
  RealMatrix  _sumKmin;
  RealMatrix  _invK;
  
}; // end of class NSchemeSysT

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NSchemeSysT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeSysT_hh
