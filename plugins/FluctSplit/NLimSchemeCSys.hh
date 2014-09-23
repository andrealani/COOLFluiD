#ifndef COOLFluiD_Numerics_FluctSplit_NLimSchemeCSys_hh
#define COOLFluiD_Numerics_FluctSplit_NLimSchemeCSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "RDS_SplitterSys.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the N limited scheme for RDS space discretization
 * based on the CRD approach
 *
 * @author Nadege Villedieu
 *
 */
class NLimSchemeCSys : public RDS_SplitterSys {
public:

  /**
   * Default constructor.
   */
  /**
   * 
   * @param name 
   */
  NLimSchemeCSys(const std::string& name);

  /**
   * Default destructor
   */
  ~NLimSchemeCSys();
 /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
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

  CFreal m_angle;

  RealVector _normal;

  RealVector _phy;

  RealVector _phyChar;

  std::vector<RealVector> _residualChar;

  RealMatrix _rightEigenVector;

  RealMatrix _leftEigenVector;

  RealVector _sumBeta;

  std::vector<RealVector> _beta;

  std::vector<RealVector> _betaLim; 

   /// temporary data for holding eignvalues
  RealVector                       _eValues;


}; // end of class NLimSchemeCSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NLimSchemeCSys_hh
