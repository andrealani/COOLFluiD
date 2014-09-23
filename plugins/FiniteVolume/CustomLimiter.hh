#ifndef COOLFluiD_Numerics_FiniteVolume_CustomLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_CustomLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Limiter.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "MathTools/FunctionParser.hh"
#include "CellCenterFVMData.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a family of user defined limiters associated to 
 * high-order reconstruction
 *
 * @author Andrea Lani
 *
 */
class CustomLimiter : public Framework::Limiter<CellCenterFVMData> {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CustomLimiter(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~CustomLimiter();
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Compute the limiter in the current face
   */
  virtual void limit(const std::vector<std::vector<Framework::Node*> >& coord,
		     Framework::GeometricEntity* const cell,
		     CFreal* limiterValue)
  {
    throw Common::NotImplementedException (FromHere(),"CustomLimiter::limit()");
  }
  
  /**
   * Apply the face limiter: entries for limiter values are ordered as follows:
   * [L0,R0,L1,R1, ...] with L=left and R=right values
   */
  virtual void limitOnFace(const RealVector& rLeft, 
			   const RealVector& rRight,
			   CFreal* limiterValue);
  
  /**
   * Apply the limiter to a scalar quantity
   */
  virtual void limitScalar(CFreal r, CFreal& limiterValue);
  
  /**
   * Set up the private data
   */
  virtual void setup();
  
protected:
  
  /// denominator of limiter argument 
  CFreal _deltaMin;
  
  /// one-component array to keep the current variable gradient ratio
  RealVector _gradRatio;

  /// limiter function parser
  MathTools::FunctionParser _functionParser;
  
  /// map that already contains some default limiter functions
  Common::CFMap<std::string,std::string> _mapNameToDefaultFunction;
  
  /// a string holding the function name or expression
  std::string _function;

  /// limiter name corrresponding to limiter functions in archive
  std::string _limiterName;
  
}; // end of class CustomLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CustomLimiter_hh
