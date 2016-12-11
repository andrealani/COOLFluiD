#ifndef COOLFluiD_FluxReconstructionMethod_VCJH_hh
#define COOLFluiD_FluxReconstructionMethod_VCJH_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {
    
    class FluxReconstructionElementData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Vincent-Castonguay-Jameson-Huyn Correction Function
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 */
class VCJH : public BaseCorrectionFunction {

public:  // methods

  /// Constructor
  VCJH(const std::string& name);

  /// Destructor
  ~VCJH();
    
  /// Compute the VCJH correction function of an instance of FluxReconstructionElementData
  void computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, const std::vector< std::vector< RealVector > > corrcts);
    
    
private :
  /// Compute the value of the VCJH 1D correction function of order p at ksi
  CFreal computeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "VCJH";
  }

  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();


}; // class VCJH

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_VCJH_hh
