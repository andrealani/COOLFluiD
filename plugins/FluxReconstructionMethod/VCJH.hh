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

  /// Defines the Config Options of this class
  static void defineConfigOptions(Config::OptionList& options);
    
  /// Configures this Method.
  virtual void configure ( Config::ConfigArgs& args );
    
  /// Constructor
  VCJH(const std::string& name);

  /// Destructor
  ~VCJH();
    
  /// Compute the VCJH correction function of an instance of FluxReconstructionElementData
  void computeCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< RealVector > >& corrcts);
    
  /// Compute the divergence of the VCJH correction function of an instance of FluxReconstructionElementData
  void computeDivCorrectionFunction(Common::SafePtr< FluxReconstructionElementData > frElemData, std::vector< std::vector< CFreal > >& corrcts);
    
private : // helper functions
  /// Compute the value of the VCJH 1D correction function of order p and cfactor at the 1D coordinate ksi
  CFreal computeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor);

  /// Compute the value of the derivative of the VCJH 1D correction function of order p and cfactor at the 1D coordinate ksi
  CFreal computeDerivativeCorrectionFunction1D(CFPolyOrder::Type solOrder, CFreal ksi, CFreal cfactor);

private : // private data
  /// Value of the C factor for VCJH 1D correction function
  CFreal  m_cfactor;
    
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
