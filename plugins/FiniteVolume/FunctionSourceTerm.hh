#ifndef COOLFluiD_Numerics_FiniteVolume_FunctionSourceTerm_hh
#define COOLFluiD_Numerics_FiniteVolume_FunctionSourceTerm_hh

//////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term given by a user-defined function
 *
 * @author Andrea Lani
 *
 */
class FunctionSourceTerm : public ComputeSourceTermFVMCC {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  FunctionSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FunctionSourceTerm();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();
  
  /**
   * Compute the source term
   * @param element   this GeometricEntity may or may not be used 
   * @param cellID    ID can always be useful in case GeometricEntity is not created
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
protected:
  
  /// RealVector holding the result of the vectorial function
  RealVector _input;
  
  /// RealVector holding the input variables
  RealVector _inputVars;
  
  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _functions;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _vars;
  
}; // end of class FunctionSourceTerm
      
//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FunctionSourceTerm_hh
