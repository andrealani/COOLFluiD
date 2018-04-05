#ifndef COOLFluiD_Numerics_FluctSplit_FunctionSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_FunctionSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
  
  namespace Framework {
    class State; 
  }
  
    namespace FluctSplit {
      
      class InwardNormalsData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term given by a user-defined function
 * 
 * @author Andrea Lani
 */
class FunctionSourceTerm : public ComputeSourceTermFSM {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  FunctionSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~FunctionSourceTerm();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /// Compute the source term
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				RealVector& source,
				const FluctSplit::InwardNormalsData& normalsData);  
  
  /// Compute the source term
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				std::vector<RealVector>& source,
				const FluctSplit::InwardNormalsData& normalsData);
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
private: // data
    
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

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_FunctionSourceTerm_hh
