#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichletFromFile_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichletFromFile_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
    namespace Common {
    template <class KEY, class VALUE> class LookUpTable;
  }
    
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Dirichlet boundary condition with data from an input file
 *
 * @author Ray Vandenhoeck
 */
class BCDirichletFromFile : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCDirichletFromFile(const std::string& name);

  /// Destructor
  ~BCDirichletFromFile();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCDirichletFromFile";
  }

  /// Setup private data
  void setup();

  /// Unsetup private data
  void unsetup();

  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
   void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                              std::vector< std::vector< RealVector* > >& ghostGrads,
                              const std::vector< RealVector >& normals,
                              const std::vector< RealVector >& coords);

protected: // data

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;
  
  /// dimensional state
  Framework::State* m_dimState;

  /// State holding the input (ghost) variables
  Framework::State* m_input;

  /// a string to hold the name of the input variables
  std::string m_updateVarStr;
  
  protected: // data
  
  /// look up table for u(y)
  std::vector<Common::LookUpTable<CFreal,CFreal>*> m_lookupState;
  
  /// input data file name
  std::string m_infile;
  
  /// a string to hold the name of the input variables
  std::string m_inputVarStr;

  protected: // method
  
  /// fill the lookup table
  void fillTable();

}; // class BCDirichlet

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichletFromFile_hh

