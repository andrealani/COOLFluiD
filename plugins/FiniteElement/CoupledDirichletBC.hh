#ifndef COOLFluiD_Numerics_FiniteElement_CoupledDirichletBC_hh
#define COOLFluiD_Numerics_FiniteElement_CoupledDirichletBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Framework/Storage.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Dirichlet Boundary condition command for Coupled SubSystems
 *
 * @author Thomas Wuilbaut
 *
 */
class CoupledDirichletBC : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  CoupledDirichletBC(const std::string& name);

  /**
   * Default destructor
   */
  ~CoupledDirichletBC();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unsetup the private data and data of the aggregated classes
   * in this command after the  processing phase
   */
  void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on the current TRS
   */
  void executeOnTrs();

protected: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for isUpdated
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// socket for boundary states neighours
  Framework::DataSocketSink<std::valarray<Framework::State*> > socket_bStatesNeighbors;

  // the socket to the data handle of arrays of flags specifying if a strong BC has been applied
  // for the given variables in boundary states
  Framework::DataSocketSink<std::vector<bool> > socket_appliedStrongBC;

  /// is this explicit or implicit
  bool _isImplicit;

  /// name of the interface
  std::string _interfaceName;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  ///Flag to say if you want to use (states-pastStates) instead of (states)
  bool _useDeltaStates;

  ///the total number of subiterations
  CFuint _nbSubIterations;

  ///the current subiteration being performed
  CFuint _currentSubIteration;

  ///store the previousIteration
  CFuint _previousIteration;

  ///Do you want to alternate with anotherBC
  bool _alternateBC;

  ///If you alternate, do you want to apply this BC at start
  bool _alternateStart;

  ///If you alternate, should you apply this BC now
  bool _currentAlternateRun;

  /// Symmetry method to apply
  std::string m_symmetryStr;

  /// ScaleDiagonal symmetry method coefficient
  CFreal m_scale;

  /// Indices of equations to apply the BC to (default: all of PhysicalModel)
  std::vector< CFuint > m_applyEqs;

  /// Backup of the values imposed at the previous iteration
  std::vector< RealVector > m_pastStates;

}; // end of class CoupledDirichletBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_CoupledDirichletBC_hh
