#ifndef COOLFluiD_Numerics_FiniteElement_DirichletBC_hh
#define COOLFluiD_Numerics_FiniteElement_DirichletBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    template< typename, typename > class DataSocketSink;
    class VectorialFunction;
  }

  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Dirichlet Boundary condition command
 *
 * @author Tiago Quintino
 *
 */
class DirichletBC : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  DirichletBC(const std::string& name);

  /// Destructor
  ~DirichletBC() {}

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unsetup the private data and data of the aggregated classes
   * in this command after the processing phase
   */
  void unsetup() {}

  /// Configures the command
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:

  /// Execute on the current TRS
  void executeOnTrs();

protected: // data

  /// Socket for RHS
  Framework::DataSocketSink<CFreal > socket_rhs;

  /// Socket for isUpdated
  Framework::DataSocketSink<bool > socket_isUpdated;

  /// Socket for states
  Framework::DataSocketSink<Framework::State* , Framework::GLOBAL> socket_states;
  
  /// Socket for neighbor states
  Framework::DataSocketSink<
    std::valarray< Framework::State* > > socket_bStatesNeighbors;

  /// the socket to the data handle of arrays of flags specifying if a strong BC has been applied
  /// for the given variables in boundary states
  Framework::DataSocketSink<std::vector<bool> > socket_appliedStrongBC;

  /// Vector of strings to hold the functions
  std::vector< std::string > m_functions;

  /// Vector of strings to hold the variables
  std::vector< std::string > m_vars;

  /// VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// Method is Implicit or explicit?
  bool m_isImplicit;

  /// Symmetry method to apply
  std::string m_symmetryStr;

  /// ScaleDiagonal symmetry method coefficient
  CFreal m_scale;

  /// Indices of equations to apply the BC to (default: all of PhysicalModel)
  std::vector< CFuint > m_applyEqs;

}; // end of class DirichletBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_DirichletBC_hh

