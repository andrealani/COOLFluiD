#ifndef COOLFluiD_UFEM_BaseBC_hh
#define COOLFluiD_UFEM_BaseBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    template< typename, typename > class DataSocketSink;
    class VectorialFunction;
  }

  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Neumann Boundary condition command
 *
 * @author Tiago Quintino
 * @author reworked for UFEM by Tamas Banyai
 *
 */
class UFEM_API BaseBC : public UFEMSolverCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BaseBC(const std::string& name);

  /// Destructor
  ~BaseBC() {}

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup() {}

  /// Configures the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected:

  /// Execute on the current TRS
  virtual void executeOnTrs();

  /// Calculating the variable values on the state given as argument
  void computeStateValuesBaseBC(const Framework::State* currState);

protected: // data

  /// Set to which equation to apply this BC for all states
  std::vector < std::vector<bool> > m_applyFlags;

  /// Sets the values to apply this BC at current State
  RealVector m_applyVars;

  /// Intermediate states
  Framework::DataSocketSink< Framework::State* > socket_interStates;

  /// States
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL > socket_states;

  /// Vector of strings to hold the functions
  std::vector< std::string > m_functions;

  /// Vector of strings to hold the variables
  std::vector< std::string > m_vars;

  /// VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// Indices of equations to apply the BC to (default: none)
  std::vector< CFuint > m_applyEqs;

  /// Assisting vector for parsing, defined outside computeStateValues function for performance reasons only.
  RealVector m_parseVector;

}; // end of class BaseBC

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_BaseBC_hh

