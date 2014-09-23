#ifndef COOLFluiD_SpectralFD_LUSGSSetup_hh
#define COOLFluiD_SpectralFD_LUSGSSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/StdSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This is a standard command to setup the SpectralFD method to be used
 * in combination with the LUSGSMethod convergence method
 *
 * @author Kris Van den Abeele
 */
class LUSGSSetup : public StdSetup {

public: // public functions

  /// Constructor
  explicit LUSGSSetup(const std::string& name);

  /// Destructor
  ~LUSGSSetup();

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

private: // private functions

public: // public data

  /// Node to state ID map
  std::vector< CFuint > m_nodeIDToStateID;

protected: // protected data

  /// socket for diagonal block Jacobian matrices
  Framework::DataSocketSource< RealMatrix > socket_diagBlockJacobMatr;

  /// socket for rhs of current set of states
  Framework::DataSocketSource< CFreal > socket_rhsCurrStatesSet;

  /// handle to the IDs of the states in each set of states
  Framework::DataSocketSource< std::vector< CFuint > > socket_statesSetStateIDs;

  /// handle to list of booleans telling whether a states set is parallel updatable
  Framework::DataSocketSource< bool > socket_isStatesSetParUpdatable;

  /// socket for the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

};  // class LUSGSSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_LUSGSSetup_hh
