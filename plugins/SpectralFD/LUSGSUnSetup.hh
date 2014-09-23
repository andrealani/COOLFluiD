#ifndef COOLFluiD_SpectralFD_LUSGSUnSetup_hh
#define COOLFluiD_SpectralFD_LUSGSUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

 //////////////////////////////////////////////////////////////////////////////

/// This is a standard command to deallocate data specific to an empty method
class LUSGSUnSetup : public StdUnSetup {

public:

  /// Constructor
  explicit LUSGSUnSetup(const std::string& name);

  /// Destructor
  ~LUSGSUnSetup() {}

  /// Execute processing actions
  void execute();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // protected data

/// socket for diagonal block Jacobian matrices
Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

/// socket for rhs of current set of states
Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

/// handle to the IDs of the states in each set of states
Framework::DataSocketSink< std::vector< CFuint > > socket_statesSetStateIDs;

//////////////////////////////////////////////////////////////////////////////

}; // class LUSGSUnSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_LUSGSUnSetup_hh
