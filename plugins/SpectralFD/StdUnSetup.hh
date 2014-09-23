#ifndef COOLFluiD_SpectralFD_StdUnSetup_hh
#define COOLFluiD_SpectralFD_StdUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

 //////////////////////////////////////////////////////////////////////////////

/// This is a standard command to deallocate data specific to an empty method
class StdUnSetup : public SpectralFDMethodCom {

public:

  /// Constructor
  explicit StdUnSetup(const std::string& name);

  /// Destructor
  ~StdUnSetup() {}

  /// Execute processing actions
  virtual void execute();

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
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected: // protected data

  /// socket for state volumes
  Framework::DataSocketSink< CFreal > socket_volumes;

  /// socket for unit normals in face flux points
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

//////////////////////////////////////////////////////////////////////////////

}; // class StdUnSetup

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_StdUnSetup_hh
