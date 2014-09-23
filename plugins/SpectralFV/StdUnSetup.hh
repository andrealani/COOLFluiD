#ifndef COOLFluiD_SpectralFV_StdUnSetup_hh
#define COOLFluiD_SpectralFV_StdUnSetup_hh

#include "MathTools/RealMatrix.hh"
#include "SpectralFV/SpectralFVMethodData.hh"
#include "SpectralFV/SpectralFVElementData.hh"

namespace COOLFluiD {
  namespace SpectralFV {

/// This is a standard command to deallocate data specific to an empty method
class StdUnSetup : public SpectralFVMethodCom {

public:

  /// Constructor
  explicit StdUnSetup(const std::string& name);

  /// Destructor
  ~StdUnSetup() {}

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

  /// socket for face normal transformation matrices
  /// (this only applies to elements with a linear transformation to a reference element!!!)
  Framework::DataSocketSink<RealMatrix > socket_faceNormTransfMatrices;

  /// socket for state volumes
  Framework::DataSocketSink< CFreal > socket_volumes;

  /// socket for face surfaces
  Framework::DataSocketSink< CFreal > socket_faceSurf;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

}; // class StdUnSetup

  }  // namespace SpectralFV
}  // namespace COOLFluiD

#endif // COOLFluiD_SpectralFV_StdUnSetup_hh
