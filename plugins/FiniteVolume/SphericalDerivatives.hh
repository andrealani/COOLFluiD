#ifndef COOLFluiD_Numerics_FiniteVolume_SphericalDerivatives_hh
#define COOLFluiD_Numerics_FiniteVolume_SphericalDerivatives_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////

/**
 * This class computes spherical from cartesian derivatives
 *
 * @author Andrea Lani
 *
 */
class SphericalDerivatives : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  SphericalDerivatives(const std::string& name);

  /**
   * Default destructor
   */
  ~SphericalDerivatives();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

/**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: //data

  /// derivatives in r
  Framework::DataSocketSource<CFreal> socket_uR;
  
  /// derivatives in theta
  Framework::DataSocketSource<CFreal> socket_uTheta;
  
  /// derivatives in phi
  Framework::DataSocketSource<CFreal> socket_uPhi;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;

  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uZ;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
}; // end of class SphericalDerivatives

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeMultiFluidMHD_SphericalDerivatives_hh



