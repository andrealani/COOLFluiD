#ifndef COOLFluiD_Numerics_FluctSplit_WeakBC2D_hh
#define COOLFluiD_Numerics_FluctSplit_WeakBC2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "WeakBC.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler2D
/// @author Andrea Lani
class FluctSplit_API WeakBC2D : public WeakBC {

public:

  /// Constructor.
  WeakBC2D(const std::string& name);

  /// Default destructor
  virtual ~WeakBC2D();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// UnSet up private data and data of the aggregated classes
  /// in this command after processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Execute on a set of dofs
  void executeOnTrs();

  /// Set the state vector in the ghost State's
  virtual void setGhostState(const Framework::State& state,
                             Framework::State& gstate) = 0;

  /// Set the average normal if moving boundary
  void setMovingFaceNormal(const CFuint faceID);

protected:

  /// socket for Past InwardNormals's
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_pastNormals;

  /// socket for Past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;

  /// socket for isUpdated
  Framework::DataSocketSink<
                            bool> socket_isUpdated;

  /// temporary states

  Framework::State* m_gstate;

  Framework::State* m_state;

  std::vector<RealVector> m_flux;

}; // end of class WeakBC2D

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakBC2D_hh
