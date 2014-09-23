#ifndef COOLFluiD_Numerics_FluctSplit_StrongMirrorEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongMirrorEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for MHD2D
 *
 * @author Andrea Lani
 *
 */

class StrongMirrorEuler2DCons : public FluctuationSplitCom {

public:

  /**
   * Constructor.
   */
  StrongMirrorEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongMirrorEuler2DCons();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

protected:

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// socket for normals
  Framework::DataSocketSink<InwardNormalsData* > socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink <
    Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > m_bcNormals;

}; // end of class StrongMirrorEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongMirrorEuler2DCons_hh
