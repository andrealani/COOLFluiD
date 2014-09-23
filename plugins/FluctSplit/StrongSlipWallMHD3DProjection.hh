#ifndef COOLFluiD_Numerics_FluctSplit_StrongSlipWallMHD3DProjection_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSlipWallMHD3DProjection_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for MHD3D
 *
 * @author Radka Keslerova
 *
 *
 *
 */
class StrongSlipWallMHD3DProjection : public FluctuationSplitCom {

public:

  /**
   * Constructor.
   */
  StrongSlipWallMHD3DProjection(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSlipWallMHD3DProjection();

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
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private:

  /// socket for rhs
  Framework::DataSocketSink<
                            CFreal> socket_rhs;

  /// socket for updateCoef
  Framework::DataSocketSink<
                            CFreal> socket_updateCoeff;

  /// socket for states
  Framework::DataSocketSink<Framework::State*,
			    Framework::GLOBAL> socket_states;
  
  /// socket for isUpdated
  Framework::DataSocketSink<
                            bool> socket_isUpdated;

  /// socket for normals
  Framework::DataSocketSink<InwardNormalsData* > socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > _bcNormals;

  /// flag telling if initialization is performed
  bool _useForInitialization;


}; // end of class StrongSlipWallMHD3DProjection

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallMHD3DProjection_hh
