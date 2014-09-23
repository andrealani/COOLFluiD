#ifndef COOLFluiD_Numerics_FluctSplit_StrongSubInletHedstrom3DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSubInletHedstrom3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "LinEuler/LinEuler3DCons.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a subsonic inlet characteristic-based
 * boundary condition for RD schemes for 3D LEE
 *
 * @author Lilla Koloszar
 *
 */

class StrongSubInletHedstrom3DCons : public FluctuationSplitCom {
public:

  /**
   * Constructor
   */
  StrongSubInletHedstrom3DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSubInletHedstrom3DCons();

   /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );


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

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// socket for normals
  Framework::DataSocketSink<InwardNormalsData* > socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::LinearizedEuler::LinEuler3DCons> _varSet;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > _bcNormals;

  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n}\f$
  RealVector _r1;

  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} \f$
  RealVector _r2;

   /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} \f$
  RealVector _r3;
  
  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} - a \f$
  RealVector _r5;

}; // end of class StrongSubInletHedstrom3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
