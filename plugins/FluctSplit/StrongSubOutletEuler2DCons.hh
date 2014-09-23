#ifndef COOLFluiD_Numerics_FluctSplit_StrongSubOutletEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSubOutletEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DCons;
    }
  }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a subsonic outlet characteristic-based
 * boundary condition for RD schemes
 *
 * @author Andrea Lani
 *
 *
 *
 */
class StrongSubOutletEuler2DCons : public FluctuationSplitCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  StrongSubOutletEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSubOutletEuler2DCons();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

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

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for normals
  Framework::DataSocketSink<InwardNormalsData* > socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n}+ a \f$
  RealVector                               _r3;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > _bcNormals;

  /// array to store the temporary computation of the derivative
  /// of pressure with respect of the conservative variables
  RealVector                               _dPdU;

  /// static pressure
  CFreal                                   _pressure;

}; // end of class StrongSubOutletEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallEuler2DCons_hh
