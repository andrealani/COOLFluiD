#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler3DCons_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler3DCons_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "NavierStokes/Euler3DCons.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler3D
 *
 * @author Thomas Wuilbaut
 *
 */
class UnsteadyWeakSlipWallEuler3DCons : public FluctuationSplitCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  UnsteadyWeakSlipWallEuler3DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyWeakSlipWallEuler3DCons();

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

  /**
   * Compute the normal flux
   */
  void computeNormalFlux(const RealVector& state,
         const RealVector& normal,
         RealVector& flux) const;


  /**
   * Compute the face normal
   */
  void setFaceNormal(const CFuint faceID,
		     RealVector& normal);

 private:

  /// socket for rhs
  Framework::DataSocketSink<
                            CFreal> socket_rhs;

  /// socket for InwardNormals's
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_normals;

  /// socket for Past InwardNormals's
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_pastNormals;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;
  

  /// socket for Past Node's
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;
  
  /// handle to the neighbor cell
  Framework::DataSocketSink<
    Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
  socket_faceNeighCell;
  
  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  
  /// socket for past State's
  Framework::DataSocketSink<
    Framework::State*> socket_pastStates;
  
  /// socket for isUpdated
  Framework::DataSocketSink<
                            bool> socket_isUpdated;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler3DCons> _varSet;

  /// temporary storage for the fluxes
  std::vector<RealVector>     _fluxes;

  /// distribution coefficient
  CFreal                         _alpha;

  /// SubSystem TimeStep
  CFreal                         _dt;

  /// Speed of the face
  RealVector _wallSpeed;

  /// say if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;

}; // end of class UnsteadyWeakSlipWallEuler3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler3DCons_hh
