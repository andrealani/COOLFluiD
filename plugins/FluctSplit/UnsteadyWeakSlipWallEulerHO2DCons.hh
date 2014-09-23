#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEulerHO2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEulerHO2DCons_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "Framework/State.hh"
#include "Framework/Node.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler2D
 *
 * @author Thomas Wuilbaut
 *
 */
class UnsteadyWeakSlipWallEulerHO2DCons : public FluctuationSplitCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  UnsteadyWeakSlipWallEulerHO2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyWeakSlipWallEulerHO2DCons();

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

  /**
   * Compute the modification of the BC for the first iteration
   */
  void docomputeFirstBC();

/**
   * Compute the modification of the BC for all each iterations except the first
   */
  void docomputeBC();

 private:

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// socket for rhs
  Framework::DataSocketSink<
                            CFreal> socket_rhs;

  /// socket for InwardNormals's
  Framework::DataSocketSink<
                            InwardNormalsData*> socket_normals;

  /// socket for Node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;


  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// socket for past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_pastStates;

  /// socket for past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_interStates;

  /// socket for isUpdated
  Framework::DataSocketSink<
                            bool> socket_isUpdated;

  /// distribution coefficient
  CFreal                         _alpha;

  /// SubSystem TimeStep
  CFreal                         _dt;

  /// Speed of the face
  RealVector _wallSpeed;

  /// say if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;


  std::vector<RealVector> m_fluxn1;
  std::vector<RealVector> m_fluxn12;
  std::vector<RealVector> m_fluxn;
  std::vector<RealVector> m_flux;

}; // end of class UnsteadyWeakSlipWallEulerHO2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEulerHO2DCons_hh
