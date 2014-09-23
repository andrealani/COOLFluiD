#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2DConsHOImpl_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2DConsHOImpl_hh

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
class UnsteadyWeakSlipWallEuler2DConsHOImpl : public FluctuationSplitCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  UnsteadyWeakSlipWallEuler2DConsHOImpl(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadyWeakSlipWallEuler2DConsHOImpl();

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

  /**
   * Compute the normal flux and the corresponding jacobian
   */
  void computeNormalFluxAndJacob(const RealVector& state,
         const RealVector& normal,
         RealVector& flux,
         RealMatrix& fluxJacob) const;


  /**
   * Compute the face normal
   */
  void setFaceNormal(const CFuint faceID,
		     RealVector& normal);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Compute the modification of the BC for the first iteration
   */
  void docomputeFirstBC();

/**
   * Compute the modification of the BC for all each iterations except the first
   */
  void docomputeBC();

 private:    
             
  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
             
  // the socket to the data handle of the state's
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// socket for past State's
  Framework::DataSocketSink<
                            Framework::State*> socket_interStates;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  // the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// distribution coefficient
  CFreal                         m_alpha;

  /// SubSystem TimeStep
  CFreal                         _dt;

  /// Speed of the face
  RealVector _wallSpeed;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im0;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in0;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im1;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in1;

  /// say if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;

  std::vector<RealVector> m_fluxn1;
  std::vector<RealVector> m_fluxn12;
  std::vector<RealVector> m_fluxn;
  std::vector<RealVector> m_flux;
  std::vector<RealMatrix> m_fluxJacobn1;
  std::vector<RealMatrix> m_fluxJacobn12;
  std::vector<RealMatrix> m_fluxJacobn;
  std::vector<RealMatrix> m_fluxJacob;

}; // end of class UnsteadyWeakSlipWallEuler2DConsHOImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyWeakSlipWallEuler2DConsHOImpl_hh
