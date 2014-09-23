#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2D_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace FluctSplit {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strong slip wall bc for Euler2D
 *
 * @author Andrea Lani
 *
 */
class WeakSlipWallEuler2D : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallEuler2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~WeakSlipWallEuler2D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:

  /**
   * Execute on a set of dofs
   */
  virtual void executeOnTrs();

  /**
   * Compute the normal flux
   */
  virtual void computeNormalFlux(const Framework::State& state,
				 const RealVector& normal,
				 RealVector& flux);
  
  /**
   * Compute the face normal
   */
  virtual void setFaceNormal(const CFuint faceID, RealVector& normal);
  
protected:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
  socket_faceNeighCell;
  
  /// say if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;
  
  /// Euler var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;
  
  /// physical data
  RealVector _physicalData;
  
  /// temporaray flux
  RealVector _tmpFlux;
  
  /// distribution coefficient
  CFreal     _alpha;
  
  // vector of fluxes of each state
  std::vector<RealVector> m_flux;
  
  // vector of states of the face
  std::vector<Framework::State*> m_states;
  
}; // end of class WeakSlipWallEuler2D
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEuler2D_hh
