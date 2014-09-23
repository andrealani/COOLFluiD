#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP2_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP2_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "FluctSplit/P2Normal.hh"

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
 * This class represents a weak slip wall bc for Euler2D
 * for P2P2 elements, it is the extention of the bc of E. Van der Weide
 * to high order discretization
 *
 * @author Nadege Villedieu
 * @author Martin Vymazal
 *
 */

class WeakSlipWallEulerHO2DIsoP2 : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallEulerHO2DIsoP2(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEulerHO2DIsoP2();

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

  /**
   * Compute the normal flux
   */
  void computeNormalFlux(const Framework::State& state,
			 const RealVector& normal,
			 RealVector& flux);

  /**
   * Compute the face normal
   */
  void setFaceNormal(const CFuint faceID,
		     RealVector& normal);

 private:

  // the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  // the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// say if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;

  /// Euler var set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;

  /// physical data
  RealVector _physicalData;

  /// Total flux across boundary face
  RealVector m_totFlux0;
  RealVector m_totFlux1;
  RealVector m_totFlux;

  /// distribution coefficient
  CFreal     _alpha;

// vector of fluxes
    std::vector<RealVector> m_flux;

  /// vector of states of the face
  std::vector<Framework::State*> m_states;

  /// vector of nodes of the face
  std::vector<Framework::Node*> m_nodes;

  ///Functor to compute normals of P2P2 triangle on the fly
  P2Normal m_CP2N;

  ///Variables for numerical integration of fluxes
  RealVector m_weights;
  RealVector m_qpPos; //Position of quadrature points
                      //Integration domain <-0.25;0.25> considered
  Framework::State* m_qdState;

  RealVector m_qdNormal;
  RealVector m_qdFlux;

}; // end of class WeakSlipWallEuler2DCurved

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP2_hh
