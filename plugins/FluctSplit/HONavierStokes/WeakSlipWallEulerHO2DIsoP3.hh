#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP3_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP3_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "FluctSplit/P3Normal.hh"

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
 * for P3P3 elements, it is the extention of the bc of E. Van der Weide
 * to high order discretization
 *
 * @author Nadege Villedieu
 * @author Martin Vymazal
 *
 */

class WeakSlipWallEulerHO2DIsoP3 : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallEulerHO2DIsoP3(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEulerHO2DIsoP3();

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
  std::vector<RealVector> m_subFaceFlux;

  /// distribution coefficient
  CFreal     _alpha;

// vector of fluxes
    std::vector<RealVector> m_nodeFlux;

  /// vector of states of the face
  std::vector<Framework::State*> m_face_states;

  /// vector of nodes of the face
  std::vector<Framework::Node*> m_face_nodes;

  ///Functor to compute normals of P3P3 triangle on the fly
  P3Normal m_CP3N;

  ///Variables for numerical integration of fluxes
  RealVector m_weights;
  RealVector m_qpPos; //Position of quadrature points
                      //Integration domain <-1/3;1/3> considered
  Framework::State* m_qdState;

  RealVector m_qdNormal;
  RealVector m_qdFlux;

  ///Number of integration points
  enum { nbQdPt = 3 };

}; // end of class WeakSlipWallEulerIsoP3

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP3_hh
