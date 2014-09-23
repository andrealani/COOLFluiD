#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP2Upwind_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP2Upwind_hh

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

class WeakSlipWallEulerHO2DIsoP2Upwind : public FluctuationSplitCom {

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWallEulerHO2DIsoP2Upwind(const std::string& name);

  /**
   * Default destructor
   */
  ~WeakSlipWallEulerHO2DIsoP2Upwind();

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


 private:

  /// Compute the transformed states
  /// @pre assumes m_cell has been set to the current geometric entity
  /// @param states the states to be transformed to the consistent states
  /// @return the transformed consistent states
  std::vector<Framework::State*>*
  computeConsistentStates (std::vector<Framework::State*> *const states);

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

  /// Flux across boundary subface

  RealVector m_subFaceFlux;

// vector of fluxes
    std::vector<RealVector> m_flux;

  /// vector of local face indexes of the big triangle
  static CFuint m_subFaceIdx0[6];
  static CFuint m_subFaceIdx1[6];
  static CFuint m_subFaceIdx2[6];

  CFuint* m_subFaceIdx;

  /// vector of local triangle indexes
  static CFuint m_subElemIdx0[2];
  static CFuint m_subElemIdx1[2];
  static CFuint m_subElemIdx2[2];

  CFuint* m_subElemIdx;

  /// vectors with definitions of subtriangles (list of vertices of subtriangles)

  static CFuint m_subTri0[6];
  static CFuint m_subTri1[6];
  static CFuint m_subTri2[6];

  CFuint* m_subTri;

  /// vector of nodes of the face
  std::vector<Framework::Node*> m_face_nodes;

  /// vector of states of the face
  std::vector<Framework::State*> m_face_states;

  ///Functor to compute normals of P2P2 triangle on the fly
  P2Normal m_CP2N;

  ///Variables for numerical integration of fluxes
  enum { nbQdPts = 3 };

  RealVector m_weights;
  RealVector m_qpPos; //Position of quadrature points
                      //Integration domain <-0.25;0.25> considered
  RealVector m_subface_center;
  Framework::State* m_qdState;

  RealVector m_qdNormal;
  RealVector m_qdFlux;

  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;

   /// solution variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_solutionVar;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;

  /// states that compose each sub element
  std::vector<Framework::State*> m_subStates;

  /// residuals for each sub element
  std::vector<RealVector> m_subresidual;

  /// temporary subcell normals
  InwardNormalsData * m_subcell_normals;

  RealMatrix matrix_face_norms;
  RealMatrix matrix_node_norms;

  RealVector vector_face_areas;
  RealVector vector_node_areas;

  RealMatrix m_betasInSubTriag;

  /// transformed linear states in current cell
  std::vector<Framework::State*> * m_linearStates;

  /// temporary storage of the cell residual
  std::vector<RealVector> m_residual;

}; // end of class WeakSlipWallEulerHO2DIsoP2Upwind


  ///This method is in fact defined in FluctuationSplitStrategy.hh
  ///Since this is a command and not a strategy ...
  ///... we copy and paste the function

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer       VectorTransformer;


  inline std::vector<Framework::State*> *
  WeakSlipWallEulerHO2DIsoP2Upwind::computeConsistentStates(std::vector<Framework::State*> *const states)
  {
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  FluctuationSplitData& mdata = getMethodData();

  SafePtr<JacobianLinearizer> jacob_lin  = mdata.getLinearizer();
  SafePtr<VectorTransformer>  vec_transf = mdata.getUpdateToLinearVecTrans();
  SafePtr<MatrixTransformer>  mat_transf = mdata.getLinearToDistribMatTrans();

  // transform the update states to linear varset
  m_linearStates = vec_transf->transform(states);

  // computes an average state in which jacobians will be linearized
  jacob_lin->linearize(*m_linearStates);

  // transformation from linear to consistent variables
  return mat_transf->transformFromRef(m_linearStates);
  }





//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWallEulerHO2DIsoP2Upwind_hh
