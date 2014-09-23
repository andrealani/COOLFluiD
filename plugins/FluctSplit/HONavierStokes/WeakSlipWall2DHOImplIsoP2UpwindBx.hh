#ifndef COOLFluiD_Numerics_FluctSplit_WeakSlipWall2DHOImplIsoP2UpwindBx_hh
#define COOLFluiD_Numerics_FluctSplit_WeakSlipWall2DHOImplIsoP2UpwindBx_hh

//////////////////////////////////////////////////////////////////////////////



#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

#include "MathTools/CFMat.hh"
#include "FluctSplit/P2Normal.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace Physics { namespace NavierStokes { class EulerTerm; } }
  namespace MathTools { class MatrixInverter; }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an extension of the implicit weak slip wall bc (E. V. Weide)
 * to P2P1 elements (meaning that the face has then 3 nodes)
 * @author Nadege Villedieu
 *
 */
class WeakSlipWall2DHOImplIsoP2UpwindBx : public FluctuationSplitCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  WeakSlipWall2DHOImplIsoP2UpwindBx(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~WeakSlipWall2DHOImplIsoP2UpwindBx();

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
   * Compute the normal flux and the corresponding jacobian
   */
  virtual void computeNormalFluxAndJacob(const Framework::State& state,
					 const RealVector& normal,
					 RealVector& flux,
					 RealMatrix& fluxJacob) = 0;

  /// Compute the transformed states
  /// @pre assumes m_cell has been set to the current geometric entity
  /// @param states the states to be transformed to the consistent states
  /// @return the transformed consistent states
  std::vector<Framework::State*>*
  computeConsistentStates (std::vector<Framework::State*> *const states);

/// Compute the upwind parameter
    ///  kplus and kmin are some output because we want to store them
    ///  For each sub triangle instead of re-computing
  void computeK(const std::vector<Framework::State*>& states,
                std::vector<RealMatrix*>& m_kPlus);
/// Distribute the residual using N scheme
void distributeN(std::vector<RealMatrix*> & m_kPlus,std::vector<RealVector> & phiN);

  /// Distribute the residual using LDA scheme
void distributeLDA(std::vector<RealMatrix*> & m_kPlus,std::vector<RealVector> & phiLDA);


void computeBlendingCoeff(const std::vector<Framework::State*>& states, CFreal & result);

protected:

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the node's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  	Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// physical data
  RealVector _physicalData;

  /// Flux across boundary subface
  RealVector m_subFaceFlux;

  /// Integral of flux Jacobian along boundary subface
  std::vector<RealMatrix> m_subFaceJacob;

  ///Matrices containing contribution from different states
  RealMatrix m_distJacobFrom0;
  RealMatrix m_distJacobFrom1;
  RealMatrix m_distJacobFrom3;
  RealMatrix m_distJacobToState;

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

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im0;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in0;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _im1;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in1;

  /// says if the states are flagged inside this TRS
  std::valarray<bool>   _flagState;

  /// temporary flux jacobian matrix
  RealMatrix _tJacob;

  /// distribution coefficient
//   CFreal m_alpha;

  /// vector of Jacobian fluxes at boundary
//   std::vector<RealMatrix> m_fluxJacob;
  ///temporary variable for integration
  RealMatrix m_qdJacob;

  ///Variables for numerical integration of fluxes
  enum { nbQdPts = 3 };

  RealVector m_subface_center;
  RealVector m_weights;
  RealVector m_qpPos; //Position of quadrature points
                      //Integration domain <-0.25;0.25> considered
  Framework::State* m_qdState;

  RealVector m_qdNormal;
  RealVector m_qdFlux;

  RealVector m_p0;
  RealVector m_p1;
  RealVector m_p3;


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

 /// fluctuation on each sub element
  std::vector<RealVector> m_phiN;

 /// fluctuation on each sub element
  std::vector<RealVector> m_phiLDA;

  RealVector m_sumKplusU;

  RealVector m_sumKU;

  RealMatrix m_sumKplus;

  RealVector m_uTemp;

  RealVector m_uMin;

  RealMatrix m_tmp;


  /// temporary data for computation of upwind parameters for 1st sub-triangle
  std::vector<RealMatrix*> m_kPlus;

/// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_k;

/// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_kMin;

  ///  eignvalues
  std::vector<RealVector*>         m_eValues;

  CFuint m_maxNbStatesInCell;

  /// adimensionalized normal vector
  RealVector               m_adimNormal;

  CFreal m_theta;
  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       m_inverter;

  RealMatrix  m_invK;
/// pmax - pmin
  CFreal _deltaP;

  /// reference length
  CFreal _length;

  /// reference speed
  CFreal _speed;

  /// name of the variable chosen for blending coefficient
  std::string _varName;

  /// temporary physical data array
  RealVector _pData;

  /// use shock-detector on all sub elements
  bool m_use_max_theta;

/// convective term
 Common::SafePtr<Physics::NavierStokes::EulerTerm> _cterm;

/// flag to see if we are dealing with an high order scheme
  bool _isHO;
  /// max number of subcells
  CFuint m_max_nbsubcells;

  /// Shock detector that we want to use
  std::string _sh_detector;

  RealVector  m_uInflow;

  /// temporary gradient
  RealVector _grad;

  /// ID of the variable chosen for blending coefficient
  CFuint _varID;

 /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
}; // end of class WeakSlipWall2DHOImplIsoP2UpwindBx


//////////////////////////////////////////////////////////////////////////////

  ///This method is in fact defined in FluctuationSplitStrategy.hh
  ///Since this is a command and not a strategy ...
  ///... we copy and paste the function

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer       VectorTransformer;


  inline std::vector<Framework::State*> *
  WeakSlipWall2DHOImplIsoP2UpwindBx::computeConsistentStates(std::vector<Framework::State*> *const states)
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

#endif // COOLFluiD_Numerics_FluctSplit_WeakSlipWall2DImpl_hh

