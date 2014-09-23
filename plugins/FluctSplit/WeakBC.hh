#ifndef COOLFluiD_Numerics_FluctSplit_WeakBC_hh
#define COOLFluiD_Numerics_FluctSplit_WeakBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
// #include "FluctSplit/P2Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler2D
/// @author Andrea Lani
class FluctSplit_API WeakBC : public FluctuationSplitCom {
public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer VectorTransformer;

  /// Constructor.
  WeakBC(const std::string& name);

  /// Default destructor
  virtual ~WeakBC();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// Execute on a set of dofs
  virtual void executeOnTrs() = 0;

  /// Compute the face normal
  void setFaceNormal(const CFuint faceID);

  /// Compute the additional flux to distribute to the boundary
  /// states
  /// @param states    state vector in the "ghost" cell
  /// @param residual  residual contribution for the ghost state
  void computeFlux(const std::vector<Framework::State*>& states,
                   RealVector& flux);

protected: // member data

  /// handle to the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

  /// handle to the normals
  Framework::DataSocketSink< InwardNormalsData*>  socket_normals;

  /// handle to the update coefficient
  Framework::DataSocketSink< CFreal>  socket_updateCoeff;

  /// handle to the neighbor cell
  Framework::DataSocketSink<
  Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        socket_faceNeighCell;

  /// handle to the states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL>  socket_states;

/// Transformer from Distribution to Solution Variables
  Common::SafePtr<MatrixTransformer> m_solutionToDistMatTrans;

  /// Transformer from Solution to Distribution Variables
  Common::SafePtr<MatrixTransformer> m_distToSolutionMatTrans;

  /// Transformer from Solution to Distribution Variables
  Common::SafePtr<MatrixTransformer> m_linearToDistMatTrans;

  /// Transformer from Update to Linear Variables
  Common::SafePtr<VectorTransformer> m_updateToLinearVecTrans;

  /// jacobian linearizer
  Common::SafePtr<Framework::JacobianLinearizer> m_linearizer;

  /// jacobian linearizer
  Common::SafePtr<Framework::ConvectiveVarSet> m_distribVar;

  /// storage for the temporary face normal
  RealVector m_faceNormal;

  /// storage for the temporary adimensional face normal
  RealVector m_adimNormal;

  /// matrix of right eigenvectors
  RealMatrix                       m_rightEv;

  /// matrix of left eigenvectors
  RealMatrix                       m_leftEv;

  /// array of the eigenvalues
  RealVector                       m_eValues;

  /// array of the eigenvalues plus
  RealVector                       m_eValuesP;

  /// kplus matrix
  RealMatrix                       m_kPlus;

  /// array for a temporary flux
  RealVector     m_tFlux;

  /// temporary storage for two states
  std::vector<Framework::State*>  m_twoStates;

  /// say if the states are flagged inside this TRS
  std::valarray<bool>          m_flagState;

  /// distribution coefficient
  CFreal                         m_alpha;

  /// is boundary composed of high-order faces?
  CFuint                            m_faceGeoOrder;

  /// temporary storage of transformed states
  std::vector<RealVector> m_tStates;

  /// linearized states
  std::vector<RealVector*> m_zStates;

  ///HACK FOR HIGH ORDER FACES:
  /// vector of nodes of the face
//   std::vector<Framework::Node*> m_nodes;

  ///Functor to compute normals of P2P2 triangle on the fly
//   P2Normal m_CP2N;

//   RealVector m_hoFaceNormal0;
//   RealVector m_hoFaceNormal1;
//   RealVector m_hoFaceNormal2;

//   CFreal J0, J1, J2;
//   CFreal totFaceLen;

  /// Order of the solution approximation
  CFuint m_solorder;

}; // end of class WeakBC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_WeakBC_hh
