#ifndef COOLFluiD_Numerics_FluctSplit_StrongSubOutletLinEuler3DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSubOutletLinEuler3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler3DCons.hh"
#include "MathTools/CFMat.hh"

#include "FluctSplit/SpaceTime_Splitter.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace MathTools { class MatrixInverter; }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * Multidimensional characteristic bc for characteristic variables
 *
 * @author Lilla Koloszar
 *
 *
 *
 */

class StrongSubOutletLinEuler3DCons : public FluctuationSplitCom {
public:

  /**
   * Constructor
   */
  StrongSubOutletLinEuler3DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSubOutletLinEuler3DCons();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  void unsetup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure (Config::ConfigArgs& args);

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /// Compute the upwind parameter in char
  void computeCharK(std::vector<Framework::State*>& states,   RealVector& _kPlus,   RealVector&  _k,   RealVector&  _kPast, std::vector<RealVector*>& normal, RealVector& faceLength, CFreal& Area);

  /// Distribute the residual using LDA scheme
  CFreal distributeLDA(RealVector& m_kPlus, CFreal m_betas, CFuint boundarystate);

  /// temporary variables for the computation of the distribution matrix
  RealVector m_kPlus;
  RealVector m_k;
  CFreal m_sumKplus;

  /// temporary storage of the past upwind matrix
  RealVector  m_kPast;

  /// dimensionalized normal vector
  std::vector<RealVector*>   m_Statenormals;

  /// adimensionalized normal vector
  RealVector               m_adimNormal;

  /// adimensionalized normal vector
  RealVector               m_adimCharNormal;

  CFreal m_betas;


private:


/////////////////////// SANDBOX //////////////////////////////////////////////////////////////////////////
  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >    socket_faceNeighCell;

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for the Nodes data
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// socket with flags to check if a state is on the boundary
  Framework::DataSocketSink<bool> socket_isBState;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > _bcNormals;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::LinearizedEuler::LinEuler3DCons> _varSet;

  /// elements around the node connectivity on the actual trs
  std::vector< std::vector<CFuint> > bndNod2Elm;

  /// boundary normals belongs to boundary nodes
  std::vector< RealVector > bndNodNorm;


  /// temporary data for holding the nb of variables
  CFuint m_nbEqs;

}; // end of class StrongSubOutletLinEuler3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallEuler3DCons_hh
