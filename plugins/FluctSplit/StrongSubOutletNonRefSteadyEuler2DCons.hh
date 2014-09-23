#ifndef COOLFluiD_Numerics_FluctSplit_StrongSubOutletNonRefSteadyEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSubOutletNonRefSteadyEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "MathTools/CFMat.hh"

//#include "FluctSplit/SpaceTime_Splitter.hh"

#include "FluctSplit/Splitter.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace MathTools { class MatrixInverter; }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * Multidimensional characteristic bc for characteristic variables for steady simulations
 * Editted from StrongSubOutletNonRefEuler2DCons
 * 
 * @author Gabriel Maher
 * @author Lilla Koloszar
 *
 *
 *
 */

class StrongSubOutletNonRefSteadyEuler2DCons : public FluctuationSplitCom {
public:

  /**
   * Constructor
   */
  StrongSubOutletNonRefSteadyEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSubOutletNonRefSteadyEuler2DCons();

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
  ///Removed Kpast from function
  
  //void computeCharK(std::vector<Framework::State*>& states,   RealVector& _kPlus,   RealVector&  _k,   RealVector&  _kPast, std::vector<RealVector*>& normal, RealVector& faceLength, CFreal& Area);
  
  void computeCharK(std::vector<Framework::State*>& states,   RealVector& _kPlus,   RealVector&  _k, std::vector<RealVector*>& normal, RealVector& faceLength, CFreal& Area);
  
  /// Distribute the residual using LDA scheme
  CFreal distributeLDA(std::vector<Framework::State*>& states, RealVector& m_kPlus, CFreal m_betas, CFuint boundarystate);

  /// temporary variables for the computation of the distribution matrix
  RealVector m_kPlus;
  RealVector m_k;
  CFreal m_sumKplus;

  /// temporary storage of the past upwind matrix
  /// Not necessary for steady calculations
  //RealVector  m_kPast;

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
  /// Not necessary in steady calculations
  //Framework::DataSocketSink< Framework::State*> socket_pastStates;

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
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// elements around the node connectivity on the actual trs
  std::vector< std::vector<CFuint> > bndNod2Elm;

  /// boundary normals belongs to boundary nodes
  std::vector< RealVector > bndNodNorm;

  /// temporary data for holding the nb of variables
  CFuint m_nbEqs;
  
}; // end of class StrongSubOutletNonRefSteadyEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallEuler2DCons_hh
