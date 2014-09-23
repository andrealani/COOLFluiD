#ifndef COOLFluiD_Numerics_FluctSplit_StrongFarFieldNonRefEuler2DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongFarFieldNonRefEuler2DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/Euler2DCons.hh"
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

class StrongFarFieldNonRefEuler2DCons : public FluctuationSplitCom {
public:

  /**
   * Constructor
   */
  StrongFarFieldNonRefEuler2DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongFarFieldNonRefEuler2DCons();

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
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

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
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> _varSet;

  /// elements around the node connectivity on the actual trs
  std::vector< std::vector<CFuint> > bndNod2Elm;

  /// boundary normals belongs to boundary nodes
  std::vector< RealVector > bndNodNorm;

  /// temporary data for holding the nb of variables
  CFuint m_nbEqs;
  
    /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n}\f$
  RealVector _r1;

  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} \f$
  RealVector _r2;

   /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} + a \f$
  RealVector _r3;
  
    /// temporry for the variables
  RealVector m_var_values;

  /// a vector of string to hold the functions
  std::vector<std::string> m_function_inflow;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars_inflow;

  /// function for the mean flow
  Framework::VectorialFunction m_function_parser_inflow;
  
  std::vector<RealVector> inflow;
  
  
  
}; // end of class StrongFarFieldNonRefEuler2DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallEuler2DCons_hh
