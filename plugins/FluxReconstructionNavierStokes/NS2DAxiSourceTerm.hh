#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_NS2DAxiSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_NS2DAxiSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/StdSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * A base command for adding the source term for 2D axisymmetric NS cases
 *
 * @author Ray Vandenhoeck
 *
 *
 */
class NS2DAxiSourceTerm : public StdSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit NS2DAxiSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~NS2DAxiSourceTerm();

  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
      needsSockets();

protected:

  /**
   * get data required for source term computation
   */
  void getSourceTermData();

  /**
   * add the source term
   */
  void addSourceTerm();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  /// the source term for one state
  RealVector m_srcTerm;

  /// dimensionality
  CFuint m_dim;
  
  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;
  
  /// corresponding diffusive variable set
  Common::SafePtr<Framework::DiffusiveVarSet> m_diffVarSet;
  
  /// variable for physical data of sol
  RealVector m_solPhysData;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> m_dummyGradients;
  
  /// the gradients in the neighbouring cell
  std::vector< std::vector< RealVector >* > m_cellGrads;
  
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

}; // class NS2DAxiSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_NS2DAxiSourceTerm_hh

