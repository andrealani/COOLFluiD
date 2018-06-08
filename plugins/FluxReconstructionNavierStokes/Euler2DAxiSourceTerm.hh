#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_Euler2DAxiSourceTerm_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_Euler2DAxiSourceTerm_hh

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
 * A base command for adding the source term for 2D axisymmetric Euler cases
 *
 * @author Ray Vandenhoeck
 *
 *
 */
class Euler2DAxiSourceTerm : public StdSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit Euler2DAxiSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~Euler2DAxiSourceTerm();

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
  void addSourceTerm(RealVector& resUpdates);

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
  
  /// variable for physical data of sol
  RealVector m_solPhysData;

}; // class Euler2DAxiSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_Euler2DAxiSourceTerm_hh

