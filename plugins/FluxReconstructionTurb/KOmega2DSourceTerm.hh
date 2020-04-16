#ifndef COOLFluiD_FluxReconstructionMethod_KOmega2DSourceTerm_hh
#define COOLFluiD_FluxReconstructionMethod_KOmega2DSourceTerm_hh

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
 * A base command for adding the source term for 2D K omega
 *
 * @author Ray Vandenhoeck
 */
class KOmega2DSourceTerm : public StdSourceTerm {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit KOmega2DSourceTerm(const std::string& name);

  /**
   * Destructor.
   */
  ~KOmega2DSourceTerm();

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
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

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
  
  /**
   * compute the production term
   */
  virtual  void computeProductionTerm(const CFuint iState,
				      const CFreal& PcoFactor,const CFreal& MUT,
				      CFreal& KProdTerm,
				      CFreal& OmegaProdTerm);
  
  /**
   * compute the destruction term
   */
  virtual void  computeDestructionTerm(const CFuint iState,
				       const CFreal& DcoFactor, CFreal& K_desterm, 
				       CFreal& Omega_desterm);


  static std::string getModuleName(); 

protected: // data
  
  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// handle to the wall distance
  Framework::DataSocketSink<CFreal> socket_wallDistance;
  
  /// storage for the wall shear stress velocity
  Framework::DataSocketSource<CFreal> socket_wallShearStressVelocity;

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
  
  /// the gradients in the current cell
  std::vector< std::vector< RealVector* > > m_cellGrads;
  
  // k production term
  CFreal  m_prodTerm_k;
  
  // omega production term
  CFreal  m_prodTerm_Omega;
  
  // omega destruction term
  CFreal  m_destructionTerm_Omega;
  
  // k destruction term
  CFreal  m_destructionTerm_k;
  
  // wall distance in current sol point
  RealVector m_currWallDist; 
  
  // boolean telling whether the case is axisymmetric
  bool m_isAxisymmetric;
  
  // boolean telling whether the source terms need to be limited for stability
  bool m_limitP;

  
}; // class KOmega2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_KOmega2DSourceTerm_hh

