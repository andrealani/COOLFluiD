#ifndef COOLFluiD_Numerics_FluctSplit_UnstP2SUPG_ArtDiffStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_UnstP2SUPG_ArtDiffStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"

#include "MathTools/CFMat.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Physics { namespace NavierStokes { class EulerTerm; } }
  namespace MathTools { class MatrixInverter; }
  

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Artificial Diffusion strategy
 *
 * @author Nadege Villedieu
 *
 */
class UnstP2SUPG_ArtDiffStrategy : public ArtificialDiffusionStrategy {

public: // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  UnstP2SUPG_ArtDiffStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~UnstP2SUPG_ArtDiffStrategy();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ArtificialDiffusionStrategy::configure(args);
  }
void addArtificialDiff(std::vector<RealVector>& residual) ;

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Set up private data and data
   */
  virtual void unsetup();

/**
   * Returns the DataSocket's that this numerical strategy needs as sinks 
   * @return a vector of SafePtr with the DataSockets 
   */ 
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets() 
  { 
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ArtificialDiffusionStrategy::needsSockets(); 
    result.push_back(&socket_pastStates); 
    result.push_back(&socket_interStates); 
    return result; 
  } 

private :
  void doFirstaddArtificialDiff(std::vector<RealVector>& residual);

  void doaddArtificialDiff(std::vector<RealVector>& residual);

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// The socket to use in this strategy for the intermediate states
  Framework::DataSocketSink<Framework::State*> socket_interStates;

  /// temporary storage of the past states
  std::vector<Framework::State*> m_pastStates;

  /// temporary storage of the states
  std::vector<Framework::State*> m_states;

  /// temporary storage of the intermediate states
  std::vector<Framework::State*> m_interStates;

  /// vector of the states that we want to linearise
  std::vector<Framework::State*> m_linstate;

  /// nabla of u at node 0 time n+1
  std::vector<RealVector> m_nabla_u0N1;

  /// nabla of u at node 1 time n+1
  std::vector<RealVector> m_nabla_u1N1;

  /// nabla of u at node 2 time n+1
  std::vector<RealVector> m_nabla_u2N1;

  /// nabla of u at node 0 time n
  std::vector<RealVector> m_nabla_u0n;

  /// nabla of u at node 1 time n
  std::vector<RealVector > m_nabla_u1n;

  /// nabla of u at node 2 time n
  std::vector<RealVector > m_nabla_u2n;

  /// nabla of u at node 0 time n-1
  std::vector<RealVector > m_nabla_u0n1;

  /// nabla of u at node 1 time n-1
  std::vector<RealVector > m_nabla_u1n1;

  /// nabla of u at node 2 time n-1
  std::vector<RealVector > m_nabla_u2n1;

  /// Residual at node 0 
  RealVector m_Res0;

  /// Residual at node 1
  RealVector m_Res1;
  
  /// Residual at node 1
  RealVector m_Res2;

  /// Temporary matrix 
  RealMatrix m_tempmat;

  /// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_kPlus;

  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*  m_inverter;

  /// Sum of the k+
  RealMatrix m_sumKplus;

  /// Scaling of the stabilisation
  RealMatrix m_tau;

  /// adimensionalized normal vector
  RealVector               m_adimNormal;

  ///  eignvalues
  RealVector         m_eValues;

  /// option for adding shock detectyion to the artificial diffusion
  bool m_with_shock_detect;


  
  /// temporary data for holding positive upwind parameter
  RealMatrix                       m_rightEv;

  /// temporary data for holding negative upwind parameter
  RealMatrix                       m_leftEv;

  /// temporary data for holding eignvaluesPlus
  RealVector                       m_eValuesP;

   /// convective term
  Common::SafePtr<Physics::NavierStokes::EulerTerm> _cterm;

  /// Vector of the shock dector
  RealVector m_theta;
  RealVector m_min_states;
  RealVector m_max_states;

  
  /// Storage for the residual of the past ant inter states
  /// at node 0
  std::vector<RealVector> m_PastRes0;

  /// Storage for the residual of the past ant inter states
  /// at node 1
  std::vector<RealVector> m_PastRes1;

  /// Storage for the residual of the past ant inter states
  /// at node 2
  std::vector<RealVector> m_PastRes2;

  /// Number of equations
  CFuint m_nbEqs;

  ///This boolean is used to know if we want a simplier sacling tau:
  ///if it is on false then we use tau=(sumkplus)^-1 otherwise we use tau=alpha*Id
  /// which is less computing demanding
  bool m_simplified_scaling;

};//End class UnstP2SUPG_ArtDiffStrategy

    } // End namespace FluctSplit

}// End namespace COOLFluiD
//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnstP2SUPG_ArtDiffStrategy_hh
