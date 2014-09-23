#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingCorrFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingCorrFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/StegerWarmingFlux.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the StegerWarmingCorr flux
 *
 * @author Andrea Lani
 *
 */
class StegerWarmingCorrFlux : public StegerWarmingFlux {
public:

  /**
   * Constructor
   */
  StegerWarmingCorrFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~StegerWarmingCorrFlux();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Build the face data
   */
  virtual void buildFaceBCData();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);
  
protected:
  
  /**
   * Compute the dissipation coefficient
   */
  void computeEps(RealVector& pData,
		  const CFreal& machL,
		  const CFreal& machR, 
		  CFuint& countSubsonic,
		  CFuint& countSupersonic,
		  CFreal& eps);
  
  /**
   * Compute upwind polynomial reconstruction
   */
  void upwindReconstruct();
  
  /**
   * Set the stencil for the reconstruction
   */
  void setStencil();

protected:
  
  /// storage of the States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

   /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// constant for 1./6.
  const CFreal _sixth;
  
  /// constant for 1./3.
  const CFreal _third;
  
  /// constant 
  const CFreal _k1;
  
  /// average left state
  Framework::State* _avStateL;
  
  /// average right state
  Framework::State* _avStateR;  
  
  /// array of flags for faces parallel to the wall
  std::vector<bool> _flagNormalFace;
  
  /// stencil storage
  std::vector<std::vector<Framework::State*> > _stencil;
  
  // update states array
  std::vector<Framework::State*> _upStates; 
  
  // wall temperature
  CFreal _TWall;
  
  // coefficient to control the jacobian dissipation
  CFreal _jacobDissip;
  
  /// coefficient to control the pressure gradient based weights
  CFreal _sigma;
  
  /// coefficient to control the pressure gradient based limiter
  CFreal _sigmaLim;
  
  /// maximum number of normal faces to neglect for carbuncle fix
  CFuint _maxNbNormalFaces;
    
  /// flag telling if the scheme must be limited
  bool _useLimiter;
  
  /// flag telling if the upwind reconstruction
  bool _useUpwindPolyRec;
  
  /// names of the wall TRSs
  std::vector<std::string> _wallTrsNames;
    
}; // end of class StegerWarmingCorrFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingCorrFlux_hh
