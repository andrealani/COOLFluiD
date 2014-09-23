#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingFluxT_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingFluxT_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the StegerWarming flux
 *
 * @author Andrea Lani
 *
 */
template <typename BASE>
class StegerWarmingFluxT : public BASE {
public:
  
  /** 
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  StegerWarmingFluxT(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~StegerWarmingFluxT();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _solutionVarSet;
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _updateVarSet;
  
  /// variable set transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> _updateToSolutionVarTrans;
  
  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;
  
  /// eigenvalues
  RealVector _eValues;
  
  /// matrix of right eigenvectors
  RealMatrix   _jacobPlus;

  /// matrix of left eigenvectors
  RealMatrix   _jacobMin;
  
  /// dummy matrix of eigenvectors
  RealMatrix   _jacobDummy;
    
  /// temporary unit normal
  RealVector   _tempUnitNormal;

  /// number of equations
  CFuint m_nbEqs;
  
  /// start ID for equations
  CFuint m_startID;
  
}; // end of class StegerWarmingFluxT

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/StegerWarmingFluxT.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingFluxT_hh
