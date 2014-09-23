#ifndef COOLFluiD_Numerics_FluctSplitNEQ_TCNEQDiffTerm_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_TCNEQDiffTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
   
  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }    
  }

  namespace FluctSplitNEQ {
      
//////////////////////////////////////////////////////////////////////////////
      
/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model with thermo-chemical NEQ
 *
 * @author Andrea Lani
 *
 */
template <typename BASE, typename UPDATEVAR>
class TCNEQDiffTerm : public BASE {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  TCNEQDiffTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~TCNEQDiffTerm();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Set the update variable set
   */
  void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar);
   
  /**
   * Compute the diffusive term flux in the current cell
   */
  virtual void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
				    std::vector<RealVector>& result,
				    const bool updateCoeff);
  /**
    * Set the update variable set
    */
  void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet>); 
  
  /**
   * Compute the cell gradient and average state and put them
   * into @see DistributionData
   */
  virtual void computeCellGradientsAndAverageState
  (Framework::GeometricEntity *const geo, const RealVector& pdata);
    
protected: // data
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVar;
  
  // average radius in the cell
  CFreal _radius;
    
  // acquaintance of the update convective var set
  Common::SafePtr<UPDATEVAR> _updateVar;
    
  // array of cell states
  std::vector<RealVector*> _states;
  
  // array of values (rho, u, v, T)
  RealMatrix _values;
  
  // unit normal
  RealVector _normal;

  /// transformation matrix from rho_s to y_s
  RealMatrix _jacobRiYi;
  
  /// gradients of rho_s
  RealMatrix _gradRi;
    
  /// gradients of y_s
  RealMatrix _gradYi;
    
}; // end of class TCNEQDiffTerm
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "TCNEQDiffTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_TCNEQDiffTerm_hh
