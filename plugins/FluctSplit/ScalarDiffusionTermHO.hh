#ifndef COOLFluiD_Numerics_FluctSplit_ScalarDiffusionTermHO_hh
#define COOLFluiD_Numerics_FluctSplit_ScalarDiffusionTermHO_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace LinearAdv {
      class AdvectionDiffusionVarSet;
      class LinearAdv2DVarSet;
    }    
  }
  


    namespace FluctSplit {
            
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model
 *
 * @author Andrea Lani
 *
 */
class ScalarDiffusionTermHO : public ComputeDiffusiveTerm {
public:

  /**
   * Constructor
   */
  ScalarDiffusionTermHO(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~ScalarDiffusionTermHO();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  void setup();
  
  /**
   * Set the update variable set
   */
  void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar);
   
  /**
   * Compute the diffusive term flux in the current cell
   */
  void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
			    std::vector<RealVector>& result,
			    const bool updateCoeff);
  /**
    * Set the update variable set
    */
  void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet>); 
// Integral of diffusive part times basis function

  void fluctuation_diff_galerkin(CFuint&, CFuint&, CFuint&);

// Integral of diffusive part times bubblefunction
  void fluctuation_diff_bubble(CFuint&, CFuint&, CFuint&);

   void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob){
    throw Common::NotImplementedException (FromHere(),"ScalarDiffTermHO::computePicardJacob()");
  }
protected: // data
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::LinearAdv::AdvectionDiffusionVarSet> _diffVar;
  
  // acquaintance of the update convective var set
  Common::SafePtr<Physics::LinearAdv::LinearAdv2DVarSet> _updateVar;

  // array of cell states
  std::vector<Framework::State*> * m_cellStates;

  // array of cell states
  std::vector<RealVector*> _states;

  // array of values (rho, u, v, T)
  RealMatrix _values;
  
  // array of gradients
  std::vector<RealVector*> _gradients;
  
  // array of average values (rho, u, v, T)
  RealVector _avValues;
  
  // unit normal
  RealVector _normal;
 
  /// fluctuation from diffusion part bubble part
  RealVector* _phi_diff_bub;

 /// fluctuation from diffusion part galerkin
  std::vector<RealVector*> _phi_diff_gal;

/// fluctuation from diffusion part bubble part splitted with the beta of LDA normally
   std::vector<RealVector*> _phi_diff_bub_split;
 
/// Pointer to the current cell being processed.
  CFuint m_cellID;

  /// quadrature points per face
  RealVector qd0;
  RealVector qd1;
  RealVector qd2;
  RealVector wqd;
/// states at quadrature points
  RealVector qdstates;
  

  /// Coeficients used to copute the diffusion residual
  MathTools::CFMat<CFreal> kappa;

// volume of the cell
  CFreal _cellVolume;
  

}; // end of class ScalarDiffusionTermHO
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ScalarDiffusionTermHO_hh
