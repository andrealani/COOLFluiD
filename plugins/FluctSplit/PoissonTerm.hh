#ifndef COOLFluiD_Numerics_FluctSplit_PoissonTerm_hh
#define COOLFluiD_Numerics_FluctSplit_PoissonTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace Poisson {
      class PoissonConvVarSet;
      class PoissonDiffVarSet;
    }
  }
  
  namespace FluctSplit {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to Poisson equation
 *
 * @author Andrea Lani
 *
 */
class PoissonTerm : public ComputeDiffusiveTerm {
public:

  /**
   * Constructor
   */
  PoissonTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~PoissonTerm();
  
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
  virtual void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet> updateVar); 
  
  /**
   * Set the update variable set
   */
  virtual void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar);
  
  /**
   * Compute the diffusive term flux in the current cell
   */
  virtual void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
				    std::vector<RealVector>& result,
				    const bool updateCoeff);
  
  /**
   * Compute the cell gradient and average state and put them
   * into @see DistributionData
   */
  virtual void computeCellGradientsAndAverageState
  (Framework::GeometricEntity *const geo, const RealVector& pdata);
  
  virtual void computePicardDiffJacob(Framework::GeometricEntity *const cell,
			      std::vector<RealMatrix*>& diffjacob);
  
protected: // data

  // acquaintance of the update convective var set
  Common::SafePtr<Physics::Poisson::PoissonConvVarSet> _upVar;
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::Poisson::PoissonDiffVarSet> _diffVar;
    
  // array of cell states
  std::vector<RealVector*> _states;
  
  // array of values (rho, u, v, T)
  RealMatrix _values;
  
  // array of average values (rho, u, v, T)
  Framework::State* _avValues;
  
  // unit normal
  RealVector _normal;
    
  /// temporary residual pointer after transformation
  RealVector * m_phi;

  CFreal cellVolume ;
  
}; // end of class PoissonTerm
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PoissonTerm_hh
