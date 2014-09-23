#ifndef COOLFluiD_Numerics_FluctSplit_NavierStokesTerm_hh
#define COOLFluiD_Numerics_FluctSplit_NavierStokesTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
      class EulerVarSet;
    }    
  }
  
  namespace FluctSplit {
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model
 *
 * @author Andrea Lani
 * @author Nadege Villedieu
 *
 */
template <typename BASE>
class NavierStokesTerm : public BASE {
public:

  /**
   * Constructor
   */
  NavierStokesTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesTerm();
  
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
  
  virtual void computeJacobNSFlux(Framework::State* state, RealVector& normal, RealMatrix& A, RealMatrix& B);
  
protected: // data
   
  // acquaintance of the update convective var set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _upVar;
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVar;
  
  // average radius in the cell
  CFreal _radius;
  
  // array of cell states
  std::vector<RealVector*> _states;
  
  // array of values (rho, u, v, T)
  RealMatrix _values;
  
  // array of average values (rho, u, v, T)
  Framework::State* _avValues;
  
  // unit normal
  RealVector _normal;
  RealVector normaljState;
  
  /// temporary residual pointer after transformation
  RealVector * m_phi;

  CFreal cellVolume ;

  RealMatrix Ax;
  RealMatrix Ay;
  
}; // end of class NavierStokesTerm
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NavierStokesTerm_hh
