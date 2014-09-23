#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadScalarRotDiffusionTerm_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadScalarRotDiffusionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace RotationAdv {
      class RotationDiffusionVarSet;
      class RotationAdv2DVarSet;
    }    
  }
  
    namespace FluctSplit {
            
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Advection diffusion
 * physical model
 *
 * @author Nadege Villedieu
 *
 */
class UnsteadScalarRotDiffusionTerm : public ComputeDiffusiveTerm {
public:

  /**
   * Constructor
   */
  UnsteadScalarRotDiffusionTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~UnsteadScalarRotDiffusionTerm();
  
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

  void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob){
    throw Common::NotImplementedException (FromHere(),"UnsteadyScalarRotDiffTerm::computePicardJacob()");
  }
protected: // data
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::RotationAdv::RotationDiffusionVarSet> _diffVar;
  
  // acquaintance of the update convective var set
  Common::SafePtr<Physics::RotationAdv::RotationAdv2DVarSet> _updateVar;
    
  // array of cell states
  std::vector<RealVector*> _states;
  
  // array of values (rho, u, v, T)
  std::vector<RealVector*> _values;
  
  // array of gradients
  std::vector<RealVector*> _gradients;
  
  // array of average values (rho, u, v, T)
  Framework::State* _avValues;
  
  // unit normal
  RealVector _normal;
  
  /// temporary residual pointer after transformation
  RealVector * m_phi;

    /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// temporary storage of the past states
  std::vector<Framework::State*> _pastStates;

}; // end of class UnsteadAdvectionDiffusionTerm
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadAdvectionDiffusionTerm_hh
