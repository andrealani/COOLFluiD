#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyWALESTerm_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyWALESTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    namespace NavierStokes {
      class EulerVarSet;
    }  
     namespace LESvki {
       class WALES2DVarSet;
   
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
class UnsteadyWALESTerm : public ComputeDiffusiveTerm {
public:

  /**
   * Constructor
   */
  UnsteadyWALESTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~UnsteadyWALESTerm();
  
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
  
  /**
   * Compute the cell gradient and average state and put them
   * into @see DistributionData
   */
  virtual void computeCellGradientsAndAverageState
  (Framework::GeometricEntity *const geo, const RealVector& pdata);
  


/**
   * Returns the DataSocket's that this numerical strategy needs as sinks 
   * @return a vector of SafePtr with the DataSockets 
   */ 
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets() 
  { 
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeDiffusiveTerm::needsSockets(); 
    result.push_back(&socket_pastStates); 
    return result; 
  } 

  
  void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob){
throw Common::NotImplementedException (FromHere(),"UnsteadyNavierStokes::computePicardJacob()");
}
protected:

  /// Compute the average state given a set of states
  void computeAverageState(CFuint nbStatesInCell, CFreal ovNbStatesInCell, 
			   CFuint nbEqs, const std::vector<RealVector*>& states,
			   RealVector& avState) 
  {
    for (CFuint j = 0; j < nbEqs; ++j) {
      (*_avState)[j] = (*_states[0])[j];
      for (CFuint i = 1; i < nbStatesInCell; ++i) {
	(*_avState)[j] += (*_states[i])[j];
      }
    }
    (*_avState) *= ovNbStatesInCell;
  }
  
protected: // data
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::LESvki::WALES2DVarSet> _diffVar;
  
  // acquaintance of the update convective var set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> _updateVar;
    
  // average radius in the cell
  CFreal _radius;
    
  // array of cell states
  std::vector<RealVector*> _states;
  
  // array of values (rho, u, v, T)
  RealMatrix _values;
  
  // array of average values (rho, u, v, T)
  Framework::State* _avValues;
  
  /// average state
  Framework::State* _avState;
  
  // unit normal
  RealVector _normal;
    
  /// gradient rho
  RealVector _gradRho;
  
  /// gradient rhoU
  RealVector _gradRhoU;
  
  /// gradient rhoV
  RealVector _gradRhoV;
   
  // physical data
  RealVector _edata;
   
  /// temporary residual pointer after transformation
  RealVector * m_phi;
  
  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;
  
}; // end of class UnsteadyNavierStokesTerm
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyNavierStokesTerm_hh
