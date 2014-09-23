#ifndef COOLFluiD_Numerics_FluctSplit_UnsteadyNavierStokesTermquad_hh
#define COOLFluiD_Numerics_FluctSplit_UnsteadyNavierStokesTermquad_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "ComputeDiffusiveTerm.hh"

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
 *
 */
class UnsteadyNavierStokesTermquad : public ComputeDiffusiveTerm {
public:

  /**
   * Constructor
   */
  UnsteadyNavierStokesTermquad(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~UnsteadyNavierStokesTermquad();
  
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
  virtual void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar);
   
  /**
   * Compute the diffusive term flux in the current cell
   */
  virtual void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
			    std::vector<RealVector>& result,
			    const bool updateCoeff);
  /**
    * Set the update variable set
    */
  virtual void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet>); 
  


/**
   * Returns the DataSocket's that this numerical strategy needs as sinks 
   * @return a vector of SafePtr with the DataSockets 
   */ 
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets() 
  { 
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeDiffusiveTerm::needsSockets(); 
    result.push_back(&socket_pastStates); 
    return result; 
  } 

  
  virtual void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob){
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

 
  void computeDiffRes(bool isPast, std::vector<RealVector>& result, bool updateCoeffFlag);
  void computeCellGradientsAndAverageStatequad ( const RealVector& pdata);
protected: // data
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVar;
  
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

  CFuint nbCellStates;
  CFuint nbEqs;
  CFreal dtCoeff;
  CFreal volume;
  CFreal ovNbStatesInCell;
  CFuint cellID;
  
}; // end of class UnsteadyNavierStokesTerm
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_UnsteadyNavierStokesTermquad_hh
