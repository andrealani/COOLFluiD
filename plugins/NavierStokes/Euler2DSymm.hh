#ifndef COOLFluiD_Physics_NavierStokes_Euler2DSymm_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DSymm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Euler2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Euler physical model 3D for conservative
   * variables
   *
   * @author Andrea Lani
   */

class Euler2DSymm : public Euler2DVarSet {
public:

  /**
   * Constructor
   * @see Euler3D
   */
  Euler2DSymm(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  ~Euler2DSymm();

  /**
   * Gets the block separator for this variable set
   */
  CFuint getBlockSeparator() const;
  
  /// Set the PhysicalData corresponding to the given State
  void computePhysicalData(const Framework::State& state, RealVector& data)
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DSymm::computePhysicalData() not implemented");
  }
  
  /// Set the State correspoding to the PhysicalData
  virtual void computeStateFromPhysicalData(const RealVector& pdata, Framework::State& state) 
  {
    throw Common::NotImplementedException(FromHere(), "Euler2DSymm::computeStateFromPhysicalData() not implemented");
  }
  
  /**
   * Set the jacobians
   */
  virtual void computeJacobians();

  /**
   * Set the scalar part of the jacobian
   */
  void computeScalarJacobian(const RealVector& normal,
                      RealVector& jacob);


  /**
   * Split the jacobian
   */
  void splitJacobian(RealMatrix& jacobPlus,
                  RealMatrix& jacobMin,
                  RealVector& eValues,
                  const RealVector& normal);
 
  
  
  /**
   * EXPERIMENT: Split the jacobian
   */
  template <typename FUN>
  void splitJacobian(RealMatrix& jacobPlus,
		     RealMatrix& jacobMin,
		     RealVector& eValues,
		     const RealVector& normal,
		     FUN& f)
  {
    const RealVector& linearData = getModel()->getPhysicalData();
    
    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];
    const CFreal avU   = linearData[EulerTerm::VX];
    const CFreal avV   = linearData[EulerTerm::VY];
    const CFreal avA   = linearData[EulerTerm::A];
    const CFreal speed = sqrt(avU*avU + avV*avV);
    const CFreal cosD = avU/speed;
    const CFreal sinD = avV/speed;
    const CFreal nCsi = nx*cosD + ny*sinD;
    const CFreal nEta = -nx*sinD + ny*cosD;
    const CFreal nEtaCsi = nEta*nCsi;
    const CFreal nCsi2 = nCsi*nCsi;
    const CFreal nEta2 = nEta*nEta;
    const CFreal Vn  = avU*nx + avV*ny;
    const CFuint sizeJacob = jacobPlus.nbRows();
    
    cf_assert(eValues.size() == 4);
    eValues[0] = eValues[1] = Vn;
    eValues[2] = Vn + avA;
    eValues[3] = Vn - avA;
    
    static RealVector eValuesCorr(4);
    // correct the eigenvalues
    f.correctEigenvalues(eValues, eValuesCorr);
    
    for (CFuint iEq = 0; iEq < 4; ++iEq) {
      _eValuesP[iEq] = std::max(0.,eValuesCorr[iEq]);
      _eValuesM[iEq] = std::min(0.,eValuesCorr[iEq]);
    }
    
    CFreal e2 = _eValuesP[1];
    CFreal e3Min4 = 0.5*(_eValuesP[2] - _eValuesP[3]);
    CFreal e3Plus4 = 0.5*(_eValuesP[2] + _eValuesP[3]);
    CFreal e3Min4nEta = e3Min4*nEta;
    CFreal e3Min4nCsi = e3Min4*nCsi;
    CFreal e234 = -e2 + e3Plus4;
    
    if (sizeJacob == 3) {
      jacobPlus(0,0) = e3Plus4;
      jacobPlus(0,1) = e3Min4*nCsi;
      jacobPlus(0,2) = e3Min4*nEta;
      jacobPlus(1,0) = e3Min4*nCsi;
      jacobPlus(1,1) = e2*nEta2 + e3Plus4*nCsi2;
      jacobPlus(1,2) = nEtaCsi*e234;
      jacobPlus(2,0) = e3Min4*nEta;
      jacobPlus(2,1) = nEtaCsi*e234;
      jacobPlus(2,2) = e2*nCsi2 + e3Plus4*nEta2;
      
      //  jacobPlus[0] = e3Plus4;
      //  jacobPlus[1] = e3Min4*nCsi;
      //  jacobPlus[2] = e3Min4*nEta;
      
      //  jacobPlus[4] = e3Min4*nCsi;
      //  jacobPlus[5] = e2*nEta2 + e3Plus4*nCsi2;
      //  jacobPlus[6] = nEtaCsi*e234;
      
      //  jacobPlus[8] = e3Min4*nEta;
      //  jacobPlus[9] = nEtaCsi*e234;
      //  jacobPlus[10] = e2*nCsi2 + e3Plus4*nEta2;
      
      e2 = _eValuesM[1];
      e3Min4 = 0.5*(_eValuesM[2] - _eValuesM[3]);
      e3Plus4 = 0.5*(_eValuesM[2] + _eValuesM[3]);
      e3Min4nEta = e3Min4*nEta;
      e3Min4nCsi = e3Min4*nCsi;
      e234 = -e2 + e3Plus4;
      
      jacobMin(0,0) = e3Plus4;
      jacobMin(0,1) = e3Min4*nCsi;
      jacobMin(0,2) = e3Min4*nEta;
      jacobMin(1,0) = e3Min4*nCsi;
      jacobMin(1,1) = e2*nEta2 + e3Plus4*nCsi2;
      jacobMin(1,2) = nEtaCsi*e234;
      jacobMin(2,0) = e3Min4*nEta;
      jacobMin(2,1) = nEtaCsi*e234;
      jacobMin(2,2) = e2*nCsi2 + e3Plus4*nEta2;
    }
    else {
      cf_assert(sizeJacob == 4);
      
      jacobPlus(0,0) = e3Plus4;
      jacobPlus(0,1) = e3Min4*nCsi;
      jacobPlus(0,2) = e3Min4*nEta;
      jacobPlus(1,0) = e3Min4*nCsi;
      jacobPlus(1,1) = e2*nEta2 + e3Plus4*nCsi2;
      jacobPlus(1,2) = nEtaCsi*e234;
      jacobPlus(2,0) = e3Min4*nEta;
      jacobPlus(2,1) = nEtaCsi*e234;
      jacobPlus(2,2) = e2*nCsi2 + e3Plus4*nEta2;
      jacobPlus(3,3) = e2; //e2 = e1
      
      e2 = _eValuesM[1];
      e3Min4 = 0.5*(_eValuesM[2] - _eValuesM[3]);
      e3Plus4 = 0.5*(_eValuesM[2] + _eValuesM[3]);
      e3Min4nEta = e3Min4*nEta;
      e3Min4nCsi = e3Min4*nCsi;
      e234 = -e2 + e3Plus4;
      
      jacobMin(0,0) = e3Plus4;
      jacobMin(0,1) = e3Min4*nCsi;
      jacobMin(0,2) = e3Min4*nEta;
      jacobMin(1,0) = e3Min4*nCsi;
      jacobMin(1,1) = e2*nEta2 + e3Plus4*nCsi2;
      jacobMin(1,2) = nEtaCsi*e234;
      jacobMin(2,0) = e3Min4*nEta;
      jacobMin(2,1) = nEtaCsi*e234;
      jacobMin(2,2) = e2*nCsi2 + e3Plus4*nEta2;
      jacobMin(3,3) = e2; //e2 = e1
    } 
  }
  
  
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  void computeEigenValuesVectors(RealMatrix& rightEv,
                             RealMatrix& leftEv,
                             RealVector& eValues,
                             const RealVector& normal);

  /**
   * Get the speed
   */
  CFreal getSpeed(const Framework::State& state) const
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DSymm::getSpeed() not implemented");
  }

  /**
   * Give dimensional values to the adimensional state variables
   */
  void setDimensionalValues(const Framework::State& state,
                            RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DSymm::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  void setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "Euler2DSymm::setDimensionalValues() not implemented");
  }
  
  /// Compute the perturbed physical data
  virtual void computePerturbedPhysicalData(const Framework::State& state,
					    const RealVector& pdataBkp,
					    RealVector& pdata,
					    CFuint iVar)
  {
    throw Common::NotImplementedException (FromHere(),"Euler2DSymm::computePerturbedPhysicalData()");
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  { 
    velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2;
  }
  
}; // end of class Euler2DSymm

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DSymm_hh
