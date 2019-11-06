#ifndef COOLFluiD_FluxReconstructionMethod_FluxData_hh
#define COOLFluiD_FluxReconstructionMethod_FluxData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {
            
//////////////////////////////////////////////////////////////////////////////

template <typename PHYS, CFuint ORDER>      
class FluxData {
public:
  /// initialize the flux data
  HOST_DEVICE void initialize() 
  {
    const CFuint nbSolPnts = (ORDER+1)*(ORDER+1);
    const CFuint nbFaceFlxPnts = ORDER+1;
    const CFuint nbFlxPnts = 2*PHYS::DIM*(ORDER+1);
    for (CFuint j = 0; j < nbSolPnts; j++) {for (CFuint i = 0; i < PHYS::NBEQS; ++i) {m_states[j][i] = 0.;}}
    for (CFuint j = 0; j < nbFlxPnts; j++) {for (CFuint i = 0; i < PHYS::NBEQS; ++i) {m_rstates[j][i] = 0.;}}
    for (CFuint j = 0; j < nbFlxPnts; j++) {for (CFuint i = 0; i < PHYS::NBEQS; ++i) {m_lstates[j][i] = 0.;}}
    for (CFuint i = 0; i < nbSolPnts; ++i) {for (CFuint j = 0; j < PHYS::DIM; ++j) {for (CFuint k = 0; k < PHYS::NBEQS; ++k) {m_flux[i][j][k] = 0.;}}}
    for (CFuint i = 0; i < nbFaceFlxPnts; ++i) {for (CFuint j = 0; j < PHYS::NBEQS; ++j) {m_fluxFlxPnt[i][j] = 0.;}}
    for (CFuint i = 0; i < nbSolPnts; ++i) {for (CFuint j = 0; j < PHYS::DIM*PHYS::DIM; ++j) {m_unitNormal[i][j] = 0.;}}
    for (CFuint i = 0; i < nbFlxPnts; ++i) {for (CFuint j = 0; j < PHYS::DIM; ++j) {m_unitFlxNormal[i][j] = 0.;}}
    for (CFuint i = 0; i < nbFaceFlxPnts; ++i) {m_faceIntegrationCoefs[i] = 0.;}
    m_updateCoeff = 0.; m_faceArea = 0.;
    m_stateID[0] = 0; m_stateID[1] = 0; m_stateID[2] = 0; m_stateID[3] = 0;
    m_nbSolPnts = nbSolPnts;
    m_isBFace = false; m_isOutward = false; m_isPerturb = false; 
  }
    
  /// get the left or right state
  HOST_DEVICE CFreal* getState(const CFuint iState) {return &m_states[iState][0];}
    
  /// get the right reconstructed state
  HOST_DEVICE CFreal* getRstate(const CFuint iState) {return &m_rstates[iState][0];}
  
  /// get the left reconstructed state
  HOST_DEVICE CFreal* getLstate(const CFuint iState) {return &m_lstates[iState][0];}
  
  /// get the flux array
  HOST_DEVICE CFreal* getFlux(const CFuint iSol, const CFuint iDim) {return &m_flux[iSol][iDim][0];}
  
  /// get the face normal scaled with jacobian
  HOST_DEVICE CFreal* getScaledNormal(const CFuint iSol) {return &m_unitNormal[iSol][0];}
  
  /// get the face normal scaled with jacobian
  HOST_DEVICE CFreal* getFlxScaledNormal(const CFuint iFlx) {return &m_unitFlxNormal[iFlx][0];}
  
  /// get the flux on a flux point
  HOST_DEVICE CFreal* getInterfaceFlux(const CFuint iFlx) {return &m_fluxFlxPnt[iFlx][0];}  
  
  /// get the face area
  HOST_DEVICE CFreal getFaceArea() const {return m_faceArea;}
  
  /// set the face area
  HOST_DEVICE void setFaceArea(const CFreal area) {m_faceArea = area;}
  
  /// get the boundary flag
  HOST_DEVICE bool isBFace() const {return m_isBFace;}
  
  /// set the boundary flag
  HOST_DEVICE void setIsBFace(const bool flag) {m_isBFace = flag;}
  
  /// get the flag telling if the normal is outward
  HOST_DEVICE bool isOutward() const {return m_isOutward;}
  
  /// set the flag telling if the normal is outward
  HOST_DEVICE void setIsOutward(const bool flag) {m_isOutward = flag;}
  
  /// get the flag telling if jacobian perturbation is applied
  HOST_DEVICE bool isPerturb() const {return m_isPerturb;}
  
  /// set the flag telling if jacobian perturbation is applied
  HOST_DEVICE void setIsPerturb(bool flag) {m_isPerturb = flag;}
  
  /// get the state IDs
  HOST_DEVICE CFuint getStateID(const CFuint iState) const {return m_stateID[iState];}
  
  /// get number of solution points
  HOST_DEVICE CFuint getNbSolPnts() const {return m_nbSolPnts;}
  
  /// set number of solution points
  HOST_DEVICE void setNbSolPnts(const CFuint nbSolPnts) {m_nbSolPnts = nbSolPnts;}
  
  /// set the state IDs
  HOST_DEVICE void setStateID(const CFuint iState, const CFuint sID) {m_stateID[iState] = sID;}
  
  /// get the update coefficient
  HOST_DEVICE CFreal* getUpdateCoeff() {return &m_updateCoeff;}
  
  /// set the update coefficient
  HOST_DEVICE void addUpdateCoeff(const CFreal coeff) {m_updateCoeff += coeff;}

  /// get the update coefficient
  HOST_DEVICE void resetUpdateCoeff() {m_updateCoeff = 0.;}

  /// get face integration coefficient
  HOST_DEVICE CFreal* getFaceIntegrationCoef() {return &m_faceIntegrationCoefs[0];}

  /// set face integration coefficient
  HOST_DEVICE void setFaceIntegrationCoef(const CFuint iFlx, const CFreal coeff) {m_faceIntegrationCoefs[iFlx] = coeff;}
    
private:
  
  /// inner states
  CFreal m_states[(ORDER+1)*(ORDER+1)][PHYS::NBEQS];
   
  /// right reconstructed states
  CFreal m_rstates[2*PHYS::DIM*(ORDER+1)][PHYS::NBEQS];
  
  /// left reconstructed states
  CFreal m_lstates[2*PHYS::DIM*(ORDER+1)][PHYS::NBEQS];
  
  /// flux in solution points
  CFreal m_flux[(ORDER+1)*(ORDER+1)][PHYS::DIM][PHYS::NBEQS];
  
  /// flux in flux points
  CFreal m_fluxFlxPnt[ORDER+1][PHYS::NBEQS];
  
  /// unit normal
  CFreal m_unitNormal[(ORDER+1)*(ORDER+1)][PHYS::DIM*PHYS::DIM];
  
  /// unit normal
  CFreal m_unitFlxNormal[2*PHYS::DIM*(ORDER+1)][PHYS::DIM];
  
  /// left and right update coefficient
  CFreal m_updateCoeff;
  
  /// face length/area
  CFreal m_faceArea;
  
  /// state IDs
  CFuint m_stateID[(ORDER+1)*(ORDER+1)];
  
  /// number of solution points
  CFuint m_nbSolPnts;
  
  /// flag telling if this is a boundary face
  bool m_isBFace;
  
  /// flag telling if the normal is outward
  bool m_isOutward;
  
  /// flag telling if jacobian perturbation has to be applied
  bool m_isPerturb;

  /// face integration coefficients
  CFreal m_faceIntegrationCoefs[ORDER+1];
  
};
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_FluxData_hh



