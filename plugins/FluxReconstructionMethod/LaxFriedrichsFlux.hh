#ifndef COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh
#define COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"
#include <stdio.h>

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/MathTypes.hh"
#include "Framework/VarSetTransformerT.hh"
#include "FluxReconstructionMethod/FluxData.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Lax-Friedrichs/Rusanov flux
 *
 * @author Ray Vandenhoeck
 */
class LaxFriedrichsFlux : public RiemannFlux {

public:  // methods

  /// Defines the Config Options of this class
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures this Method.
  virtual void configure ( Config::ConfigArgs& args );
    
  #ifdef CF_HAVE_CUDA
  /// nested class defining local options
  template <typename P = NOTYPE>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() : epsilon(0.5) {}
    
    /// destructor
    HOST_DEVICE ~DeviceConfigOptions() {}
    
    /// initialize with another object of the same kind
    HOST_DEVICE void init(DeviceConfigOptions<P> *const in) 
    {
      epsilon = in->epsilon;
    }
    
    /// diffusion reduction coefficient for Rusanov flux
    CFreal epsilon;
  };
  
  /// nested class defining a functor
  template <DeviceType DT, typename VS, CFuint ORDER>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef LaxFriedrichsFlux BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
 
    /// Compute needed private variables that do not depend on the physical state
    HOST_DEVICE void prepareComputation(FluxData<VS,ORDER>* data, VS* model) {}
   
    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS,ORDER>* data, VS* model, const CFuint iFlxPnt, const CFuint flxIdx, const CFreal intCoeff, const bool isLEFT, CFreal &waveSpeedUpd);

    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS,ORDER>* data, VS* model, const CFuint iSol); 
    
  private:
    DeviceConfigOptions<NOTYPE>* m_dco;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmp;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmp2;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmpState;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmpState2;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC m_pdata;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC m_pdata2;
    typename MathTypes<CFreal, DT, VS::DIM*VS::DIM>::VEC m_tempUnitNormal;
    typename MathTypes<CFreal, DT, VS::DIM>::VEC m_tempFlxUnitNormal;
  };
  
  /// copy the local configuration options to the device
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
    dco->epsilon = m_epsilon;
  }  
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    dco->epsilon = m_epsilon;
  }   
#endif

  /// Constructor
  LaxFriedrichsFlux(const std::string& name);

  /// Destructor
  ~LaxFriedrichsFlux();
 
  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState,
				  Framework::State& rState,
				  const RealVector& normal);

  /// Compute the Riemann flux, from left and right states and a normal vector
  virtual RealVector& computeFlux(Framework::State& lState, RealVector& lExtraVars,
                                  Framework::State& rState, RealVector& rExtraVars,
                                  const RealVector& normal);


  /// Gets the Class name
  static std::string getClassName()
  {
    return "LaxFriedrichsFlux";
  }

  /// Set up private data and data
  virtual void setup();
  
  /// Unset up private data and data
  virtual void unsetup();

private: // data

  /// array storing the sum of the right and left flux
  RealVector  m_sumFlux;

  /// pre-allocated left solution state (avoids heap allocation per flux call)
  RealVector  m_lSolState;

  /// pre-allocated right solution state (avoids heap allocation per flux call)
  RealVector  m_rSolState;

  /// diffusion reduction coefficient (epsilon) for tuning the dissipative part of the Rusanov flux
  CFreal m_epsilon;

}; // class LaxFriedrichsFlux

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA
/// nested class defining the flux
template <DeviceType DT, typename VS, CFuint ORDER>
void LaxFriedrichsFlux::DeviceFunc<DT, VS, ORDER>::operator()(FluxData<VS,ORDER>* data, VS* model, const CFuint iFlxPnt, const CFuint flxIdx, const CFreal intCoeff, const bool isLEFT, CFreal &waveSpeedUpd)
{
    typename VS::UPDATE_VS* updateVS = model->getUpdateVS();

    // right physical data, flux
    updateVS->computePhysicalData(data->getLstate(flxIdx), &m_pdata[0]);
    updateVS->computePhysicalData(data->getRstate(flxIdx), &m_pdata2[0]);

    typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getFlxScaledNormal(flxIdx));
    m_tempFlxUnitNormal = unitNormal;

    CFreal nJacob2 = 0.0;

    for (CFuint iDim = 0; iDim < VS::DIM; ++iDim)
    {
      nJacob2 += m_tempFlxUnitNormal[iDim]*m_tempFlxUnitNormal[iDim];
    }

    const CFreal nJacob = pow(nJacob2,0.5);

    for (CFuint iDim = 0; iDim < VS::DIM; ++iDim)
    {
      m_tempFlxUnitNormal[iDim] *= 1.0/nJacob;
    }
    
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getInterfaceFlux(flxIdx));
        
    updateVS->getFlux(&m_pdata[0], &m_tempFlxUnitNormal[0], &m_tmp[0]); 
    updateVS->getFlux(&m_pdata2[0], &m_tempFlxUnitNormal[0], &m_tmp2[0]); 

    

    const CFreal lMaxAbsEVal = updateVS->getMaxAbsEigenValue(&m_pdata[0], &m_tempFlxUnitNormal[0]);
    const CFreal rMaxAbsEVal = updateVS->getMaxAbsEigenValue(&m_pdata2[0], &m_tempFlxUnitNormal[0]);

    const CFreal absA = 0.5*(lMaxAbsEVal+rMaxAbsEVal);
    
    Framework::VarSetTransformerT<typename VS::UPDATE_VS, typename VS::SOLUTION_VS, NOTYPE>* up2Sol = model->getUpdateToSolution();
  
    // transform to solution variables
    up2Sol->transform(data->getLstate(flxIdx), &m_tmpState[0]);
    up2Sol->transform(data->getRstate(flxIdx), &m_tmpState2[0]);
  
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC stateL(&m_tmpState[0]);
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC stateR(&m_tmpState2[0]); 

    //CFreal* intCoeff = data->getFaceIntegrationCoef();
    //const CFreal intCoeff = data->getFaceIntegrationCoef(iFlxPnt);

    if (isLEFT)
    { 
      for (CFuint iEq = 0; iEq < VS::NBEQS; ++iEq)
      {
        // Rusanov flux: F = 0.5*(FL+FR) - epsilon*|A|*(UR-UL)
        flux[iEq] = (0.5*(m_tmp[iEq]+m_tmp2[iEq]) - m_dco->epsilon*absA*(stateR[iEq] - stateL[iEq]))*nJacob;
//if (cellID == 768) printf("flux %f, var %d, lState %f, rState %f, lflux %f, rflux %f, lEV %f, rEV %f, normalX %f, normalY %f, jacob %f\n",flux[iEq],iEq,stateL[iEq],stateR[iEq],m_tmp[iEq],m_tmp2[iEq],lMaxAbsEVal,rMaxAbsEVal,m_tempFlxUnitNormal[0],m_tempFlxUnitNormal[1],nJacob);
      }
    }
    else
    {
      for (CFuint iEq = 0; iEq < VS::NBEQS; ++iEq)
      {
        // Rusanov flux: F = 0.5*(FL+FR) - epsilon*|A|*(UL-UR)
        flux[iEq] = (0.5*(m_tmp[iEq]+m_tmp2[iEq]) - m_dco->epsilon*absA*(stateL[iEq] - stateR[iEq]))*nJacob;
//if (cellID == 768) printf("flux %f, var %d, lState %f, rState %f, lflux %f, rflux %f, lEV %f, rEV %f, normalX %f, normalY %f, jacob %f\n",flux[iEq],iEq,stateL[iEq],stateR[iEq],m_tmp[iEq],m_tmp2[iEq],lMaxAbsEVal,rMaxAbsEVal,m_tempFlxUnitNormal[0],m_tempFlxUnitNormal[1],nJacob);
      }
    }
//if (cellID == 768) printf("upBefore: %f\n", waveSpeedUpd);
    // compute the wave speed updates
    //*(data->getUpdateCoeff()) = *(data->getUpdateCoeff()) + nJacob * 1.0 * lMaxAbsEVal;
    waveSpeedUpd += nJacob * intCoeff * lMaxAbsEVal;
    
//printf("intCoeff: %f\n", *intCoeff);
//if (cellID == 11) printf("cellID: %d, upd: %f\n",cellID,*(data->getUpdateCoeff()));
//if (cellID == 768) printf("iFlx: %d, maxAbs: %f, intCoeff: %f, jacob: %f\n", flxIdx, lMaxAbsEVal, 1.0, nJacob);
//if (cellID == 768) printf("upAfter: %f\n", waveSpeedUpd);
    //data->addUpdateCoeff(waveSpeedUpd);
}
#endif

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA
/// nested class defining the flux
template <DeviceType DT, typename VS, CFuint ORDER>
void LaxFriedrichsFlux::DeviceFunc<DT, VS, ORDER>::operator()(FluxData<VS,ORDER>* data, VS* model, const CFuint iSol)
{
  typename VS::UPDATE_VS* updateVS = model->getUpdateVS();

  // right physical data, flux
  updateVS->computePhysicalData(data->getState(iSol), &m_pdata[0]);

  typename MathTypes<CFreal,DT,VS::DIM*VS::DIM>::SLICEVEC unitNormal(data->getScaledNormal(iSol));
  m_tempUnitNormal = unitNormal;
  for (CFuint iDim = 0; iDim < VS::DIM; ++iDim)
  {
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getFlux(iSol,iDim));
        
    updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[iDim*VS::DIM], &m_tmp[0]); 
      
    for (CFuint iEq = 0; iEq < VS::NBEQS; ++iEq)
    {
      flux[iEq] = m_tmp[iEq];
    }
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh
