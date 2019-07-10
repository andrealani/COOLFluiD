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
    
  #ifdef CF_HAVE_CUDA
  /// nested class defining local options
  template <typename P = NOTYPE>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() {}
    
    /// destructor
    HOST_DEVICE ~DeviceConfigOptions() {}
    
    /// initialize with another object of the same kind
    HOST_DEVICE void init(DeviceConfigOptions<P> *const in) 
    {
    }
  };
  
  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef LaxFriedrichsFlux BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
 
    /// Compute needed private variables that do not depend on the physical state
    HOST_DEVICE void prepareComputation(FluxData<VS>* data, VS* model) {}
   
    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS>* data, VS* model, bool isInterface); 
    
  private:
    DeviceConfigOptions<NOTYPE>* m_dco;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmp;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmp2;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC m_pdata;
    typename MathTypes<CFreal, DT, VS::DIM>::VEC m_tempUnitNormal;
  };
  
  /// copy the local configuration options to the device
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
  }  
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
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


}; // class LaxFriedrichsFlux

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA
/// nested class defining the flux
template <DeviceType DT, typename VS>
void LaxFriedrichsFlux::DeviceFunc<DT, VS>::operator()(FluxData<VS>* data, VS* model, bool isInterface) 
{
  if (isInterface)
  {
    typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getScaledNormal());
    const CFreal coeff = (data->isOutward()) ? 1. : -1.;
    m_tempUnitNormal = coeff*unitNormal;
  
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getFlux());
  
    typename VS::UPDATE_VS* updateVS = model->getUpdateVS();
    // right physical data, flux and eigenvalues
    updateVS->computePhysicalData(data->getRstate(1), &m_pdata[0]);
    updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
    flux = 0.5*m_tmp;
  
    const CFreal rMaxAbsEVal = updateVS->getMaxAbsEigenValue(&m_pdata[0], &m_tempUnitNormal[0]);
  
    // left physical data, flux and eigenvalues
    updateVS->computePhysicalData(data->getRstate(0), &m_pdata[0]);
    updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
    flux += 0.5*m_tmp;
  
    const CFreal lMaxAbsEVal = updateVS->getMaxAbsEigenValue(&m_pdata[0], &m_tempUnitNormal[0]);

    const CFreal absA = 0.5*(lMaxAbsEVal+rMaxAbsEVal);
  
    Framework::VarSetTransformerT<typename VS::UPDATE_VS, typename VS::SOLUTION_VS, NOTYPE>* up2Sol = model->getUpdateToSolution();
  
    // transform to solution variables
    up2Sol->transform(data->getRstate(LEFT), &m_tmp[0]);
    up2Sol->transform(data->getRstate(RIGHT), &m_tmp2[0]);
  
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC stateL(&m_tmp[0]);
    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC stateR(&m_tmp2[0]);
  
    flux -= 0.5*absA*(stateR - stateL);
  }
  else
  {
    typename VS::UPDATE_VS* updateVS = model->getUpdateVS();

    typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getFlux());

    // right physical data, flux
    updateVS->computePhysicalData(data->getState(0), &m_pdata[0]);
    if (data->getStateID(0) == 0) printf("state: %d \n",m_pdata[0]);
    typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getScaledNormal());
    m_tempUnitNormal = unitNormal;
    updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]); 
    
    flux = m_tmp;
  }
}
#endif

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_LaxFriedrichsFlux_hh
