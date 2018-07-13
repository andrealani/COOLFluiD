#ifndef COOLFluiD_Numerics_FiniteVolume_LaxFriedFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_LaxFriedFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"

#include "FiniteVolume/FVMCC_FluxSplitter.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/MathTypes.hh"
#include "Framework/VarSetTransformerT.hh"
#include "FiniteVolume/FluxData.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Lax-Friedrichs flux
 *
 * @author Andrea Lani
 * @author Sarp Yalim
 *
 */
class LaxFriedFlux : public FVMCC_FluxSplitter {
public:
  
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
      currentDiffRedCoeff = in->currentDiffRedCoeff;
    }
    
    /// diffusion reduction coefficient
    CFreal currentDiffRedCoeff;
  };
  
  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef LaxFriedFlux BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
 
    /// Compute needed private variables that do not depend on the physical state
    HOST_DEVICE void prepareComputation(FluxData<VS>* data, VS* model) {}
   
    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS>* data, VS* model); 
    
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
    // consider to copy to constant memory
    CFreal currentDiffRedCoeff = getReductionCoeff(); 
    CudaEnv::copyHost2Dev(&dco->currentDiffRedCoeff, &currentDiffRedCoeff, 1);
  }  
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    // consider to copy to constant memory
    dco->currentDiffRedCoeff = getReductionCoeff();
  }   
#endif
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  LaxFriedFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LaxFriedFlux();
  
  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Configuration
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
protected:
  
  /**
   * Compute the artificial diffusion reduction coefficient
   */
  virtual CFreal getReductionCoeff();
  
protected:
  
  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
  /// array storing the sum of the right and left flux
  RealVector   _sumFlux;

  /// array storing the temporary right eigenvalues
  RealVector    _rightEv;

  /// array storing the temporary left eigenvalues
  RealVector    _leftEv;

  /// temporary unit normal
  RealVector    _tempUnitNormal;
  
  /// diffusion reduction coefficient defined interactively
  CFreal _currentDiffRedCoeff;
  
}; // end of class LaxFriedFlux

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA
/// nested class defining the flux
template <DeviceType DT, typename VS>
void LaxFriedFlux::DeviceFunc<DT, VS>::operator()(FluxData<VS>* data, VS* model) 
{
  typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getUnitNormal());
  const CFreal coeff = (data->isOutward()) ? 1. : -1.;
  m_tempUnitNormal = coeff*unitNormal;
  
  typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getResidual());
  
  typename VS::UPDATE_VS* updateVS = model->getUpdateVS();
  // right physical data, flux and eigenvalues
  updateVS->computePhysicalData(data->getRstate(1), data->getRnode(1), &m_pdata[0]);
  updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
  flux = 0.5*m_tmp;
  
  updateVS->computeEigenValues(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
  CFreal aR = 0.0;
  for (CFuint i = 0; i < VS::NBEQS; ++i) {
    aR = max(aR, abs(m_tmp[i]));
  }
  
  // left physical data, flux and eigenvalues
  updateVS->computePhysicalData(data->getRstate(0), data->getRnode(0), &m_pdata[0]);
  updateVS->getFlux(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
  flux += 0.5*m_tmp;
  
  updateVS->computeEigenValues(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
    
  // compute update coefficient
  if (!data->isPerturb()) {    
    const CFreal k = max(m_tmp.max(), 0.)*data->getFaceArea();
    data->setUpdateCoeff(k);
  }
  
  CFreal aL = 0.0;
  for (CFuint i = 0; i < VS::NBEQS; ++i) {
    aL = max(aL, abs(m_tmp[i]));
  }
  
  const CFreal a = fmax(aR,aL);
  const CFreal aDiff = a*m_dco->currentDiffRedCoeff;
  
  Framework::VarSetTransformerT<typename VS::UPDATE_VS, typename VS::SOLUTION_VS, NOTYPE>* up2Sol = 
    model->getUpdateToSolution();
  
  // transform to solution variables
  up2Sol->transform(data->getRstate(LEFT), &m_tmp[0]);
  up2Sol->transform(data->getRstate(RIGHT), &m_tmp2[0]);
  
  typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC stateL(&m_tmp[0]);
  typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC stateR(&m_tmp2[0]);
  
  flux -= (0.5*aDiff)*(stateR - stateL);
  
  // NOTE THE AREA HERE !!!!!!!!!!!!!!!!
  flux *= data->getFaceArea();
}
#endif
  
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LaxFriedFlux_hh
