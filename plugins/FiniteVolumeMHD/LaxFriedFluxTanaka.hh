#ifndef COOLFluiD_Numerics_FiniteVolume_LaxFriedFluxTanaka_hh
#define COOLFluiD_Numerics_FiniteVolume_LaxFriedFluxTanaka_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/LaxFriedFlux.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents the Lax-Friedrichs flux according to Tanaka's
 * approach (magnetic field splitting in case of a background potential field) 
 * corresponding to the MHD physical model (corresponding to the update variable 
 * set UPDATEVAR) according to the Powell99 paper
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
template <typename UPDATEVAR>
class LaxFriedFluxTanaka : public LaxFriedFlux {
public:
  
#ifdef CF_HAVE_CUDA
  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef LaxFriedFluxTanaka<UPDATEVAR> BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(LaxFriedFlux::DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}

   
    /// Compute needed private variables that do not depend on the physical state
    HOST_DEVICE void prepareComputation(FluxData<VS>* data, VS* model) {}

    
    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS>* data, VS* model)
    {
      typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getUnitNormal());
      const CFreal coeff = (data->isOutward()) ? 1. : -1.;
      m_tempUnitNormal = coeff*unitNormal;
      
      typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getResidual());
      
      typename VS::UPDATE_VS* updateVS = model->getUpdateVS();
      // right physical data, flux and eigenvalues
      updateVS->computePhysicalData(data->getRstate(1), data->getRnode(1), &m_pdata[0]);
      updateVS->computeTanakaFluxPowell99Formulation(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
      flux = 0.5*m_tmp;
      
      updateVS->computeEigenValues(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
      CFreal aR = 0.0;
      for (CFuint i = 0; i < VS::NBEQS; ++i) {
    	aR = max(aR, abs(m_tmp[i]));
      }
      
      // left physical data, flux and eigenvalues
      updateVS->computePhysicalData(data->getRstate(0), data->getRnode(0), &m_pdata[0]);
      updateVS->computeTanakaFluxPowell99Formulation(&m_pdata[0], &m_tempUnitNormal[0], &m_tmp[0]);
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
    
  private:
    LaxFriedFlux::DeviceConfigOptions<NOTYPE>* m_dco;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmp;
    typename MathTypes<CFreal, DT, VS::NBEQS>::VEC m_tmp2;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC m_pdata;
    typename MathTypes<CFreal, DT, VS::DIM>::VEC m_tempUnitNormal;
  };
#endif
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  LaxFriedFluxTanaka(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LaxFriedFluxTanaka();

  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
private:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> _updateVarSet;
  
}; // end of class LaxFriedFluxTanaka

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#include "LaxFriedFluxTanaka.ci"

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LaxFriedFluxTanaka_hh
