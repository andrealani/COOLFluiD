#ifndef COOLFluiD_Numerics_FiniteVolume_StegerWarmingMaxwellProjection2D_hh
#define COOLFluiD_Numerics_FiniteVolume_StegerWarmingMaxwellProjection2D_hh

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
 * This class computes the Steger Warming flux
 *
 * @author Andrea Lani
 * @author Alejandro Alvarez Laguna
 * @author Isaac Alonso
 *
 */
template <class UPDATEVAR>    
class StegerWarmingMaxwellProjection2D : public FVMCC_FluxSplitter {
public:
  
//New code
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
    
    /// This scheme does not have configurable parameters
  };

  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef StegerWarmingMaxwellProjection2D BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
    
    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS>* data, VS* model);


  private:
    DeviceConfigOptions<NOTYPE>* m_dco;

    typename MathTypes<CFreal, DT, VS::DIM>::VEC d_tempUnitNormal;

    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_pdata;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_lData;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_rData;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::MAT d_Aplus;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::MAT d_Aminus;


  };
  
  /// copy the local configuration
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
  }  
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
  }  


#endif



  /**
   * Constructor
   */
  StegerWarmingMaxwellProjection2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~StegerWarmingMaxwellProjection2D();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FVMCC_FluxSplitter::configure(args);
  }
  
  /**
   * Set up private data
   */
  virtual void setup();
  
  /**
   * Compute the flux : implementation
   */
  virtual void compute(RealVector& result);
  
   /**
   * Compute the Aplus Matrix
   */
  virtual void computeMatrixAplus();   
  
   /**
   * Compute the Aminus Matrix
   */
  virtual void computeMatrixAminus();  
    
  /**
   * Compute the left flux jacobian
   */
  virtual void computeLeftJacobian();
  
  /**
   * Compute the right flux jacobian
   */
  virtual void computeRightJacobian();
  
protected:
  
  /**
   * Compute the update coefficient in the standard way
   */
  void computeUpdateCoeff();

protected:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> m_updateVarSet;
  
  /// left data  
  RealVector* m_lData;
  
  /// right data
  RealVector* m_rData;
  
  /// temporary unit normal
  RealVector m_tempUnitNormal;
  
  /// array of physical data 
  RealVector m_pdata;
  
  /// matrix of right eigenvectors
  RealMatrix   _rightEv;

  /// matrix of left eigenvectors
  RealMatrix   _leftEv;

  /// vector of eigenvalues
  RealVector   _eValues;
  
  /// vector of eigenvalues
  RealVector   _absEvalues;

  /// abs of the jacobian matrix
  RealMatrix   _absJacob;
  
  /// right jacobian matrix
  RealMatrix   _jRight;
  
  /// left jacobian matrix
  RealMatrix   _jLeft;
  
  /// jacobian matrix
  RealMatrix   _jacob;
  
  /// jacobian matrix
  RealMatrix   _jacobLeftTransf;
  
  /// jacobian matrix
  RealMatrix   _jacobRightTransf;
  
  /// dummy jacobian matrix
  RealMatrix   _jacobDummy;
  
  /// vector with the electromagnetic field variables (LEFT)
  RealVector _EMField_l;

  /// vector with the electromagnetic field variables (RIGHT)
  RealVector _EMField_r;
  
  /// A plus Matrix
  RealMatrix   _Aplus;
  
  /// A minus Matrix
  RealMatrix   _Aminus;   
   
  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;
  
}; // end of class StegerWarmingMaxwellProjection2D

//////////////////////////////////////////////////////////////////////////////





#ifdef CF_HAVE_CUDA

template <class UPDATEVAR>
template <DeviceType DT, typename VS>
void StegerWarmingMaxwellProjection2D<UPDATEVAR>::DeviceFunc<DT, VS>::operator()(FluxData<VS>* data, VS* model) 
{

      typename VS::UPDATE_VS* updateVS = model->getUpdateVS(); 
      typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getResidual());

      // right physical data
      updateVS->computePhysicalData(data->getRstate(1), &d_rData[0]); 


      // left physical data
      updateVS->computePhysicalData(data->getRstate(0), &d_lData[0]);
           
     typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getUnitNormal());
     const CFreal coeff = (data->isOutward()) ? 1. : -1.;    
     d_tempUnitNormal = coeff*unitNormal;                              
     
 
      const CFreal gamma = updateVS->getDivBCleaningConst();	
      const CFreal chi = updateVS->getDivECleaningConst();
      const CFreal c_e = updateVS->getLightSpeed();
   //   printf("gamma = %f \t chi = %f \t c_e = %f \n", gamma, chi, c_e);
      
      
  d_Aplus(0,0) = (d_tempUnitNormal[1]*d_tempUnitNormal[1] + gamma*d_tempUnitNormal[0]*d_tempUnitNormal[0])*c_e;
  d_Aplus(0,1) = (gamma - 1)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aplus(0,2) = 0;
  d_Aplus(0,3) = 0;
  d_Aplus(0,4) = 0;
  d_Aplus(0,5) = d_tempUnitNormal[1];
  d_Aplus(0,6) = gamma*gamma*d_tempUnitNormal[0]; 
  d_Aplus(0,7) = 0; 

  d_Aplus(1,0) = (gamma - 1)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aplus(1,1) = (d_tempUnitNormal[0]*d_tempUnitNormal[0] + gamma*d_tempUnitNormal[1]*d_tempUnitNormal[1])*c_e;
  d_Aplus(1,2) = 0;
  d_Aplus(1,3) = 0;
  d_Aplus(1,4) = 0;
  d_Aplus(1,5) = -d_tempUnitNormal[0];
  d_Aplus(1,6) = gamma*gamma*d_tempUnitNormal[1]; 
  d_Aplus(1,7) = 0; 

  d_Aplus(2,0) = 0;
  d_Aplus(2,1) = 0;
  d_Aplus(2,2) = c_e;
  d_Aplus(2,3) = -d_tempUnitNormal[1];
  d_Aplus(2,4) = d_tempUnitNormal[0];
  d_Aplus(2,5) = 0;
  d_Aplus(2,6) = 0; 
  d_Aplus(2,7) = 0; 

  d_Aplus(3,0) = 0;
  d_Aplus(3,1) = 0;
  d_Aplus(3,2) = -d_tempUnitNormal[1]*c_e*c_e;
  d_Aplus(3,3) = (d_tempUnitNormal[1]*d_tempUnitNormal[1] + chi*d_tempUnitNormal[0]*d_tempUnitNormal[0])*c_e;
  d_Aplus(3,4) = (chi - 1)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aplus(3,5) = 0;
  d_Aplus(3,6) = 0;  
  d_Aplus(3,7) = chi*chi*d_tempUnitNormal[0]*c_e*c_e;  
  
  d_Aplus(4,0) = 0;
  d_Aplus(4,1) = 0;
  d_Aplus(4,2) = d_tempUnitNormal[0]*c_e*c_e;
  d_Aplus(4,3) = (chi -1)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aplus(4,4) = (d_tempUnitNormal[0]*d_tempUnitNormal[0] + chi*d_tempUnitNormal[1]*d_tempUnitNormal[1])*c_e;
  d_Aplus(4,5) = 0;
  d_Aplus(4,6) = 0; 
  d_Aplus(4,7) = chi*chi*d_tempUnitNormal[1]*c_e*c_e;  

  d_Aplus(5,0) = d_tempUnitNormal[1]*c_e*c_e;
  d_Aplus(5,1) = -d_tempUnitNormal[0]*c_e*c_e;
  d_Aplus(5,2) = 0;
  d_Aplus(5,3) = 0;
  d_Aplus(5,4) = 0;
  d_Aplus(5,5) = c_e; 
  d_Aplus(5,6) = 0; 
  d_Aplus(5,7) = 0; 
  
  d_Aplus(6,0) = d_tempUnitNormal[0]*c_e*c_e; 
  d_Aplus(6,1) = d_tempUnitNormal[1]*c_e*c_e; 
  d_Aplus(6,2) = 0;
  d_Aplus(6,3) = 0;
  d_Aplus(6,4) = 0;
  d_Aplus(6,5) = 0; 
  d_Aplus(6,6) = gamma*c_e; 
  d_Aplus(6,7) = 0;
  
  d_Aplus(7,0) = 0;
  d_Aplus(7,1) = 0;
  d_Aplus(7,2) = 0;
  d_Aplus(7,3) = d_tempUnitNormal[0];
  d_Aplus(7,4) = d_tempUnitNormal[1];
  d_Aplus(7,5) = 0; 
  d_Aplus(7,6) = 0; 
  d_Aplus(7,7) = chi*c_e;  

  d_Aminus(0,0) = -(d_tempUnitNormal[1]*d_tempUnitNormal[1] + gamma*d_tempUnitNormal[0]*d_tempUnitNormal[0])*c_e;
  d_Aminus(0,1) = (1 - gamma)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aminus(0,2) = 0;
  d_Aminus(0,3) = 0;
  d_Aminus(0,4) = 0;
  d_Aminus(0,5) = d_tempUnitNormal[1];
  d_Aminus(0,6) = gamma*gamma*d_tempUnitNormal[0]; 
  d_Aminus(0,7) = 0; 

  d_Aminus(1,0) = (1 - gamma)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aminus(1,1) = -(d_tempUnitNormal[0]*d_tempUnitNormal[0] + gamma*d_tempUnitNormal[1]*d_tempUnitNormal[1])*c_e;
  d_Aminus(1,2) = 0;
  d_Aminus(1,3) = 0;
  d_Aminus(1,4) = 0;
  d_Aminus(1,5) = -d_tempUnitNormal[0];
  d_Aminus(1,6) = gamma*gamma*d_tempUnitNormal[1]; 
  d_Aminus(1,7) = 0; 

  d_Aminus(2,0) = 0;
  d_Aminus(2,1) = 0;
  d_Aminus(2,2) = -c_e;
  d_Aminus(2,3) = -d_tempUnitNormal[1];
  d_Aminus(2,4) = d_tempUnitNormal[0];
  d_Aminus(2,5) = 0;
  d_Aminus(2,6) = 0;
  d_Aminus(2,7) = 0; 

  d_Aminus(3,0) = 0;
  d_Aminus(3,1) = 0;
  d_Aminus(3,2) = -d_tempUnitNormal[1]*c_e*c_e;
  d_Aminus(3,3) = -(d_tempUnitNormal[1]*d_tempUnitNormal[1] + chi*d_tempUnitNormal[0]*d_tempUnitNormal[0])*c_e;
  d_Aminus(3,4) = (1 - chi)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e;
  d_Aminus(3,5) = 0;
  d_Aminus(3,6) = 0;  
  d_Aminus(3,7) = chi*chi*d_tempUnitNormal[0]*c_e*c_e;  
  
  d_Aminus(4,0) = 0;
  d_Aminus(4,1) = 0;
  d_Aminus(4,2) = d_tempUnitNormal[0]*c_e*c_e;
  d_Aminus(4,3) = (1 - chi)*d_tempUnitNormal[0]*d_tempUnitNormal[1]*c_e ;
  d_Aminus(4,4) = -(d_tempUnitNormal[0]*d_tempUnitNormal[0] + chi*d_tempUnitNormal[1]*d_tempUnitNormal[1])*c_e;
  d_Aminus(4,5) = 0;
  d_Aminus(4,6) = 0; 
  d_Aminus(4,7) = chi*chi*d_tempUnitNormal[1]*c_e*c_e;  

  d_Aminus(5,0) = d_tempUnitNormal[1]*c_e*c_e;
  d_Aminus(5,1) = -d_tempUnitNormal[0]*c_e*c_e;
  d_Aminus(5,2) = 0;
  d_Aminus(5,3) = 0;
  d_Aminus(5,4) = 0;
  d_Aminus(5,5) = -c_e; 
  d_Aminus(5,6) = 0; 
  d_Aminus(5,7) = 0; 
  
  d_Aminus(6,0) = d_tempUnitNormal[0]*c_e*c_e; 
  d_Aminus(6,1) = d_tempUnitNormal[1]*c_e*c_e; 
  d_Aminus(6,2) = 0;
  d_Aminus(6,3) = 0;
  d_Aminus(6,4) = 0;
  d_Aminus(6,5) = 0; 
  d_Aminus(6,6) = -gamma*c_e; 
  d_Aminus(6,7) = 0;
  
  d_Aminus(7,0) = 0;
  d_Aminus(7,1) = 0;
  d_Aminus(7,2) = 0;
  d_Aminus(7,3) = d_tempUnitNormal[0];
  d_Aminus(7,4) = d_tempUnitNormal[1];
  d_Aminus(7,5) = 0; 
  d_Aminus(7,6) = 0; 
  d_Aminus(7,7) = -chi*c_e;  

 
      flux = 0.5*d_Aplus*d_lData + 0.5*d_Aminus*d_rData;
 
      flux *= data->getFaceArea();
      

      //Compute update coefficient
      if (!data->isPerturb()) {

        const CFreal FaceArea = data->getFaceArea();
        const CFreal MaxeValue = updateVS->getMaxEigenValue(&d_pdata[0], &d_tempUnitNormal[0]); 
        const CFreal k = MaxeValue*FaceArea;
        data->setUpdateCoeff(k);
      }

}
#endif

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "StegerWarmingMaxwellProjection2D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StegerWarmingMaxwellProjection2D_hh
