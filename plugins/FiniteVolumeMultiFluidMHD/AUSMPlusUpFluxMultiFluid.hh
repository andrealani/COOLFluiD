#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFluxMultiFluid_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFluxMultiFluid_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeMultiFluidMHD/AUSMFluxMultiFluid.hh"
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
 * This class computes the AUSM flux
 *
 * @author Alejandro Alvarez Laguna
 * @author Isaac Alonso
 *
 */
template <class UPDATEVAR>
class AUSMPlusUpFluxMultiFluid : public AUSMFluxMultiFluid<UPDATEVAR> {
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
      coeff = in->coeff;
      useMacCormackScaling = in->useMacCormackScaling;
      coeffKu = in->coeffKu;
      coeffKp = in->coeffKp;
      coeffSigma = in->coeffSigma;
      machInf[0] = in->machInf[0];       
      machInf[1] = in->machInf[1];
      beta = in->beta;
      fa = in->fa;
    }

    CFreal fa;
    CFreal coeff;
    bool useMacCormackScaling;
    CFreal coeffKu;
    CFreal coeffKp;
    CFreal coeffSigma;
    CFreal machInf[2];   //VS is not accesible from here, so this size is fixed (to run more than two species, change to the number desired)
    CFreal beta;
  };
  
  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef AUSMPlusUpFluxMultiFluid BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
    
    /// Compute the flux : implementation
    HOST_DEVICE void operator()(FluxData<VS>* data, VS* Model); 

    HOST_DEVICE CFreal mach2Plus(const CFreal mach) {return 0.25*pow(mach + 1.0, 2.0);}
   
    HOST_DEVICE CFreal mach2Min(const CFreal mach) {return -0.25*pow(mach - 1.0, 2.0);}

    HOST_DEVICE CFreal mach1Plus(const CFreal mach) {return 0.5*(mach + abs(mach));}
   
    HOST_DEVICE CFreal mach1Min(const CFreal mach) {return 0.5*(mach - abs(mach));}

    HOST_DEVICE virtual CFreal correctMachInfT(CFreal oldMach) const
    {
      return oldMach;
    }

  private:
    DeviceConfigOptions<NOTYPE>* m_dco;

    //VS* model;           //Model that can be accesed by all the functions inside the class
    //FluxData<VS>* data;  //The same for the data, in this case is not needed


    typename MathTypes<CFreal, DT, VS::DIM>::VEC d_normal;

    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_pdata;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_lData;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_psi_l;
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_psi_r;    
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_rData;
    
    //(to run more than two species, change to the number desired)
    typename MathTypes<CFreal, DT, 2>::VEC d_a12Vec;		//VS::NBSPECIES, but does not work because VarSetTransformerT does not have 
    typename MathTypes<CFreal, DT, 2>::VEC d_p12;		//NBSPECIES defined..
    typename MathTypes<CFreal, DT, 2>::VEC d_unL;		//In the MultiFluidMHDVarSetT DIM, DATASIZE and even NBSPECIES are defined
    typename MathTypes<CFreal, DT, 2>::VEC d_unR;		
    typename MathTypes<CFreal, DT, 2>::VEC d_mL;
    typename MathTypes<CFreal, DT, 2>::VEC d_mR;
    typename MathTypes<CFreal, DT, 2>::VEC d_machInf;
    typename MathTypes<CFreal, DT, 2>::VEC d_mflux12;

    typename MathTypes<CFreal, DT, 8>::VEC d_EMField_l;
    typename MathTypes<CFreal, DT, 8>::VEC d_EMField_r;
    typename MathTypes<CFreal, DT, 8>::MAT d_Aplus;
    typename MathTypes<CFreal, DT, 8>::MAT d_Aminus;
    typename MathTypes<CFreal, DT, 8>::VEC d_fluxEM;

  };
  
  /// copy the local configuration options to the device
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {


    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::copyConfigOptionsToDevice START9 \n \n");
    CFreal coeff = this->getMacCormackCoeff();
    bool useMacCormackScaling = this->getUseMacCormackScaling();
    CFreal coeffKu = getCoeffKu();
    CFreal coeffKp = getCoeffKp();
    CFreal coeffSigma = getCoeffSigma();

    CFreal* machInf = getMachInf();
    CFreal beta = getBeta();
    CFreal fa = getFa();

    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions coeff = " << coeff  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions useMacCormackScaling = " << useMacCormackScaling  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions coeffKu = " << coeffKu  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions coeffKp = " << coeffKp  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions coeffSigma = " << coeffSigma  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions machInf[0] = " << machInf[0]  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions machInf[1] = " << machInf[1]  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions beta = " << beta  << "\n");
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::DeviceConfigOptions fa = " << fa  << "\n");

    CudaEnv::copyHost2Dev(&dco->coeff, &coeff, 1);
    CudaEnv::copyHost2Dev(&dco->useMacCormackScaling, &useMacCormackScaling, 1);
    CudaEnv::copyHost2Dev(&dco->coeffKu, &coeffKu, 1);
    CudaEnv::copyHost2Dev(&dco->coeffKp, &coeffKp, 1);
    CudaEnv::copyHost2Dev(&dco->coeffSigma, &coeffSigma, 1);
    CudaEnv::copyHost2Dev(&dco->machInf[0], &machInf[0], 2);  //number of species
    CudaEnv::copyHost2Dev(&dco->beta, &beta, 1);
    CudaEnv::copyHost2Dev(&dco->fa, &fa, 1);
    

    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::copyConfigOptionsToDevice END \n \n");
  }  
  

  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    CFLog(VERBOSE, "AUSMPlusUpFluxMultiFluid::copyConfigOptions \n");
    // consider to copy to constant memory
    dco->useMacCormackScaling = this->getUseMacCormackScaling();
    dco->coeff = this->getMacCormackCoeff();
    dco->coeffKu = getCoeffKu();
    dco->coeffKp = getCoeffKp();
    dco->coeffSigma = getCoeffSigma();
    CFreal* machInf = getMachInf();
    dco->machInf[0] = machInf[0];
    dco->machInf[1] = machInf[1];      //number of species
    dco->beta = getBeta();
    dco->fa = getFa();
  }   

#endif







  /**
   * Constructor
   */
  AUSMPlusUpFluxMultiFluid(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusUpFluxMultiFluid() ;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the data from the supplied arguments.
   * @param args configuration arguments
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data
   */
  virtual void setup();

protected:
  
  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux();
  
  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux();
  
  /// Correct the given Mach number
  virtual CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach;
  }

  virtual CFreal getCoeffKu() {return m_coeffKu;}
  virtual CFreal getCoeffKp() {return m_coeffKp;}
  virtual CFreal getCoeffSigma() {return m_coeffSigma;}
  virtual CFreal* getMachInf() {return &m_machInf[0];}
  virtual CFreal getBeta() {return m_beta;}
  virtual CFreal getFa() {return m_fa;}
  

private:
  
  /// preconditioning coefficient 
  CFreal m_fa;
  
  /// user defined coefficient for Ku
  CFreal m_coeffKu;
  
  /// user defined coefficient for Kp
  CFreal m_coeffKp;
  
  /// user defined coefficient for sigma
  CFreal m_coeffSigma;

  /// mach infinity
  std::vector<CFreal> m_machInf;   //std::vector<CFreal>
  
  /// beta  coefficient
  CFreal m_beta;
  
}; // end of class AUSMPlusUpFluxMultiFluid

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA
/// functor that computes the flux
template <class UPDATEVAR>
template <DeviceType DT, typename VS>
HOST_DEVICE void AUSMPlusUpFluxMultiFluid<UPDATEVAR>::DeviceFunc<DT, VS>::operator()(FluxData<VS>* data, VS* model) 
{    
  
  typename VS::UPDATE_VS* updateVS = model->getUpdateVS();
  typename MathTypes<CFreal,DT,VS::NBEQS>::SLICEVEC flux(data->getResidual());   
   
  //Right physical data
  updateVS->computePhysicalData(data->getRstate(1), &d_rData[0]); 
  //Left physical data
  updateVS->computePhysicalData(data->getRstate(0), &d_lData[0]);

  typename MathTypes<CFreal,DT,VS::DIM>::SLICEVEC unitNormal(data->getUnitNormal());
  const CFreal coeff = (data->isOutward()) ? 1. : -1.;     
  d_normal = coeff*unitNormal;    


  CFuint dim = updateVS->getDim();
  bool is2DHalf = updateVS->getIs2DHalf();
  const CFuint firstTemperature = updateVS->getFirstTemperature(); 
  const CFuint firstVelocity = updateVS->getFirstVelocity();
  const CFuint firstSpecies = updateVS->getFirstSpecies();
  const CFuint nbSpecies = updateVS->getNbSpecies();

 //Gamma of the fluid
  const CFreal gamma = updateVS->getGamma(); 

 //Variables needed for the EM part
  const CFreal gammaB = updateVS->getDivBCleaningConst();   //Not get confused by the gamma of the fluid!!!
  const CFreal chi = updateVS->getDivECleaningConst();
  const CFreal c_e = updateVS->getLightSpeed();   
  CFreal factor1, factor2, factor3;
  const CFuint endEM = 8;

// DEBUG
//
//  printf("normal = %f \t %f \t %f \n", d_normal[0], d_normal[1], d_normal[2]);
//  printf("dim = %d \n", dim);
//  printf("firstTemperature = %d \n", firstTemperature);
//  printf("firstVelocity = %d \n", firstVelocity);
//  printf("firstSpecies = %d \n", firstSpecies);
//  printf("nbSpecies = %d \n", nbSpecies);
//  printf("gammaB = %f \n", gammaB);
//  printf("gamma = %f \n", gamma);
//  printf("chi = %f \n", chi);
//  printf("c_e = %f \n", c_e);
  
  
//  printf("useMacCormackScaling = %d \n", m_dco->useMacCormackScaling);
//  printf("coeff = %f \n", m_dco->coeff);
//  printf("coeffKu = %f \n", m_dco->coeffKu);
//  printf("coeffKp = %f \n", m_dco->coeffKp);
//  printf("coeffSigma = %f \n", m_dco->coeffSigma);
//  printf("machInf = %f \t %f \n", m_dco->machInf[0], m_dco->machInf[1]);
//  printf("beta = %f \n", m_dco->beta);
//  printf("fa = %f \n", m_dco->fa);


    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
      d_unL[ie] = 0.0;
      d_unR[ie] = 0.0;
      if(!is2DHalf){
        for (CFuint idim = 0; idim < dim; ++idim) {
          const CFuint currXX = firstVelocity + dim*ie + idim;
          d_unL[ie] += d_lData[currXX]*d_normal[idim];    
          d_unR[ie] += d_rData[currXX]*d_normal[idim];
        }
      }
      const CFuint dim2DHalf = 3;
      if(is2DHalf){
        for (CFuint idim = 0; idim < dim2DHalf; ++idim) {
          
          if( idim != 2) { // x,y-direction
            const CFuint currXX = firstVelocity + dim2DHalf*ie + idim;
            d_unL[ie] += d_lData[currXX]*d_normal[idim];   
            d_unR[ie] += d_rData[currXX]*d_normal[idim];
          }
        }
      }
    }
    
  //ComputeMassFlux
  for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
    
    d_machInf[ie] = m_dco->machInf[ie];    
    if (d_machInf[ie] <= 0.) {
      printf("AUSMPlusUpFluxMultiFluid requires machInf > 0.!!: change your input file!! \n");
    }
  }
  
  // the user must choose one of the following 3 ways of calculation of the
  // interface speed of sound, a12  
  d_a12Vec = 0.0;
  //switch(m_choiceA12) {
  //case 1:
  
  const CFreal gammaMinus1 = gamma - 1.0;  


  for (CFuint ie = 0; ie < nbSpecies; ++ie){

    const CFreal hL = d_lData[firstTemperature + 4*ie + 3];   
    const CFreal hR = d_rData[firstTemperature + 4*ie + 3];  

    // 1st Way of calculation of the interface speed of sound, a12,
    // suggested by Liou in his AIAA 2003-4116 paper

    const CFreal aCritL = sqrt( (2.0*gammaMinus1/(gamma+1.0))*hL);
    const CFreal aCritR = sqrt( (2.0*gammaMinus1/(gamma+1.0))*hR);

    const CFreal acL = (aCritL*aCritL)/max(aCritL, d_unL[ie]);
    const CFreal acR = (aCritR*aCritR)/max(aCritR, -d_unR[ie]);
    d_a12Vec[ie] = acL < acR ? acL : acR ; //min(acL, acR);		      //Array with the speed of sound of the different species
  }
  /*  break;     TODO: Implement the other cases (easy)
  case 2:
    computeSoundSpeed2();
    break;
  case 3:
    computeSoundSpeed3();
    break;
  case 4:
    computeSoundSpeed4();
    break;
  case 5:
    computeSoundSpeed5();
    break;
  default:
    computeSoundSpeed1();
    break;
  }
  */
  
  // calculation of the Mach number for the left and the right states
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {          
    d_mL[ie] = d_unL[ie]/d_a12Vec[ie];
    d_mR[ie] = d_unR[ie]/d_a12Vec[ie];
    const CFreal mL = d_mL[ie];
    const CFreal mR = d_mR[ie];
    const CFreal mBarSq = 0.5*(mL*mL + mR*mR);
      
    const CFreal mInf = correctMachInfT(d_machInf[ie]);   
    const CFreal A = mBarSq;  
     
    const CFreal B = mInf*mInf;
    const CFreal minAB = A < B ? A : B ;
    const CFreal Max1Min = 1.0 > minAB ? 1.0 : minAB;
    const CFreal mZero = sqrt(Max1Min);
    //const CFreal mZero = sqrt(min(1.0,    mBarSq < Inf*Inf ? mBarSq : Inf*Inf  ));             //max(mBarSq,  mInf*mInf)));
    //cf_assert(mZero <= 1.0);
    
    m_dco->fa = mZero * (2.0-mZero);
   
    //cf_assert(m_fa > 0.0);
    
    const CFreal M4Plus = (abs(mL) >= 1.0) ? mach1Plus(mL) :
      mach2Plus(mL)*(1.0 - 16.*m_dco->beta*mach2Min(mL));
    
    const CFreal M4Minus = (abs(mR) >= 1.0) ? mach1Min(mR) :
      mach2Min(mR)*(1.0 + 16.*m_dco->beta*mach2Plus(mR));

//   CFreal M4Plus = 0.0;
//   if (std::abs(mL) >= 1.0) {
//     M4Plus = 0.5 * (mL + std::abs(mL));
//   }
//   else {
//     M4Plus = 0.25*pow(mL + 1.0, 2.0) + beta*pow(mL*mL - 1.0, 2.0);
//   }

//   CFreal M4Minus = 0.0;
//   if (std::abs(mR) >= 1.0) {
//     M4Minus = 0.5 * (mR - std::abs(mR));
//   }
//   else {
//     M4Minus = -0.25*pow(mR - 1.0, 2.0) - beta*pow(mR*mR - 1.0, 2.0);
//   }

    const CFreal rhoL = d_lData[UPDATEVAR::PTERM::RHO]*d_lData[firstSpecies + ie];   //(*d_lData)
    const CFreal rhoR = d_rData[UPDATEVAR::PTERM::RHO]*d_rData[firstSpecies + ie];
    const CFreal pL = d_lData[firstTemperature + 4*ie + 1];
    const CFreal pR = d_rData[firstTemperature + 4*ie + 1];
    const CFreal rhoa2 = 0.5*(rhoL + rhoR)*d_a12Vec[ie]*d_a12Vec[ie];
    const CFreal mP = (-m_dco->coeffKp/m_dco->fa) * max(1.0 - m_dco->coeffSigma*mBarSq, 0.0)*
      (pR-pL)/rhoa2;
  
    // calculation of the Mach number at the interface
    const CFreal m12 = M4Plus + M4Minus + mP;
    // calculation of the mass flux at the interface
    d_mflux12[ie] = (m12 > 0.0) ? d_a12Vec[ie] * m12 * rhoL : d_a12Vec[ie] * m12 * rhoR;   
    
  }


  //This part is for the MultiFluid euler equations

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    
    const CFreal alpha = (3.0/16.0) * (-4.0 + 5.0*m_dco->fa*m_dco->fa);
    const CFreal mL = d_mL[ie];
    const CFreal mR = d_mR[ie];
    const CFreal P5Plus = (abs(mL) >= 1.0) ? mach1Plus(mL)/mL :
      mach2Plus(mL)*((2.0-mL)-16.*alpha*mL*mach2Min(mL));

    const CFreal P5Minus = (abs(mR) >= 1.0) ? mach1Min(mR)/mR :
      mach2Min(mR)*((-2.0-mR)+16.*alpha*mR*mach2Plus(mR));
   
  // CFreal P5Plus = 0.0;
//   if (std::abs(mL) >= 1.0) {
//     P5Plus = 0.5 * (mL + std::abs(mL))/mL;
//   }
//   else {
//     P5Plus = 0.25*pow(mL + 1.0, 2.0)*(2.0-mL) + alpha*mL*pow(mL*mL - 1.0, 2.0);
//   }
  
//   // CFreal P5Minus = getP5Min(mR);
//   CFreal P5Minus = 0.0;
//   if (std::abs(mR) >= 1.0) {
//     P5Minus = 0.5 * (mR - std::abs(mR))/mR;
//   }
//   else {
//     P5Minus = 0.25*pow(mR - 1.0, 2.0)*(2.0+mR) - alpha*mR*pow(mR*mR - 1.0, 2.0);
//   }
  
    const CFreal rhoL = d_lData[UPDATEVAR::PTERM::RHO]*d_lData[firstSpecies + ie];   //(*d_lData)
    const CFreal rhoR = d_rData[UPDATEVAR::PTERM::RHO]*d_rData[firstSpecies + ie];
    const CFreal pL = d_lData[firstTemperature + 4*ie + 1];
    const CFreal pR = d_rData[firstTemperature + 4*ie + 1];
    const CFreal pU = -m_dco->coeffKu * P5Plus * P5Minus *
      (rhoL + rhoR) * m_dco->fa * d_a12Vec[ie]*(d_unR[ie]-d_unL[ie]);

  // calculation of the pressure flux at the interface
    d_p12[ie] = P5Plus*pL + P5Minus*pR + pU;

//    cout << "PressureFlux ie = "<< ie << endl;
//    cout << "PressureFlux d_p12[ie] = "<< d_p12[ie] << endl;

  }


  // calculation of the dimensional numerical fluxes at the interface
  // computeMassFluxImpl
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    
    if (d_mflux12[ie] > 0.0) { // LEFT CASE

      d_psi_l[ie] = 1;      
      flux[endEM + ie] = d_mflux12[ie];
      
      //loop to set the velocities inside Psi
      if(!is2DHalf) {
        for (CFuint i = 0; i < dim; ++i) {
          d_psi_l[nbSpecies + ie*dim + i] =  d_lData[firstVelocity + ie*dim + i];
          flux[endEM + nbSpecies + ie*dim + i] = d_mflux12[ie]*d_psi_l[nbSpecies + ie*dim + i] + d_p12[ie]*d_normal[i];
        }
      }
      else { // is 2DHalf
        const CFuint dim2DHalf = 3;
        for (CFuint i = 0; i < dim2DHalf; ++i) {
          if( i != 2) { // x,y-direction
            d_psi_l[nbSpecies + ie*dim2DHalf + i] =  d_lData[firstVelocity + ie*dim2DHalf + i];
            flux[endEM + nbSpecies + ie*dim2DHalf + i] = d_mflux12[ie]*d_psi_l[nbSpecies + ie*dim2DHalf + i] + d_p12[ie]*d_normal[i];
          }
          else { // z-direction
            CFreal nz = 0;
            d_psi_l[nbSpecies + ie*dim2DHalf + i] =  d_lData[firstVelocity + ie*dim2DHalf + i];
            flux[endEM + nbSpecies + ie*dim2DHalf + i] = d_mflux12[ie]*d_psi_l[nbSpecies + ie*dim2DHalf + i] + d_p12[ie]*nz;
          }
        }
      }
      //enthalpy
      if(!is2DHalf) { //Default case
        d_psi_l[nbSpecies + nbSpecies*dim + ie] = d_lData[firstTemperature + 4*ie + 3];
        flux[endEM + nbSpecies + nbSpecies*dim + ie] = d_mflux12[ie]*d_psi_l[nbSpecies + nbSpecies*dim + ie];
      }
      else { // when it is 2D half
        d_psi_l[nbSpecies + nbSpecies*3 + ie] = d_lData[firstTemperature + 4*ie + 3];
        flux[endEM + nbSpecies + nbSpecies*3 + ie] = d_mflux12[ie]*d_psi_l[nbSpecies + nbSpecies*3 + ie];
      }
    }

  
    else { // RIGHT CASE
      
      d_psi_r[ie] = 1;      
      flux[endEM + ie] = d_mflux12[ie];
      
      //loop to set the velocities inside Psi
      if(!is2DHalf) {
        for (CFuint i = 0; i < dim; ++i) {
          d_psi_r[nbSpecies + ie*dim + i] =  d_rData[firstVelocity + ie*dim + i];
          flux[endEM + nbSpecies + ie*dim + i] = d_mflux12[ie]*d_psi_r[nbSpecies + ie*dim + i] + d_p12[ie]*d_normal[i];
        }
      }
      else { // is 2DHalf
        const CFuint dim2DHalf = 3;
        for (CFuint i = 0; i < dim2DHalf; ++i) {
          if( i != 2) { // x,y-direction
            d_psi_r[nbSpecies + ie*dim2DHalf + i] =  d_rData[firstVelocity + ie*dim2DHalf + i];
            flux[endEM + nbSpecies + ie*dim2DHalf + i] = d_mflux12[ie]*d_psi_r[nbSpecies + ie*dim2DHalf + i] + d_p12[ie]*d_normal[i];
          }
          else { // z-direction
            CFreal nz = 0;
            d_psi_r[nbSpecies + ie*dim2DHalf + i] =  d_rData[firstVelocity + ie*dim2DHalf + i];
            flux[endEM + nbSpecies + ie*dim2DHalf + i] = d_mflux12[ie]*d_psi_r[nbSpecies + ie*dim2DHalf + i] + d_p12[ie]*nz;
          }
        }
      }
      //total enthalpy
      if(!is2DHalf) { //Default case
        d_psi_r[nbSpecies + nbSpecies*dim + ie] = d_rData[firstTemperature + 4*ie + 3];
        flux[endEM + nbSpecies + nbSpecies*dim + ie] = d_mflux12[ie]*d_psi_r[nbSpecies + nbSpecies*dim + ie] ;
      }
      else { // when it is 2D half
        d_psi_r[nbSpecies + nbSpecies*3 + ie] = d_rData[firstTemperature + 4*ie + 3];
        flux[endEM + nbSpecies + nbSpecies*3 + ie] = d_mflux12[ie]*d_psi_r[nbSpecies + nbSpecies*3 + ie] ;
      }
    }

  }
  
  


  ///flux splitting scheme for Maxwell's Equations
  //loop to set the electromagnetic variables 
  for (CFuint i = 0; i < endEM; ++i){
    d_EMField_l[i] = d_lData[i];
    d_EMField_r[i] = d_rData[i];    
  }


  if(m_dco->useMacCormackScaling){
    factor1 = 1;
    factor2 = c_e*c_e;
    factor3 = m_dco->coeff*c_e*c_e;
  }
  else{
    factor1 = c_e;
    factor2 = c_e;
    factor3 = c_e;
  }
  
  if(dim == 2){
    d_Aplus(0,0) = (d_normal[1]*d_normal[1] + gammaB*d_normal[0]*d_normal[0])*factor1; //flag2
    d_Aplus(0,1) = (gammaB - 1)*d_normal[0]*d_normal[1]*factor1; //flag2
    d_Aplus(0,2) = 0;
    d_Aplus(0,3) = 0;
    d_Aplus(0,4) = 0;
    d_Aplus(0,5) = d_normal[1];
    d_Aplus(0,6) = gammaB*gammaB*d_normal[0]; 
    d_Aplus(0,7) = 0; 

    d_Aplus(1,0) = (gammaB - 1)*d_normal[0]*d_normal[1]*factor1; //flag2;
    d_Aplus(1,1) = (d_normal[0]*d_normal[0] + gammaB*d_normal[1]*d_normal[1])*factor1; //flag2
    d_Aplus(1,2) = 0;
    d_Aplus(1,3) = 0;
    d_Aplus(1,4) = 0;
    d_Aplus(1,5) = -d_normal[0];
    d_Aplus(1,6) = gammaB*gammaB*d_normal[1]; 
    d_Aplus(1,7) = 0; 

    d_Aplus(2,0) = 0;
    d_Aplus(2,1) = 0;
    d_Aplus(2,2) = factor1; //flag2
    d_Aplus(2,3) = -d_normal[1];
    d_Aplus(2,4) = d_normal[0];
    d_Aplus(2,5) = 0;
    d_Aplus(2,6) = 0; 
    d_Aplus(2,7) = 0; 

    d_Aplus(3,0) = 0;
    d_Aplus(3,1) = 0;
    d_Aplus(3,2) = -d_normal[1]*c_e*c_e;
    d_Aplus(3,3) = (d_normal[1]*d_normal[1] + chi*d_normal[0]*d_normal[0])*factor2; //flag
    d_Aplus(3,4) = (chi - 1)*d_normal[0]*d_normal[1]*factor2; //flag
    d_Aplus(3,5) = 0;
    d_Aplus(3,6) = 0;  
    d_Aplus(3,7) = chi*chi*d_normal[0]*c_e*c_e;  
    
    d_Aplus(4,0) = 0;
    d_Aplus(4,1) = 0;
    d_Aplus(4,2) = d_normal[0]*c_e*c_e;
    d_Aplus(4,3) = (chi -1)*d_normal[0]*d_normal[1]*factor2; //flag
    d_Aplus(4,4) = (d_normal[0]*d_normal[0] + chi*d_normal[1]*d_normal[1])*factor2; //flag
    d_Aplus(4,5) = 0;
    d_Aplus(4,6) = 0; 
    d_Aplus(4,7) = chi*chi*d_normal[1]*c_e*c_e;  

    d_Aplus(5,0) = d_normal[1]*c_e*c_e;
    d_Aplus(5,1) = -d_normal[0]*c_e*c_e;
    d_Aplus(5,2) = 0;
    d_Aplus(5,3) = 0;
    d_Aplus(5,4) = 0;
    d_Aplus(5,5) = factor2; //flag 
    d_Aplus(5,6) = 0; 
    d_Aplus(5,7) = 0; 
    
    d_Aplus(6,0) = d_normal[0]*c_e*c_e;
    d_Aplus(6,1) = d_normal[1]*c_e*c_e;
    d_Aplus(6,2) = 0;
    d_Aplus(6,3) = 0;
    d_Aplus(6,4) = 0;
    d_Aplus(6,5) = 0;  
    d_Aplus(6,6) = gammaB*factor3;
    d_Aplus(6,7) = 0;
    
    d_Aplus(7,0) = 0;
    d_Aplus(7,1) = 0;
    d_Aplus(7,2) = 0;
    d_Aplus(7,3) = d_normal[0];
    d_Aplus(7,4) = d_normal[1];
    d_Aplus(7,5) = 0; 
    d_Aplus(7,6) = 0; 
    d_Aplus(7,7) = chi*factor1; //flag2  
  }
  if(dim == 3){
    d_Aplus(0,0) = (1 + d_normal[0]*d_normal[0]*(gammaB - 1));
    d_Aplus(0,1) = (gammaB - 1)*d_normal[0]*d_normal[1];
    d_Aplus(0,2) = (gammaB - 1)*d_normal[0]*d_normal[2];
    d_Aplus(0,3) = 0;
    d_Aplus(0,4) = -d_normal[2];
    d_Aplus(0,5) = d_normal[1];
    d_Aplus(0,6) = gammaB*gammaB*d_normal[0]; 
    d_Aplus(0,7) = 0; 

    d_Aplus(1,0) = (gammaB - 1)*d_normal[0]*d_normal[1];
    d_Aplus(1,1) = (1 + d_normal[1]*d_normal[1]*(gammaB - 1));
    d_Aplus(1,2) = (gammaB - 1)*d_normal[1]*d_normal[2];
    d_Aplus(1,3) = d_normal[2];
    d_Aplus(1,4) = 0;
    d_Aplus(1,5) = -d_normal[0];
    d_Aplus(1,6) = gammaB*gammaB*d_normal[1]; 
    d_Aplus(1,7) = 0; 

    d_Aplus(2,0) = (gammaB - 1)*d_normal[0]*d_normal[2];
    d_Aplus(2,1) = (gammaB - 1)*d_normal[1]*d_normal[2];
    d_Aplus(2,2) = (1 + d_normal[2]*d_normal[2]*(gammaB - 1));
    d_Aplus(2,3) = -d_normal[1];
    d_Aplus(2,4) = d_normal[0];
    d_Aplus(2,5) = 0;
    d_Aplus(2,6) = gammaB*gammaB*d_normal[2]; 
    d_Aplus(2,7) = 0; 

    d_Aplus(3,0) = 0; 
    d_Aplus(3,1) = d_normal[2]*c_e*c_e; 
    d_Aplus(3,2) = -d_normal[1]*c_e*c_e;
    d_Aplus(3,3) = (1 + d_normal[0]*d_normal[0]*(chi - 1))*c_e*c_e;
    d_Aplus(3,4) = (chi - 1)*d_normal[0]*d_normal[1]*c_e*c_e;
    d_Aplus(3,5) = (chi - 1)*d_normal[0]*d_normal[2]*c_e*c_e; 
    d_Aplus(3,6) = 0; 
    d_Aplus(3,7) = chi*chi*d_normal[0]*c_e*c_e;

    d_Aplus(4,0) = -d_normal[2]*c_e*c_e;
    d_Aplus(4,1) = 0;
    d_Aplus(4,2) = d_normal[0]*c_e*c_e;
    d_Aplus(4,3) = (chi - 1)*d_normal[0]*d_normal[1]*c_e*c_e;
    d_Aplus(4,4) = (1 + d_normal[1]*d_normal[1]*(chi - 1))*c_e*c_e;
    d_Aplus(4,5) = (chi - 1)*d_normal[1]*d_normal[2]*c_e*c_e;
    d_Aplus(4,6) = 0; 
    d_Aplus(4,7) = chi*chi*d_normal[1]*c_e*c_e; 

    d_Aplus(5,0) = d_normal[1]*c_e*c_e;
    d_Aplus(5,1) = -d_normal[0]*c_e*c_e;
    d_Aplus(5,2) = 0;
    d_Aplus(5,3) = (chi - 1)*d_normal[0]*d_normal[2]*c_e*c_e;
    d_Aplus(5,4) = (chi - 1)*d_normal[1]*d_normal[2]*c_e*c_e;
    d_Aplus(5,5) = (1 + d_normal[2]*d_normal[2]*(chi - 1))*c_e*c_e;
    d_Aplus(5,6) = 0; 
    d_Aplus(5,7) = chi*chi*d_normal[2]*c_e*c_e; 
    
    d_Aplus(6,0) = d_normal[0]*c_e*c_e;
    d_Aplus(6,1) = d_normal[1]*c_e*c_e;
    d_Aplus(6,2) = d_normal[2]*c_e*c_e;
    d_Aplus(6,3) = 0;
    d_Aplus(6,4) = 0;
    d_Aplus(6,5) = 0; 
    d_Aplus(6,6) = gammaB*c_e*c_e; 
    d_Aplus(6,7) = 0;
    
    d_Aplus(7,0) = 0;
    d_Aplus(7,1) = 0;
    d_Aplus(7,2) = 0;
    d_Aplus(7,3) = d_normal[0];
    d_Aplus(7,4) = d_normal[1];
    d_Aplus(7,5) = d_normal[2]; 
    d_Aplus(7,6) = 0; 
    d_Aplus(7,7) = chi;
  }



  


  if(m_dco->useMacCormackScaling){
    factor1 = 1;
    factor2 = c_e*c_e;
    factor3 = m_dco->coeff*c_e*c_e;
  }
  else{
    factor1 = c_e;
    factor2 = c_e;
    factor3 = c_e;
  }

  if(dim == 2){
    d_Aminus(0,0) = -(d_normal[1]*d_normal[1] + gammaB*d_normal[0]*d_normal[0])*factor1; //flag2
    d_Aminus(0,1) = (1 - gammaB)*d_normal[0]*d_normal[1]*factor1; //flag2
    d_Aminus(0,2) = 0;
    d_Aminus(0,3) = 0;
    d_Aminus(0,4) = 0;
    d_Aminus(0,5) = d_normal[1];
    d_Aminus(0,6) = gammaB*gammaB*d_normal[0]; 
    d_Aminus(0,7) = 0; 

    d_Aminus(1,0) = (1 - gammaB)*d_normal[0]*d_normal[1]*factor1; //flag2
    d_Aminus(1,1) = -(d_normal[0]*d_normal[0] + gammaB*d_normal[1]*d_normal[1])*factor1; //flag2
    d_Aminus(1,2) = 0;
    d_Aminus(1,3) = 0;
    d_Aminus(1,4) = 0;
    d_Aminus(1,5) = -d_normal[0];
    d_Aminus(1,6) = gammaB*gammaB*d_normal[1]; 
    d_Aminus(1,7) = 0; 

    d_Aminus(2,0) = 0;
    d_Aminus(2,1) = 0;
    d_Aminus(2,2) = factor1; //flag2
    d_Aminus(2,3) = -d_normal[1];
    d_Aminus(2,4) = d_normal[0];
    d_Aminus(2,5) = 0;
    d_Aminus(2,6) = 0;
    d_Aminus(2,7) = 0; 

    d_Aminus(3,0) = 0;
    d_Aminus(3,1) = 0;
    d_Aminus(3,2) = -d_normal[1]*c_e*c_e;
    d_Aminus(3,3) = -(d_normal[1]*d_normal[1] + chi*d_normal[0]*d_normal[0])*factor2; //flag
    d_Aminus(3,4) = (1 - chi)*d_normal[0]*d_normal[1]*factor2; //flag
    d_Aminus(3,5) = 0;
    d_Aminus(3,6) = 0;  
    d_Aminus(3,7) = chi*chi*d_normal[0]*c_e*c_e;  
    
    d_Aminus(4,0) = 0;
    d_Aminus(4,1) = 0;
    d_Aminus(4,2) = d_normal[0]*c_e*c_e;
    d_Aminus(4,3) = (1 - chi)*d_normal[0]*d_normal[1]*factor2; //flag
    d_Aminus(4,4) = -(d_normal[0]*d_normal[0] + chi*d_normal[1]*d_normal[1])*factor2; //flag
    d_Aminus(4,5) = 0;
    d_Aminus(4,6) = 0; 
    d_Aminus(4,7) = chi*chi*d_normal[1]*c_e*c_e;  

    d_Aminus(5,0) = d_normal[1]*c_e*c_e;
    d_Aminus(5,1) = -d_normal[0]*c_e*c_e;
    d_Aminus(5,2) = 0;
    d_Aminus(5,3) = 0;
    d_Aminus(5,4) = 0;
    d_Aminus(5,5) = -1*factor2; //flag 
    d_Aminus(5,6) = 0; 
    d_Aminus(5,7) = 0; 
    
    d_Aminus(6,0) = d_normal[0]*c_e*c_e;
    d_Aminus(6,1) = d_normal[1]*c_e*c_e;
    d_Aminus(6,2) = 0;
    d_Aminus(6,3) = 0;
    d_Aminus(6,4) = 0;
    d_Aminus(6,5) = 0;
    d_Aminus(6,6) = -gammaB*factor3;
    d_Aminus(6,7) = 0;
    
    d_Aminus(7,0) = 0;
    d_Aminus(7,1) = 0;
    d_Aminus(7,2) = 0;
    d_Aminus(7,3) = d_normal[0];
    d_Aminus(7,4) = d_normal[1];
    d_Aminus(7,5) = 0; 
    d_Aminus(7,6) = 0; 
    d_Aminus(7,7) = -chi*factor1; //flag2  
  }
  if(dim == 3){
    d_Aminus(0,0) = -(1 + d_normal[0]*d_normal[0]*(gammaB - 1));
    d_Aminus(0,1) = -(gammaB - 1)*d_normal[0]*d_normal[1];
    d_Aminus(0,2) = -(gammaB - 1)*d_normal[0]*d_normal[2];
    d_Aminus(0,3) = 0;
    d_Aminus(0,4) = -d_normal[2];
    d_Aminus(0,5) = d_normal[1];
    d_Aminus(0,6) = gammaB*gammaB*d_normal[0]; 
    d_Aminus(0,7) = 0; 

    d_Aminus(1,0) = -(gammaB - 1)*d_normal[0]*d_normal[1];
    d_Aminus(1,1) = -(1 + d_normal[1]*d_normal[1]*(gammaB - 1));
    d_Aminus(1,2) = -(gammaB - 1)*d_normal[1]*d_normal[2];
    d_Aminus(1,3) = d_normal[2];
    d_Aminus(1,4) = 0;
    d_Aminus(1,5) = -d_normal[0];
    d_Aminus(1,6) = gammaB*gammaB*d_normal[1]; 
    d_Aminus(1,7) = 0; 

    d_Aminus(2,0) = -(gammaB - 1)*d_normal[0]*d_normal[2];
    d_Aminus(2,1) = -(gammaB - 1)*d_normal[1]*d_normal[2];
    d_Aminus(2,2) = -(1 + d_normal[2]*d_normal[2]*(gammaB - 1));
    d_Aminus(2,3) = -d_normal[1];
    d_Aminus(2,4) = d_normal[0];
    d_Aminus(2,5) = 0;
    d_Aminus(2,6) = gammaB*gammaB*d_normal[2]; 
    d_Aminus(2,7) = 0; 

    d_Aminus(3,0) = 0; 
    d_Aminus(3,1) = d_normal[2]*c_e*c_e; 
    d_Aminus(3,2) = -d_normal[1]*c_e*c_e;
    d_Aminus(3,3) = -(1 + d_normal[0]*d_normal[0]*(chi - 1))*c_e*c_e;
    d_Aminus(3,4) = -(chi - 1)*d_normal[0]*d_normal[1]*c_e*c_e;
    d_Aminus(3,5) = -(chi - 1)*d_normal[0]*d_normal[2]*c_e*c_e; 
    d_Aminus(3,6) = 0; 
    d_Aminus(3,7) = chi*chi*d_normal[0]*c_e*c_e;

    d_Aminus(4,0) = -d_normal[2]*c_e*c_e;
    d_Aminus(4,1) = 0;
    d_Aminus(4,2) = d_normal[0]*c_e*c_e;
    d_Aminus(4,3) = -(chi - 1)*d_normal[0]*d_normal[1]*c_e*c_e;
    d_Aminus(4,4) = -(1 + d_normal[1]*d_normal[1]*(chi - 1))*c_e*c_e;
    d_Aminus(4,5) = -(chi - 1)*d_normal[1]*d_normal[2]*c_e*c_e;
    d_Aminus(4,6) = 0; 
    d_Aminus(4,7) = chi*chi*d_normal[1]*c_e*c_e; 

    d_Aminus(5,0) = d_normal[1]*c_e*c_e;
    d_Aminus(5,1) = -d_normal[0]*c_e*c_e;
    d_Aminus(5,2) = 0;
    d_Aminus(5,3) = -(chi - 1)*d_normal[0]*d_normal[2]*c_e*c_e;
    d_Aminus(5,4) = -(chi - 1)*d_normal[1]*d_normal[2]*c_e*c_e;
    d_Aminus(5,5) = -(1 + d_normal[2]*d_normal[2]*(chi - 1))*c_e*c_e;
    d_Aminus(5,6) = 0; 
    d_Aminus(5,7) = chi*chi*d_normal[2]*c_e*c_e; 
    
    d_Aminus(6,0) = d_normal[0]*c_e*c_e;
    d_Aminus(6,1) = d_normal[1]*c_e*c_e;
    d_Aminus(6,2) = d_normal[2]*c_e*c_e;
    d_Aminus(6,3) = 0;
    d_Aminus(6,4) = 0;
    d_Aminus(6,5) = 0; 
    d_Aminus(6,6) = -gammaB*c_e*c_e; 
    d_Aminus(6,7) = 0;
    
    d_Aminus(7,0) = 0;
    d_Aminus(7,1) = 0;
    d_Aminus(7,2) = 0;
    d_Aminus(7,3) = d_normal[0];
    d_Aminus(7,4) = d_normal[1];
    d_Aminus(7,5) = d_normal[2]; 
    d_Aminus(7,6) = 0; 
    d_Aminus(7,7) = -chi;  
  }  

  d_fluxEM = 0.5*d_Aplus*d_EMField_l + 0.5*d_Aminus*d_EMField_r;   
    
  for(CFuint iem = 0; iem < endEM; ++iem){
    flux[iem] = d_fluxEM[iem];   
  }    
  
  



  const CFreal FaceArea = data->getFaceArea();
  const CFreal MaxeValue = updateVS->getMaxEigenValue(&d_pdata[0], &d_normal[0]);  
  const CFreal k = MaxeValue*FaceArea;     
  if (!data->isPerturb()) {
    //printf("UpdateCoeff k = %f \n", k);
    data->setUpdateCoeff(k);
  }
  
 
  flux *= FaceArea;   
   
  
/*
  if(is2DHalf){
    printf("%.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %e \n",
     flux[0],flux[1], flux[2], flux[3],flux[4],flux[5],flux[6],flux[7],flux[8],flux[9],flux[10],flux[11],flux[12],flux[13],flux[14],flux[15], k);
  }

  else{
    printf("%.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e \t %e \n",
     flux[0],flux[1], flux[2], flux[3],flux[4],flux[5],flux[6],flux[7],flux[8],flux[9],flux[10],flux[11],flux[12],flux[13],flux[14],flux[15], k);
  }
*/  

}


#endif //CF_HAVE_CUDA

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMPlusUpFluxMultiFluid.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFluxMultiFluid_hh
