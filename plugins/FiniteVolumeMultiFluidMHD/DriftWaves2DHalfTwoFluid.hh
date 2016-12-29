//***************************************************************************
//
// Title: Template to create a Source term in multi-fluid MHD module
//
// Description: This class is composed of three files:
//   1.- DriftWaves2DHalfTwoFluid.hh: Is the header. Only declaration of the methods (functions) and members (variables used inside the class).
//   2.- DriftWaves2DHalfTwoFluid.cxx: Only implemented the provider of the class. It is used to be refered in the .CFcase file.
//   3.- DriftWaves2DHalfTwoFluid.ci: Is the main implementation in here.
//
//  This class implements a source term. In CoolFluid, we solve equations of the form:
//
//  \frac{\partial U}{\partial t} + \nabla \cdot F^c = \nabla \cdot F^d + S
//
//  where U is the array of variables, F^c is the convective fluxes, F^d is the diffusive fluxes and S the source.
//
//  In here we implement the source in the array source[] of function computeSource().
//
//  NOTE: Don't forget to multiply the source times the volume since it is a Finite Volume code
//
//
// In the following we will assume:
//
//  - The order of the equations are:
//     0 -- Bx
//     1 -- By
//     2 -- Bz
//     3 -- Ex
//     4 -- Ey
//     5 -- Ez
//     6 -- Psi
//     7 -- Phi
//     8 -- massConservation_electron
//     9 -- massConservation_ion
//     10 -- massConservation_neutral
//     11 -- momConservation_electron
//     12 -- momConservation_ion
//     13 -- momConservation_neutral
//     14 -- energyConservation_electron
//     15 -- energyConservation_ion
//     16 -- energyConservation_neutral
//
//  e.g., the source for the Bx equation is: source[0].
//
//  - The data is stored in an array called: _physicalData[] as follows:
//     _physicalData[0] -- Bx
//     _physicalData[1] -- By
//     _physicalData[2] -- Bz
//     _physicalData[3] -- Ex
//     _physicalData[4] -- Ey
//     _physicalData[5] -- Ez
//     _physicalData[6] -- Psi
//     _physicalData[7] -- Phi
//     _physicalData[8] -- RHO (total density)
//     _physicalData[9] -- x (coordinate of the state)
//     _physicalData[10] -- y (coordinate of the state)
//     _physicalData[11] -- y_electrons (partial density of ions) ============= firstDensity
//     _physicalData[12] -- y_ions (partial density of electron)
//     _physicalData[13] -- y_neutral (partial density of neutral)
//     _physicalData[14] -- U_electron ======================================= firstVelocity
//     _physicalData[15] -- V_electron
//     _physicalData[16] -- U_ion
//     _physicalData[17] -- V_ion
//     _physicalData[18] -- U_neutral
//     _physicalData[19] -- V_neutral
//     _physicalData[20] -- T_electron ======================================= firstTemperature
//     _physicalData[21] -- p_electron (pressure)
//     _physicalData[22] -- a_electron (speed of sound)
//     _physicalData[23] -- H_electron (total enthalpy)
//     _physicalData[24] -- T_ion
//     _physicalData[25] -- p_ion (pressure)
//     _physicalData[26] -- a_ion (speed of sound)
//     _physicalData[27] -- H_ion (total enthalpy)
//     _physicalData[28] -- T_neutral
//     _physicalData[29] -- p_neutral (pressure)
//     _physicalData[30] -- a_neutral (speed of sound)
//     _physicalData[31] -- H_neutral (total enthalpy)
//
//
// Author: Alejandro Alvarez
// Rev History: 2/2015
//
//***************************************************************************



#ifndef COOLFluiD_Numerics_FiniteVolume_DriftWaves2DHalfTwoFluid_hh
#define COOLFluiD_Numerics_FiniteVolume_DriftWaves2DHalfTwoFluid_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Common/SafePtr.hh"
#include "Framework/DataSocketSource.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/MathTypes.hh"
#include "Framework/VarSetTransformerT.hh"
#include "FiniteVolume/FluxData.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Source term for Three-fluid model considering fluids: ions + electrons + neutrals
 * variables
 *
 * @author Alejandro Alvarez
 * @author Isaac Alonso
 */
template <class UPDATEVAR>
class DriftWaves2DHalfTwoFluid : public ComputeSourceTermFVMCC {

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
      electricCharge = in->electricCharge;
      isCollisional = in->isCollisional;
      //parameters needed
    }
    
    //parameters needed
    bool isCollisional;
    CFreal electricCharge;
  };
  
  /// nested class defining a functor
  template <DeviceType DT, typename VS>
  class DeviceFunc {
  public:
    typedef VS MODEL;
    typedef DriftWaves2DHalfTwoFluid BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
    
    /// Compute the source term : implementation
    HOST_DEVICE void operator()(CFreal* state, VS* model, CFreal* source); 
   

  private:
    DeviceConfigOptions<NOTYPE>* m_dco;

    
    typename MathTypes<CFreal, DT, VS::DATASIZE>::VEC d_physicalData;
    typename MathTypes<CFreal, DT, 6>::VEC d_NonInducedEMField;
    typename MathTypes<CFreal, DT, 3>::VEC d_Btotal;
    typename MathTypes<CFreal, DT, 3>::VEC d_Etotal;

  };
  
  /// copy the local configuration options to the device
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {

    CFLog(VERBOSE, "DriftWaves2DHalfTwoFluid::copyConfigOptionsToDevice START1 \n \n");
    
    //CopyParameters
    CFreal electricCharge = getElectricCharge();
    bool isCollisional = getIsCollisional();
    
    CFLog(VERBOSE, "DriftWaves2DHalfTwoFluid::DeviceConfigOptions electricCharge = " << electricCharge  << "\n");
    CFLog(VERBOSE, "DriftWaves2DHalfTwoFluid::DeviceConfigOptions isCollisional = " << isCollisional  << "\n");
    CudaEnv::copyHost2Dev(&dco->electricCharge, &electricCharge, 1);
    CudaEnv::copyHost2Dev(&dco->isCollisional, &isCollisional, 1);


    CFLog(VERBOSE, "DriftWaves2DHalfTwoFluid::copyConfigOptionsToDevice END \n \n");
  }  
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    CFLog(VERBOSE, "DriftWaves2DHalfTwoFluid::copyConfigOptions \n");
    // consider to copy to constant memory
    dco->electricCharge = getElectricCharge();
    dco->isCollisional = getIsCollisional();

  } 



#endif

  /**
   * Constructor
   * @see ComputeSourceTermFVMCC
   */
  DriftWaves2DHalfTwoFluid(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DriftWaves2DHalfTwoFluid();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ComputeSourceTermFVMCC::configure(args);
    _sockets.template createSocketSink<RealVector>("nstates");
  }
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */
  virtual void computeSource(Framework::GeometricEntity *const element,
			     RealVector& source,
			     RealMatrix& jacobian);
  
  /**
   * Returns the DataSocket's that this command provides as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
  /**
  * Compute the electric Current
  */
  void computeEMField();

  
protected: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// handle to nodal states
  Framework::DataHandle<RealVector> _nstates;  
  
  /// socket for storing the Charge
  Framework::DataSocketSource<CFreal> socket_Qtot;
  
  /// socket for storing the Pressure gradient for the ions in y direction
  Framework::DataSocketSource<CFreal> socket_GradPyi;
  
  /// socket for storing the Pressure gradient for the elects in y direction
  Framework::DataSocketSource<CFreal> socket_GradPye;
  /// socket for storing the Ionization Rate
  //Framework::DataSocketSource<CFreal> socket_GammaRec;  
    
  /// pointer to the physical-chemical library
  //Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// Euler physical data
  RealVector _physicalData;

  /// vector to store temporary result
  RealVector _temp;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;
  
  /// array of temporary values
  RealMatrix _values;

  ///Non Induced Part of the electrocmagnetic Field
  RealVector _NonInducedEMField;
  
  ///Dummy vector for the gradients
  std::vector<RealVector*> _dummyGradients; 
  
  ///Vector storing the total magnetic Field
  RealVector _Btotal;
  
  ///Vector storing the total electric Field
  RealVector _Etotal;

  /// Option to change the electric charge
  CFreal _electricCharge;

  bool _isCollisional;

  CFreal getElectricCharge(){return _electricCharge;}
  bool getIsCollisional(){return _isCollisional;}

private:

  //options here
  
}; // end of class DriftWaves2DHalfTwoFluid

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA

template <class UPDATEVAR>
template <DeviceType DT, typename VS>
HOST_DEVICE void DriftWaves2DHalfTwoFluid<UPDATEVAR>::DeviceFunc<DT, VS>::operator()(CFreal* state, VS* model, CFreal* source) 
                                                           
{
 
  typename VS::UPDATE_VS* updateVS = model->getUpdateVS();
  updateVS->computePhysicalData(&state[0], &d_physicalData[0]);    
 								   

  for (CFuint i=0; i<6; i++){
    d_NonInducedEMField[i] = updateVS->getNonInducedEMField()[i]; 
  }
   
  
  //AAL: Here call all the functions needed to compute the source of Maxwell equations
  d_Etotal = 0;
  d_Btotal = 0;

  d_Btotal[XX] = d_physicalData[UPDATEVAR::PTERM::BX] + d_NonInducedEMField[0];
  d_Btotal[YY] = d_physicalData[UPDATEVAR::PTERM::BY] + d_NonInducedEMField[1];
  d_Btotal[ZZ] = d_physicalData[UPDATEVAR::PTERM::BZ] + d_NonInducedEMField[2];
  d_Etotal[XX] = d_physicalData[UPDATEVAR::PTERM::EX] + d_NonInducedEMField[3];
  d_Etotal[YY] = d_physicalData[UPDATEVAR::PTERM::EY] + d_NonInducedEMField[4];
  d_Etotal[ZZ] = d_physicalData[UPDATEVAR::PTERM::EZ] + d_NonInducedEMField[5];  



//Physical constants:
  const CFreal kB = updateVS->getK(); 				//Framework::PhysicalConsts::Boltzmann(); 
  const CFuint firstDensity = updateVS->getFirstSpecies(); 	//_varSet->getModel()->getFirstScalarVar(0);
  const CFreal qe = m_dco->electricCharge*(-1);                 // charge of electrons in Coulombs
  const CFreal qi = qe*(-1);                                    // charge of ions in Coulombs
  const CFreal mi = updateVS->getMolecularMass2(); 		//_varSet->getModel()->getMolecularMass2();  // Proton's mass [kg] 
								             //source:Standart Handbook for Electrical Engineerings
  const CFreal me = updateVS->getMolecularMass1();              // Electron's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal rho = d_physicalData[UPDATEVAR::PTERM::RHO];
  const CFreal rhoe = rho*d_physicalData[firstDensity]; 	//electrons density
  const CFreal rhoi = rho*d_physicalData[firstDensity + 1];     //ions density
 // printf("rho %e d_physicalData[firstDensity] %e d_physicalData[firstDensity+1] %e firstDensity %d \n", rho, d_physicalData[firstDensity],d_physicalData[firstDensity]+1, firstDensity);
  //std::cout << "rho  = " << rho  << "\n";getFirstDensity
  //std::cout << "rhoe = " << rhoe << "\n";
  //std::cout << "rhoi = " << rhoi << std::endl;
  const CFreal Qtot = qe*rhoe/me + qi*rhoi/mi;
  //printf("Qtot %e qe %e rhoe %e me %e qi %e rhoi %e mi %e \n", Qtot, qe, rhoe, me, qi, rhoi, mi);
  //std::cout << "Qtot = " << Qtot << std::endl;
  const CFuint firstVelocity = updateVS->getFirstVelocity(); 	
  const CFreal ue = d_physicalData[firstVelocity];
  const CFreal ve = d_physicalData[firstVelocity + 1];
  const CFreal we = d_physicalData[firstVelocity + 2];
  const CFreal ui = d_physicalData[firstVelocity + 3];
  const CFreal vi = d_physicalData[firstVelocity + 4];
  const CFreal wi = d_physicalData[firstVelocity + 5];
  //const CFreal un = d_physicalData[firstVelocity + 6];
  //const CFreal vn = d_physicalData[firstVelocity + 7];
  //const CFreal wn = d_physicalData[firstVelocity + 8];

  // Computing the electric current
  const CFreal Jx = qe*(rhoe/me)*ue + qi*(rhoi/mi)*ui;
  const CFreal Jy = qe*(rhoe/me)*ve + qi*(rhoi/mi)*vi;
  const CFreal Jz = qe*(rhoe/me)*we + qi*(rhoi/mi)*wi;
//  std::cout <<"Jx= " <<Jx <<"\n";
//  std::cout <<"ue= " <<ue <<"\n";
//  std::cout <<"ui= " <<ui <<"\n";
//  abort();

//AAL: Here goes the source of Maxwell equations
    /// MAXWELL
    const CFreal c_e = updateVS->getLightSpeed();
    const CFreal mu0 = updateVS->getPermeability();  
    const CFreal ovEpsilon = c_e*c_e*mu0;
   // printf("Jx %e Jy %e Jz %e ovEpsilon %e Qtot %e \n", Jx, Jy, Jz, ovEpsilon, Qtot);
    source[0] = 0.;             //x-Faraday's Law       	
    source[1] = 0.;          	//y-Faraday's Law		
    source[2] = 0.;          	//z-Faraday's Law
    source[3] = -Jx*ovEpsilon;	//x-Ampere's Law
    source[4] = -Jy*ovEpsilon;	//y-Ampere's Law
    source[5] = -Jz*ovEpsilon;  //z-Ampere's Law
    source[6] = 0.;		//divB
    source[7] = Qtot*ovEpsilon; //divE
     
  //AAL: Here the source for three-fluid continuity, momentum and energy equations
    //AAL: The following should be changed for the 3 Fluid case
    /// FLUID EQUATIONS
    //AAL: CONTINUITY
    source[8] = 0.;					// Electrons continuity equation
    source[9] = 0.;					// Ions continuity equation
      
    //AAL: MOMENTUM
    const CFuint firstTemperature = updateVS->getFirstTemperature(); 	//_varSet->getModel()->getFirstScalarVar(2);
    CFreal Te = d_physicalData[firstTemperature]; 			//electron temperature
    CFreal Ti = d_physicalData[firstTemperature + 4]; 			//ion temperature
    //printf("d_physicalData[%d] %e \n", firstTemperature+4, Ti);


    //electron, ion properties
    const CFreal ne = rhoe/me;		   	// number density electrons [m^-3]
    const CFreal ni = rhoi/mi;		   	// number density ions [m^-3]
    const CFreal pi = 3.141592653589793238462643383279502; //MathTools::MathConsts::CFrealPi(); 			//Pi number    
    const CFreal Epsilon0 = updateVS->getPermittivity();
    //std::cout <<"Epsilon0e = "<<Epsilon0<< "\n";
    //to calculate collision Frequency
 //   const CFreal r_di = sqrt(Epsilon0*kB*Ti/(ni*qe*qe)); //Debye length of ions [m]
 //   const CFreal r_de = sqrt(Epsilon0*kB*Te/(ne*qe*qe)); //Debye length of electrons [m]
 //   const CFreal r_deb = r_de*r_di/(sqrt(r_de*r_de + r_di*r_di)); //Debye length for plasma with several species

    // Collisional terms (Alex' way)
    // General terms
    const CFreal gamma_e = me/(kB*Te); 
    const CFreal gamma_i = mi/(kB*Ti);
    const CFreal mu_ie   = mi*me/(mi + me);
    const CFreal gamma_ie   = gamma_i*gamma_e/(gamma_i + gamma_e);

    //momentum equations - NO removed const
    CFreal emMomentumXe = qe*ne*(d_Etotal[XX] + ve*d_Btotal[ZZ] - we*d_Btotal[YY]);	//Electromagnetic momentum for electrons in X
    CFreal emMomentumXi = qi*ni*(d_Etotal[XX] + vi*d_Btotal[ZZ] - wi*d_Btotal[YY]);	//Electromagnetic momentum for ions in X
    CFreal emMomentumYe = qe*ne*(d_Etotal[YY] + we*d_Btotal[XX] - ue*d_Btotal[ZZ]);	//Electromagnetic momentum for electrons in Y
    CFreal emMomentumYi = qi*ni*(d_Etotal[YY] + wi*d_Btotal[XX] - ui*d_Btotal[ZZ]);	//Electromagnetic momentum for ions in Y
    CFreal emMomentumZe = qe*ne*(d_Etotal[ZZ] + ue*d_Btotal[YY] - ve*d_Btotal[XX]);	//Electromagnetic momentum for electrons in Z
    CFreal emMomentumZi = qi*ni*(d_Etotal[ZZ] + ui*d_Btotal[YY] - vi*d_Btotal[XX]);	//Electromagnetic momentum for ions in Z


    
    //DEBUG
    //printf("Epsilon0 %e \t kB %e \t ne %e \t qe %e \t ni %e \t qi %e \t Te %e \t Ti %e \n",
    //        Epsilon0,      kB,      ne,      qe,      ni,      qi,      Te,      Ti);
    //printf("ne %f \t mu_ie %f \t tau_minusOne_ei %f \t tau_minusOne_ie %f \t C_log %f \t Lambda_ie %f \t Debye_minusTwo %f \n",
    //        ne,      mu_ie,      tau_minusOne_ei,      tau_minusOne_ie,      C_log,      Lambda_ie,      Debye_minusTwo);
    //printf("emMomentumXe %f \t emMomentumXi %f \t emMomentumYe %f \t emMomentumYi %f \t emMomentumZe %f \t emMomentumZi %f \n",
    //        emMomentumXe,      emMomentumXi,      emMomentumYe,      emMomentumYi,      emMomentumZe,      emMomentumZi);
    //printf("collMomentumXe %f \t collMomentumYe %f \t collMomentumZe %f \t collMomentumXi %f \t collMomentumYi %f \t collMomentumZi %f \n",
    //        collMomentumXe,      collMomentumYe,      collMomentumZe,      collMomentumXi,      collMomentumYi,      collMomentumZi);


if(m_dco->isCollisional) {	

// Coulomb Collisions
    const CFreal Debye_minusTwo = ne*qe*qe/(Epsilon0*kB*Te) + ni*qi*qi/(Epsilon0*kB*Ti); 
    const CFreal Debye = sqrt(1/Debye_minusTwo);
    const CFreal Lambda_ie = 12*pi*Epsilon0/abs(qe*qi)*mu_ie/gamma_ie*Debye;
    const CFreal C_log = log(Lambda_ie); //Coulomb logarithm 
    const CFreal tau_minusOne_ie = 16*sqrt(pi)/3*ne*pow(gamma_ie/2,3./2.)*pow(qi*qe/(4*pi*Epsilon0*mu_ie),2.)*C_log;//ion collision frequency for collisions with electrons (target)
    const CFreal tau_minusOne_ei = 16*sqrt(pi)/3*ni*pow(gamma_ie/2,3./2.)*pow(qi*qe/(4*pi*Epsilon0*mu_ie),2.)*C_log;//electron collision frequency for collisions with ions (target)

    //collisional momentum:
    const CFreal collMomentumXe = -(ne*mu_ie*tau_minusOne_ei*(ue - ui));
    const CFreal collMomentumYe = -(ne*mu_ie*tau_minusOne_ei*(ve - vi));
    const CFreal collMomentumZe = -(ne*mu_ie*tau_minusOne_ei*(we - wi));
    const CFreal collMomentumXi = -(ni*mu_ie*tau_minusOne_ie*(ui - ue));
    const CFreal collMomentumYi = -(ni*mu_ie*tau_minusOne_ie*(vi - ve));
    const CFreal collMomentumZi = -(ni*mu_ie*tau_minusOne_ie*(wi - we));

      source[10] = emMomentumXe + collMomentumXe;   //Electrons X momentum
      source[11] = emMomentumYe + collMomentumYe;   //Electrons Y momentum
      source[12] = emMomentumZe + collMomentumZe;   //Electrons Z momentum

      source[13] = emMomentumXi + collMomentumXi;   //Ions X momentum
      source[14] = emMomentumYi + collMomentumYi;   //Ions Y momentum
      source[15] = emMomentumZi + collMomentumZi;   //Ions Z momentum


    //AAL: ENERGY
    // Computation of hydrodynamic pressure
//    const CFreal u = (rhoe*ue + rhoi*ui)/rho; 
//    const CFreal v = (rhoe*ve + rhoi*vi)/rho; 
//    const CFreal w = (rhoe*we + rhoi*wi)/rho;

    const CFreal workColle = (ue - ui)*collMomentumXe + (ve - vi)*collMomentumYe + (we - wi)*collMomentumZe; //Joule heating
    const CFreal workColli = (ui - ue)*collMomentumXi + (vi - ve)*collMomentumYi + (wi - we)*collMomentumZi;

    const CFreal heatColle = -3*kB*ne*(mu_ie/(me + mi))*tau_minusOne_ei*(Te - Ti); // collisional energy transfer
    const CFreal heatColli = -3*kB*ni*(mu_ie/(me + mi))*tau_minusOne_ie*(Ti - Te);

    const CFreal emEnergye = qe*(rhoe/me)*(ue*d_Etotal[XX] + ve*d_Etotal[YY] + we*d_Etotal[ZZ]); //electrons
    const CFreal emEnergyi = qi*(rhoi/mi)*(ui*d_Etotal[XX] + vi*d_Etotal[YY] + wi*d_Etotal[ZZ]); //ions

      source[16] = emEnergye + workColle + heatColle; //Electrons Energy Source terms
      source[17] = emEnergyi + workColli + heatColli; //Ions Energy

    }
    else { 

      source[10] = emMomentumXe;   //Electrons X momentum
      source[11] = emMomentumYe;   //Electrons Y momentum
      source[12] = emMomentumZe;   //Electrons Z momentum

      source[13] = emMomentumXi;   //Ions X momentum
      source[14] = emMomentumYi;   //Ions Y momentum
      source[15] = emMomentumZi;   //Ions Z momentum 

    const CFreal emEnergye = qe*(rhoe/me)*(ue*d_Etotal[XX] + ve*d_Etotal[YY] + we*d_Etotal[ZZ]); //electrons
    const CFreal emEnergyi = qi*(rhoi/mi)*(ui*d_Etotal[XX] + vi*d_Etotal[YY] + wi*d_Etotal[ZZ]); //ions

      source[16] = emEnergye;
      source[17] = emEnergyi;   
    }

      
    //AAL: ENERGY
    // Computation of hydrodynamic pressure
    //const CFreal u = (rhoe*ue + rhoi*ui)/rho;      
    //const CFreal v = (rhoe*ve + rhoi*vi)/rho; 
    //const CFreal w = (rhoe*we + rhoi*wi)/rho;

// DEBUG
/*
  printf("firstTemperature = %d \n", firstTemperature); //20
  printf("firstVelocity = %d \n", firstVelocity);  //14
  printf("firstSpecies = %d \n", firstDensity);  //12
  printf("k_B = %e \n", kB);
  printf("qe = %e \n", qe);
  printf("qi = %e \n", qi);
  printf("mi = %e \n", mi);
  printf("me = %e \n", me);
  printf("c_e = %e \n", c_e);
  printf("mu0 = %e \n", mu0); 
  printf("Epsilon0 = %e \n", Epsilon0);
  printf("NonInducedEMField = %e \t %e \t %e \t %e \t %e \t %e \n ", d_NonInducedEMField[0], d_NonInducedEMField[1], d_NonInducedEMField[2], d_NonInducedEMField[3], d_NonInducedEMField[4], d_NonInducedEMField[5]);
*/


  //DEBUG
  //for (CFuint i=0; i<VS::NBEQS; ++i){
  //  printf("source[%d] = %e \n", i, source[i]);
  //}					


}

#endif

/////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "DriftWaves2DHalfTwoFluid.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DriftWaves2DHalfTwoFluid_hh
