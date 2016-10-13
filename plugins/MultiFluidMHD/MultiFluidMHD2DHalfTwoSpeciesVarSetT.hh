#ifndef COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHD2DHalfTwoSpeciesVarSetT_hh
#define COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHD2DHalfTwoSpeciesVarSetT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MultiScalarVarSetBase.hh"
#include "EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an MultiFluidMHD physical model. With 2 species in 2DHalf
 *
 * @author Andrea Lani
 */

template <typename BASE>      
class MultiFluidMHD2DHalfTwoSpeciesVarSetT : public BASE {
public: // classes

  //Hard-coded for each VarSet
  enum {DIM=2, NBEQS=18, DATASIZE=32, NBSPECIES=2};




  typedef Framework::MultiScalarTerm<EulerMFMHDTerm> PTERM;  
  typedef MultiFluidMHD2DHalfTwoSpeciesVarSetT EULERSET;
  /**
   * Constructor
   * @see MultiFluidMHDModel
   */
  

  HOST_DEVICE MultiFluidMHD2DHalfTwoSpeciesVarSetT(EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* dco) :
    m_dco(dco) {} 


  HOST_DEVICE MultiFluidMHD2DHalfTwoSpeciesVarSetT() {}


  HOST_DEVICE void setModelData(EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* dco) {m_dco = dco;}  


  /**
   * Default destructor
   */
  HOST_DEVICE virtual ~MultiFluidMHD2DHalfTwoSpeciesVarSetT() {};

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  //virtual void setup();
  
  /**
   * Gets the block separator for this variable set
   */
  //virtual CFuint getBlockSeparator() const = 0;
   
  /**
   * Set the jacobians
   */
  //virtual void computeJacobians()
  //{
  //  throw Common::NotImplementedException (FromHere(),"MultiFluidMHDVarSet::computeJacobians()");
  //}
  
  /**
   * Compute the pressure derivative
   */
  //virtual void computePressureDerivatives(const Framework::State& state, RealVector& dp)
  //{
  //  throw Common::NotImplementedException (FromHere(),"MultiFluidMHDVarSet::computePressureDerivatives()");
  //}
  
  /**
   * Split the jacobian
   */
  //virtual void splitJacobian(RealMatrix& jacobPlus,
  //			     RealMatrix& jacobMin,
  //			     RealVector& eValues,
  //			     const RealVector& normal) = 0;
  //
  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  //virtual void computeEigenValuesVectors(RealMatrix& rightEv,
  //					 RealMatrix& leftEv,
  //					 RealVector& eValues,
  //					 const RealVector& normal);
  
  


  /// Set the PhysicalData corresponding to the given State
  //virtual void computePhysicalData (const Framework::State& state, RealVector& pdata) = 0;
  
  /// Get the speed
  //virtual CFreal getSpeed(const Framework::State& state) const = 0;
  
  /// Get the normal speed
  //CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const;
  
  /**
   * Get some data corresponding to the subset of equations related with
   * this variable set
   * @pre The most concrete ConvectiveVarSet will have to set these data
   */
  //static std::vector<Framework::EquationSetData>& getEqSetData()    
  //{
  //  static std::vector<Framework::EquationSetData> eqSetData;
  //  return eqSetData;
  //}




 
  /// Set the vector of the eigenValues
  HOST_DEVICE void computeEigenValues (CFreal* pdata, CFreal* normal, CFreal* result)
  {
    //using namespace std;
    //using namespace COOLFluiD::Framework;
    

    
    //const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
    //const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();	// Maxwell's Eqs.+ Multifluid NavierStokes Eqs.
    //const CFuint nbEqs = eqSS.getNbEqsSS();  				
    //const CFuint iEqSS = eqSS.getEqSS();
   
    const CFuint dim = m_dco->dim;
    const bool is2DHalf = m_dco->is2DHalf;


    //if (nbEqs == totalNbEqs || iEqSS == 0) {
      
      const CFreal ce = BASE::getLightSpeed();
      const CFreal chi = BASE::getDivECleaningConst();
      const CFreal gamma = BASE::getDivBCleaningConst();

      result[0] = ce; 
      result[1] = ce;
      result[2] = -ce;
      result[3] = -ce;
      result[4] = gamma*ce;
      result[5] = -gamma*ce;
      result[6] = chi*ce;
      result[7] = -chi*ce;
    //} 
    //if (nbEqs == totalNbEqs || iEqSS == 1) {
      //   EquationSetData& eqSetSpecies = MultiFluidMHDVarSet::getEqSetData()[0];
      // const vector<CFuint>& varIDsSpecies = eqSetSpecies.getEqSetVarIDs();
      const CFuint nbSpecies = m_dco->nbSpecies;             //return 0 en EULERMFMHDTerm  
      // const CFuint firstSpecies = this->getModel()->getFirstScalarVar(0);     
      // EquationSetData& eqSetVelocities = MultiFluidMHDVarSet::getEqSetData()[1];
      const CFuint firstVelocity = m_dco->firstVelocity;      //6 (ConvMaxwellTerm) + 4 (EulerMFMHDTerm)
      const CFuint firstTemperature = m_dco->firstTemperature;   //6 (ConvMaxwellTerm) + 4 (EulerMFMHDTerm)
      const CFuint endEM = 8;
      if (dim == 2 && !is2DHalf){
      
        const CFreal nx = normal[XX];
        const CFreal ny = normal[YY];  
    
        for (CFuint i = 0; i < nbSpecies; ++i) {  
          const CFreal a_i = pdata[firstTemperature + 4*i + 2]; //i species
	  const CFreal u_i = pdata[firstVelocity + dim*i];
	  const CFreal v_i = pdata[firstVelocity + 1 + dim*i];
	  const CFreal un_i = u_i*nx + v_i*ny; 
	
	  result[endEM + i] = un_i;
	  result[endEM + nbSpecies + i] = un_i;
	  result[endEM + 2*nbSpecies + i] = un_i + a_i;
	  result[endEM + 3*nbSpecies + i] = un_i - a_i;    
        } 
      }
      else if (dim == 3 || is2DHalf){

      const CFreal nx = normal[XX];
      const CFreal ny = normal[YY];  
      CFreal nz = 0.;

      (is2DHalf) ? nz = 0 : nz = normal[ZZ];
    
        for (CFuint i = 0; i < nbSpecies; ++i) {  
	  const CFreal a_i = pdata[firstTemperature + 4*i + 2];				//i species
  	  const CFreal u_i = pdata[firstVelocity + dim*i];
	  const CFreal v_i = pdata[firstVelocity + 1+ dim*i];
	  const CFreal w_i = pdata[firstVelocity + dim*i + 2];
	  const CFreal un_i = u_i*nx + v_i*ny + w_i*nz; 
	
	  result[endEM + i] = un_i;
	  result[endEM + nbSpecies + i] = un_i;
	  result[endEM + 2*nbSpecies + i] = un_i;
	  result[endEM + 3*nbSpecies + i] = un_i + a_i;
	  result[endEM + 4*nbSpecies + i] = un_i - a_i;    
        } 
      } 
    //}
  
    //BASE::computeEigenValues(pdata, normals, result); 
  }
  


  /// Get the maximum eigenvalue  NECESITO DECIR SI DE MAXWELL O DEL FLUIDO, O EL MAXIMO DE AMBOS!!!!!
  HOST_DEVICE CFreal getMaxEigenValue(CFreal* pdata, CFreal* normal)
  {
  
      //if (this->getEqSS() == 0) {
      const CFreal ce = m_dco->lightSpeed;
      const CFreal chi = m_dco->divECleaningConst;
      const CFreal gamma = m_dco->divBCleaningConst;
      //printf("ce = %f \t chi = %f \t gamma = %f \n", ce, chi, gamma);
      if (gamma > 1.) {
        if (chi > gamma){
	  return chi*ce;
        }
        else{
	  return gamma*ce;
        }
      }
      else{
        if (chi > 1.){
 	  return chi*ce;
        }
        else {
  	  return ce;
        }
      }
    //}
/* IAA: As the eigenvalues of the Maxwell part are higher, this can be ommited 
    //if (this->getEqSS() == 1) {
      CFreal maxEv;
      maxEv = 0.;
      const CFuint dim = m_dco->dim;  //Fijar o no?
      const CFuint nbSpecies = m_dco->nbSpecies;
      const CFuint firstVelocity = m_dco->firstVelocity;
      const CFuint firstTemperature = m_dco->firstTemperature;
      maxEv = 0;
      for (CFuint i = 0; i < nbSpecies; ++i) {  
        const CFreal a_i = pdata[firstTemperature + 4*i + 2];
        CFreal un_i = 0.0;
        for (CFuint d = 0; d < dim; ++d) {
  	  un_i += pdata[firstVelocity + dim*i+d]*normal[d];
        }
        //maxEv = std::max(0., un_i + a_i);
        if (maxEv < un_i+a_i){maxEv=un_i+a_i;}
      }
      //std::cout << "MultiFluidMHDVarSet<BASE>::getMaxEigenValue" << maxEv << std::endl;
      return maxEv;
    //}
*/   
  
  }





  /// Get the maximum absolute eigenvalue  
  HOST_DEVICE CFreal getMaxAbsEigenValue(CFreal* pdata, CFreal* normal)
  {
    //if (this->getEqSS() == 0) {

      const CFreal ce = m_dco->lightSpeed;
      const CFreal chi = m_dco->divECleaningConst;
      const CFreal gamma = m_dco->divBCleaningConst;
    
      if (gamma > 1){
        if (chi > gamma){
	  return chi*ce;
        }
        else{
	  return gamma*ce;
        }
      }
      else{
        if (chi > 1){
	  return chi*ce;
        }
        else{
 	  return ce;
        }
      }
    //}
  
/* IAA: Same as before
    //if (this->getEqSS() == 1) {
      CFreal maxEv;
      maxEv = 0.;
      const CFuint dim = m_dco->dim;
      const CFuint nbSpecies = m_dco->nbSpecies;
      const CFuint firstVelocity = m_dco->firstVelocity;
      const CFuint firstTemperature = m_dco->firstTemperature;
      for (CFuint i = 0; i < nbSpecies; ++i) {  
        const CFreal a_i = pdata[firstTemperature + 4*i + 2];
        CFreal un_i = 0.0;
        for (CFuint d = 0; d < dim; ++d) {
  	  un_i += pdata[firstVelocity + dim*i+d]*normal[d];
        }
        //maxEv = std::max(maxEv, std::abs(un_i) + a_i);
        if (maxEv < un_i+a_i){maxEv=un_i+a_i;}
      }
      return maxEv;
      //std::cout << "MultiFluidMHDVarSet<BASE>::getMaxAbsEigenValue" << maxEv << std::endl;
    //} 
*/
  }


  /// Computes the convective flux projected on a normal
  HOST_DEVICE void computeFlux(CFreal* pdata, CFreal* normals, CFreal* flux)
  {

    CFuint dim = m_dco->dim;
    const CFuint endEM = 8; 
    const bool is2DHalf = m_dco->is2DHalf;
  
    const CFreal gamma = BASE::getDivBCleaningConst();
    const CFreal ce = BASE::getLightSpeed();
    const CFreal chi = BASE::getDivECleaningConst();


    if (dim == 2 && !is2DHalf){
      const CFreal nx = normals[XX];
      const CFreal ny = normals[YY];
      this->flux[0] = gamma*gamma*pdata[PTERM::PSI]*nx + pdata[PTERM::EZ]*ny;
      this->flux[1] = - pdata[PTERM::EZ]*nx + gamma*gamma*pdata[PTERM::PSI]*ny;
      this->flux[2] = pdata[PTERM::EY]*nx - pdata[PTERM::EX]*ny;
      this->flux[3] = (chi*chi*pdata[PTERM::PHI]*nx- pdata[PTERM::BZ]*ny)*ce*ce;
      this->flux[4] = (pdata[PTERM::BZ]*nx + chi*chi*pdata[PTERM::PHI]*ny)*ce*ce;
      this->flux[5] = (- pdata[PTERM::BY]*nx + pdata[PTERM::BX]*ny)*ce*ce;
      this->flux[6] = (pdata[PTERM::BX]*nx + pdata[PTERM::BY]*ny)*ce*ce;
      this->flux[7] = pdata[PTERM::EX]*nx + pdata[PTERM::EY]*ny;  
    }
    if (dim == 3 || is2DHalf){
      const CFreal nx = normals[XX];
      const CFreal ny = normals[YY];
      CFreal nz = 0;

      (is2DHalf) ? nz = 0 : nz = normals[ZZ];


      this->flux[0] = gamma*gamma*pdata[PTERM::PSI]*nx + pdata[PTERM::EZ]*ny - pdata[PTERM::EY]*nz ;
      this->flux[1] = - pdata[PTERM::EZ]*nx + gamma*gamma*pdata[PTERM::PSI]*ny + pdata[PTERM::EX]*nz;
      this->flux[2] = pdata[PTERM::EY]*nx - pdata[PTERM::EX]*ny + gamma*gamma*pdata[PTERM::PSI]*nz;
      this->flux[3] = chi*chi*ce*ce*pdata[PTERM::PHI]*nx - ce*ce*pdata[PTERM::BZ]*ny +  ce*ce*pdata[PTERM::BY]*nz;
      this->flux[4] = ce*ce*pdata[PTERM::BZ]*nx + chi*chi*ce*ce*pdata[PTERM::PHI]*ny - ce*ce*pdata[PTERM::BX]*nz ;
      this->flux[5] = - ce*ce*pdata[PTERM::BY]*nx + ce*ce*pdata[PTERM::BX]*ny + ce*ce*chi*chi*pdata[PTERM::PHI]*nz;
      this->flux[6] = (pdata[PTERM::BX]*nx + pdata[PTERM::BY]*ny + pdata[PTERM::BZ]*nz)*ce*ce;
      this->flux[7] = pdata[PTERM::EX]*nx + pdata[PTERM::EY]*ny + pdata[PTERM::EZ]*nz;
    }
      
    const CFuint nbSpecies = m_dco->nbSpecies;
    
    const CFuint firstSpecies = m_dco->firstSpecies;
    const CFuint firstVelocity = m_dco->firstVelocity;
    const CFuint firstTemperature = m_dco->firstTemperature;

  

    if (dim == 2 && !is2DHalf){

      const CFreal nx = normals[XX];
      const CFreal ny = normals[YY]; 

      
      for (CFuint i = 0; i < nbSpecies; ++i) {  
        const CFreal rho_i = pdata[PTERM::RHO]*pdata[firstSpecies + i];				//i species
	const CFreal u_i = pdata[firstVelocity + dim*i];
	const CFreal v_i = pdata[firstVelocity + 1 + dim*i];
	const CFreal un_i = u_i*nx + v_i*ny;
	const CFreal rhoVn_i = rho_i*un_i;
	const CFreal P_i = pdata[firstTemperature + 4*i + 1];
	//const CFreal T_i = pdata[firstTemperature + 4*i];      
       	const CFreal H_i = pdata[firstTemperature + 4*i + 3]; 
	
       	this->flux[endEM + i] = rhoVn_i;	
	this->flux[endEM + nbSpecies + i*dim] = rhoVn_i*u_i + P_i*nx;
	this->flux[endEM + nbSpecies + i*dim + 1] = rhoVn_i*v_i + P_i*ny;      
	this->flux[endEM + nbSpecies + nbSpecies*dim + i] = rho_i*un_i*H_i;
      }  
    }
    else if (dim == 3 || is2DHalf){
	
      const CFreal nx = normals[XX];
      const CFreal ny = normals[YY];
      CFreal nz = 0;

      (is2DHalf) ? nz = 0 : nz = normals[ZZ];
	  
      for (CFuint i = 0; i < nbSpecies; ++i) {  
	const CFreal rho_i = pdata[PTERM::RHO]*pdata[firstSpecies + i];				//i species
	const CFreal u_i = pdata[firstVelocity + dim*i];
	const CFreal v_i = pdata[firstVelocity + dim*i + 1];
	const CFreal w_i = pdata[firstVelocity + dim*i + 2];
	const CFreal un_i = u_i*nx + v_i*ny + w_i*nz;
	const CFreal rhoVn_i = rho_i*un_i;      
	const CFreal P_i = pdata[firstTemperature + 4*i + 1]; 
	//const CFreal T_i = pdata[firstTemperature + 4*i];      
	const CFreal H_i = pdata[firstTemperature + 4*i + 3];      
	

        this->flux[endEM + i] = rhoVn_i;
        this->flux[endEM + nbSpecies + i*dim] = rhoVn_i*u_i +  P_i*nx;
        this->flux[endEM + nbSpecies + i*dim + 1] = rhoVn_i*v_i + P_i*ny;
        this->flux[endEM + nbSpecies + i*dim + 2] = rhoVn_i*w_i + P_i*nz;
        this->flux[endEM + nbSpecies + nbSpecies*dim + i] = rho_i*un_i*H_i;
  	//cout <<"Entering Function \n";
      }
    }
  }  




  /// Computes the physical convective flux
  HOST_DEVICE void computeStateFlux(CFreal* pdata)
  {  

    CFuint dim = m_dco->dim;
    const bool is2DHalf = m_dco->is2DHalf;

    //std::cout << "MultiFluidMHDVarSet<BASE>::computeStateFlux" <<"\n";
  
    const CFuint nbSpecies = m_dco->nbSpecies;

    const CFuint firstSpecies = m_dco->firstSpecies;
    const CFuint firstVelocity = m_dco->firstVelocity;
    const CFuint firstTemperature = m_dco->firstTemperature;

    const CFreal ce = BASE::LightSpeed;
    const CFreal gamma = BASE::divBCleaningConst;
    const CFreal chi = BASE::divECleaningConst; 
    const CFuint endEM = 8;


    if (dim == 2 && !is2DHalf){
    
    
      this->_physFlux(0,XX) = gamma*gamma*pdata[PTERM::PSI];
      this->_physFlux(0,YY) = pdata[PTERM::EZ];

      this->_physFlux(1,XX) = - pdata[PTERM::EZ];
      this->_physFlux(1,YY) = gamma*gamma*pdata[PTERM::PSI];

      this->_physFlux(2,XX) = pdata[PTERM::EY];
      this->_physFlux(2,YY) = - pdata[PTERM::EX];

      this->_physFlux(3,XX) = chi*chi*pdata[PTERM::PHI]*ce*ce;
      this->_physFlux(3,YY) = - pdata[PTERM::BZ]*ce*ce;

      this->_physFlux(4,XX) = pdata[PTERM::BZ]*ce*ce;
      this->_physFlux(4,YY) = chi*chi*pdata[PTERM::PHI]*ce*ce;

      this->_physFlux(5,XX) = - pdata[PTERM::BY]*ce*ce;
      this->_physFlux(5,YY) = pdata[PTERM::BX]*ce*ce;
    
      this->_physFlux(6,XX) = pdata[PTERM::BX]*ce*ce;
      this->_physFlux(6,YY) = pdata[PTERM::BY]*ce*ce;

      this->_physFlux(7,XX) = pdata[PTERM::EX];
      this->_physFlux(7,YY) = pdata[PTERM::EY]; 

      for (CFuint i = 0; i < nbSpecies; ++i) {  
        const CFreal u_i = pdata[firstVelocity + dim*i];
	const CFreal v_i = pdata[firstVelocity + dim*i + 1];
	const CFreal rho_i = pdata[PTERM::RHO]*pdata[firstSpecies + i];
	const CFreal rhoU_i = rho_i*u_i;
	const CFreal rhoV_i = rho_i*v_i;
	const CFreal rhoUU_i = rhoU_i*u_i;
	const CFreal rhoVV_i = rhoV_i*v_i;
	const CFreal rhoUV_i = rhoU_i*v_i;
	const CFreal P_i = pdata[firstTemperature + 4*i + 1]; 
	const CFreal H_i = pdata[firstTemperature + 4*i + 3]; 
	
	this->_physFlux(endEM + i,XX) = rhoU_i;
	this->_physFlux(endEM + i,YY) = rhoV_i;
	
	this->_physFlux(endEM + nbSpecies + i*dim,XX) = P_i + rhoUU_i;
	this->_physFlux(endEM + nbSpecies + i*dim,YY) = rhoUV_i;      
	
	this->_physFlux(endEM + nbSpecies + i*dim + 1,XX) = rhoUV_i;
  	this->_physFlux(endEM + nbSpecies + i*dim + 1,YY) = P_i + rhoVV_i;
	
	this->_physFlux(endEM + nbSpecies + nbSpecies*dim + i,XX) = rho_i*u_i*H_i;      
	this->_physFlux(endEM + nbSpecies + nbSpecies*dim + i,YY) = rho_i*v_i*H_i;         
      }
    }
    else if (dim == 3 || is2DHalf){
      this->_physFlux(0,XX) = gamma*gamma*pdata[PTERM::PSI];
      this->_physFlux(0,YY) = pdata[PTERM::EZ];
      this->_physFlux(0,ZZ) = - pdata[PTERM::EY]; 

      this->_physFlux(1,XX) = - pdata[PTERM::EZ];
      this->_physFlux(1,YY) = gamma*gamma*pdata[PTERM::PSI];
      this->_physFlux(1,ZZ) = pdata[PTERM::EX];  

      this->_physFlux(2,XX) = pdata[PTERM::EY];
      this->_physFlux(2,YY) = - pdata[PTERM::EX];
      this->_physFlux(2,ZZ) = gamma*gamma*pdata[PTERM::PSI];  

      this->_physFlux(3,XX) = ce*ce*chi*chi*pdata[PTERM::PHI];
      this->_physFlux(3,YY) = - ce*ce*pdata[PTERM::BZ];
      this->_physFlux(3,ZZ) = ce*ce*pdata[PTERM::BY];  

      this->_physFlux(4,XX) = ce*ce*pdata[PTERM::BZ];
      this->_physFlux(4,YY) = ce*ce*chi*chi*pdata[PTERM::PHI];
      this->_physFlux(4,ZZ) = - ce*ce*pdata[PTERM::BX];  

      this->_physFlux(5,XX) = - ce*ce*pdata[PTERM::BY];
      this->_physFlux(5,YY) = ce*ce*pdata[PTERM::BX];
      this->_physFlux(5,ZZ) = ce*ce*chi*chi*pdata[PTERM::PHI];  
      
      this->_physFlux(6,XX) = pdata[PTERM::BX]*ce*ce;
      this->_physFlux(6,YY) = pdata[PTERM::BY]*ce*ce;
      this->_physFlux(6,ZZ) = pdata[PTERM::BZ]*ce*ce;    

      this->_physFlux(7,XX) = pdata[PTERM::EX];
      this->_physFlux(7,YY) = pdata[PTERM::EY]; 
      this->_physFlux(8,ZZ) = pdata[PTERM::EZ];  

      for (CFuint i = 0; i < nbSpecies; ++i) {  
    
	const CFreal u_i = pdata[firstVelocity + dim*i];
	const CFreal v_i = pdata[firstVelocity + dim*i + 1];
	const CFreal w_i = pdata[firstVelocity + dim*i + 2];      
	const CFreal rho_i = pdata[PTERM::RHO]*pdata[firstSpecies + i];
	const CFreal rhoU_i = rho_i*u_i;
	const CFreal rhoV_i = rho_i*v_i;
	const CFreal rhoW_i = rho_i*w_i;
	const CFreal rhoUU_i = rhoU_i*u_i;
	const CFreal rhoVV_i = rhoV_i*v_i;
	const CFreal rhoWW_i = rhoW_i*w_i;      
	const CFreal rhoUV_i = rhoU_i*v_i;
	const CFreal rhoUW_i = rhoU_i*w_i;
	const CFreal rhoVW_i = rhoV_i*w_i; 
	const CFreal P_i = pdata[firstTemperature + 4*i + 1];
  	const CFreal H_i = pdata[firstTemperature + 4*i + 3];     
	
        this->_physFlux(endEM + i,XX) = rhoU_i;
        this->_physFlux(endEM + i,YY) = rhoV_i;
        this->_physFlux(endEM + i,ZZ) = rhoW_i;
	
        this->_physFlux(endEM + nbSpecies + i*dim,XX) = P_i + rhoUU_i;
        this->_physFlux(endEM + nbSpecies + i*dim,YY) = rhoUV_i;
        this->_physFlux(endEM + nbSpecies + i*dim,ZZ) = rhoUW_i;
  	
        this->_physFlux(endEM + nbSpecies + i*dim + 1,XX) = rhoUV_i;
        this->_physFlux(endEM + nbSpecies + i*dim + 1,YY) = P_i + rhoVV_i;
        this->_physFlux(endEM + nbSpecies + i*dim + 1,ZZ) = rhoVW_i;
	 
        this->_physFlux(endEM + nbSpecies + i*dim + 2,XX) = rhoUW_i;
        this->_physFlux(endEM + nbSpecies + i*dim + 2,YY) = rhoVW_i;
        this->_physFlux(endEM + nbSpecies + i*dim + 2,ZZ) = P_i + rhoWW_i;
	
        this->_physFlux(endEM + nbSpecies + nbSpecies*dim + i,XX) = rho_i*u_i*H_i;
        this->_physFlux(endEM + nbSpecies + nbSpecies*dim + i,YY) = rho_i*v_i*H_i;
        this->_physFlux(endEM + nbSpecies + nbSpecies*dim + i,ZZ) = rho_i*w_i*H_i;
      }

    }      
  }


  HOST_DEVICE CFuint getDim() { return m_dco->dim;}
  HOST_DEVICE bool getIs2DHalf() {return m_dco->is2DHalf;}
  HOST_DEVICE CFuint getNbSpecies() {return m_dco->nbSpecies;}
  HOST_DEVICE CFuint getFirstVelocity() {return m_dco->firstVelocity;}
  HOST_DEVICE CFuint getFirstTemperature() {return m_dco->firstTemperature;}
  HOST_DEVICE CFuint getFirstSpecies() {return m_dco->firstSpecies;}
  HOST_DEVICE CFreal getGamma() {return m_dco->gamma;}
  HOST_DEVICE CFreal getDivECleaningConst() {return m_dco->divECleaningConst;}
  HOST_DEVICE CFreal getDivBCleaningConst() {return m_dco->divBCleaningConst;}
  HOST_DEVICE CFreal getLightSpeed() {return m_dco->lightSpeed;}
  HOST_DEVICE CFreal getPermeability() {return m_dco->mu;}
  HOST_DEVICE CFreal getPermittivity() {return m_dco->epsilon;}
  HOST_DEVICE CFreal getK() {return m_dco->K;}
  HOST_DEVICE CFreal getMolecularMass1() {return m_dco->molecularMass1;}
  HOST_DEVICE CFreal getMolecularMass2() {return m_dco->molecularMass2;}
  HOST_DEVICE CFreal getMolecularMass3() {return m_dco->molecularMass3;}
  HOST_DEVICE CFreal* getNonInducedEMField() {return &m_dco->NonInducedEMField[0];}
  //HOST_DEVICE CFreal getIsLeake() {return m_dco->isLeake;}



protected:
  

  EulerMFMHDTerm::DeviceConfigOptions<NOTYPE>* m_dco;
  CFreal _m_i[3];
 

      
}; // end of class MultiFluidMHD2DHalfTwoSpeciesVarSetT

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD


//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_MultiFluidMHD_MultiFluidMHD2DHalfTwoSpeciesVarSetT_hh
