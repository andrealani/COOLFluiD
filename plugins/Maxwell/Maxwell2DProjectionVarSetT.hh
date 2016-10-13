#ifndef COOLFluiD_Physics_Maxwell_Maxwell2DProjectionVarSetT_hh
#define COOLFluiD_Physics_Maxwell_Maxwell2DProjectionVarSetT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Maxwell/MaxwellProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a GPU-enabled variable set for projection scheme
 * for 2D Maxwell physical model.
 *
 * @author Alejandro Alvarez Laguna
 * @author Andrea Lani
 * @author Isaac Alonso
 */
class Maxwell2DProjectionVarSetT  {

public: // classes

  enum {DIM=2, NBEQS=8, DATASIZE = 8}; 
  typedef MaxwellProjectionTerm PTERM;
  typedef Maxwell2DProjectionVarSetT EULERSET;
  
  /**
   * Constructor
   * @see MaxwellModel
   */
  HOST_DEVICE Maxwell2DProjectionVarSetT(MaxwellProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) :
    m_dco(dco) {} 
  

  HOST_DEVICE Maxwell2DProjectionVarSetT() {}



  HOST_DEVICE void setModelData(MaxwellProjectionTerm::DeviceConfigOptions<NOTYPE>* dco) {m_dco = dco;}


  /**
   * Default destructor
   */
  HOST_DEVICE virtual ~Maxwell2DProjectionVarSetT() {};
  
  
  /// Set the vector of the eigenValues 

  HOST_DEVICE void computeEigenValues(CFreal* data,
				       CFreal* normal,
				       CFreal* result)
  {

    const CFreal ce = m_dco->LightSpeed;
    const CFreal chi = m_dco->divECleaningConst;
    const CFreal gamma = m_dco->divBCleaningConst;

    result[0] = ce; 
    result[1] = ce;
    result[2] = -ce;
    result[3] = -ce;
    result[4] = gamma*ce;
    result[5] = -gamma*ce;
    result[6] = chi*ce;
    result[7] = -chi*ce;
  
  }



  /// Get the maximum eigenvalue  
  HOST_DEVICE CFreal getMaxEigenValue(CFreal* data, CFreal* normal)
  {

    const CFreal ce = m_dco->LightSpeed;
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

  }



  /// Get the maximum absolute eigenvalue
  HOST_DEVICE CFreal getMaxAbsEigenValue(CFreal* data, CFreal* normal)
  {

    const CFreal ce = m_dco->LightSpeed;
    const CFreal chi = m_dco->divECleaningConst;
    const CFreal gamma = m_dco->divBCleaningConst;
  
    if (gamma > 1.){
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
      else{
        return ce;
      }
    }

  }
  


  /// Computes the convective flux projected on a normal
  HOST_DEVICE void computeFlux(CFreal* data, CFreal* normals, CFreal* flux)
  {    

    //CFLog(NOTICE, "computeFlux()\n");

    const CFreal nx = normals[XX];
    const CFreal ny = normals[YY];
    const CFreal gamma = m_dco->divBCleaningConst;
    const CFreal ce = m_dco->LightSpeed;
    const CFreal chi = m_dco->divECleaningConst;

    flux[0] = gamma*gamma*data[MaxwellProjectionTerm::PSI]*nx + data[ConvMaxwellTerm::EZ]*ny;
    flux[1] = - data[ConvMaxwellTerm::EZ]*nx + gamma*gamma*data[MaxwellProjectionTerm::PSI]*ny;
    flux[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny;
    flux[3] = (chi*chi*data[MaxwellProjectionTerm::PHI]*nx- data[ConvMaxwellTerm::BZ]*ny)*ce*ce;
    flux[4] = (data[ConvMaxwellTerm::BZ]*nx + chi*chi*data[MaxwellProjectionTerm::PHI]*ny)*ce*ce;
    flux[5] = (- data[ConvMaxwellTerm::BY]*nx + data[ConvMaxwellTerm::BX]*ny)*ce*ce;
    flux[6] = (data[ConvMaxwellTerm::BX]*nx + data[ConvMaxwellTerm::BY]*ny)*ce*ce;
    flux[7] = data[ConvMaxwellTerm::EX]*nx + data[ConvMaxwellTerm::EY]*ny;
   
  }


  
  HOST_DEVICE CFreal getDivBCleaningConst() { return m_dco->divBCleaningConst;}
  HOST_DEVICE CFreal getDivECleaningConst() {return m_dco->divECleaningConst;}
  HOST_DEVICE CFreal getLightSpeed() {return m_dco->LightSpeed;}


   /// @returns the number of equations of this VarSet
  HOST_DEVICE CFuint getNbEqs() const { return 8;}

protected:

  MaxwellProjectionTerm::DeviceConfigOptions<NOTYPE>* m_dco;

}; // end of class Maxwell2DProjectionVarSetT

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Maxwell_Maxwell2DProjectionVarSetT_hh
