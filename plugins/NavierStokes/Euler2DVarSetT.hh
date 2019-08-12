#ifndef COOLFluiD_Physics_NavierStokes_Euler2DVarSetT_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DVarSetT_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/EulerTerm.hh"

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for the 2D Euler physical model ported to GPU.
 *
 * @author Ray Vandenhoeck
 */
class Euler2DVarSetT {
  
public: // classes

  enum {DIM=2, NBEQS=4, DATASIZE = 14};
  typedef EulerTerm PTERM;
  
  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  HOST_DEVICE Euler2DVarSetT(EulerTerm::DeviceConfigOptions<NOTYPE>* dco) :
    m_dco(dco) {} 
  
  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  HOST_DEVICE Euler2DVarSetT() {}
  
  /**
   * Set the model data
   */
  HOST_DEVICE void setModelData(EulerTerm::DeviceConfigOptions<NOTYPE>* dco) {m_dco = dco;}
  
  /**
   * Default destructor
   */
  HOST_DEVICE virtual ~Euler2DVarSetT() {}
  
  /// Computes the convective flux projected on a normal
  HOST_DEVICE void getFlux(CFreal* data, CFreal* normal, CFreal* flux) 
  {
    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];
    const CFreal rho = data[EulerTerm::RHO];
    const CFreal u = data[EulerTerm::VX];
    const CFreal v = data[EulerTerm::VY];
    const CFreal un = u*nx  + v*ny;
    const CFreal rhoVn = data[EulerTerm::RHO]*un;
    const CFreal p = data[EulerTerm::P];

    flux[0] = rhoVn;
    flux[1] = p*nx + u*rhoVn;
    flux[2] = p*ny + v*rhoVn;
    flux[3] = rhoVn*data[EulerTerm::H];
  }
  
  
  /// Set the vector of the eigenValues
  HOST_DEVICE  void computeEigenValues (CFreal* data, CFreal* normal, CFreal* eValues)
  {
    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];
    const CFreal u = data[EulerTerm::VX];
    const CFreal v = data[EulerTerm::VY];
    const CFreal un = u*nx  + v*ny;
    const CFreal a = data[EulerTerm::A];

    eValues[0] = un;
    eValues[1] = un;
    eValues[2] = un + a;
    eValues[3] = un - a;
  }
  
  /// Compute the maximum eigenvalue
  HOST_DEVICE  CFreal getMaxEigenValue(CFreal* data, CFreal* normal)
  {
    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];
    const CFreal u = data[EulerTerm::VX];
    const CFreal v = data[EulerTerm::VY];
    const CFreal un = u*nx  + v*ny;
    return un + data[EulerTerm::A];
  }
  
  /// Compute the maximum absolute value eigenvalue
  HOST_DEVICE  CFreal getMaxAbsEigenValue(CFreal* data, CFreal* normal)
  {
    const CFreal nx = normal[XX];
    const CFreal ny = normal[YY];
    const CFreal u = data[EulerTerm::VX];
    const CFreal v = data[EulerTerm::VY];
    const CFreal un = u*nx  + v*ny;
    return fabs(un) + data[EulerTerm::A];
  }
  
protected:
  
  /// configurable options
  EulerTerm::DeviceConfigOptions<NOTYPE>* m_dco;
  
}; // end of class Euler2DVarSetT

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DVarSetT_hh

