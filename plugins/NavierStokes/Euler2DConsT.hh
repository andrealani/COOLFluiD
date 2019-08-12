#ifndef COOLFluiD_Physics_NavierStokes_Euler2DConsT_hh
#define COOLFluiD_Physics_NavierStokes_Euler2DConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Euler2DVarSetT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative variables
 *
 * @author Ray Vandenhoeck
 */
class Euler2DConsT : public Euler2DVarSetT {
public: // function
  
  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  HOST_DEVICE Euler2DConsT(EulerTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    Euler2DVarSetT(dco) {}
 
  /**
   * Constructor
   * @see EulerPhysicalModel
   */
  HOST_DEVICE Euler2DConsT() : Euler2DVarSetT() {}
  
  /**
   * Default destructor
   */
  HOST_DEVICE ~Euler2DConsT() {}
  
  /// Compute the physical data starting ffrom the corresponding state variables
  HOST_DEVICE void computePhysicalData(CFreal* state, CFreal* data) 
  { 
    // we assume that if conservative variables are used, the flow is compressible
    // p = static pressure
    const CFreal rho  = state[0];
    const CFreal ovRho = 1./rho;
    const CFreal u = state[1]*ovRho;
    const CFreal v = state[2]*ovRho;
    const CFreal V2 = u*u + v*v;
    const CFreal gamma = m_dco->gamma;
    const CFreal gammaMinus1 = gamma - 1.;
    const CFreal rhoE = state[3];
  
    data[EulerTerm::RHO] = rho;
    data[EulerTerm::P] = gammaMinus1*(rhoE - 0.5*rho*V2);
  
    const CFreal pOvRho = data[EulerTerm::P]*ovRho;
    data[EulerTerm::E] = rhoE*ovRho;
    data[EulerTerm::H] = data[EulerTerm::E] + pOvRho;
    data[EulerTerm::A] = sqrt(gamma*pOvRho);
    data[EulerTerm::T] = pOvRho/m_dco->R;
    data[EulerTerm::V] = sqrt(V2);
    data[EulerTerm::VX] = u;
    data[EulerTerm::VY] = v;
    data[EulerTerm::GAMMA] = gamma;
  }
  
}; // end of class Euler2DConsT

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_Euler2DConsT_hh

