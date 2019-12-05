#ifndef COOLFluiD_Physics_NavierStokes_NavierStokes2DConsT_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokes2DConsT_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DVarSetT.hh"
#include "Euler2DConsT.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NavierStokes physical model 2D for conservative variables
 *
 * @author Ray Vandenhoeck
 */
class NavierStokes2DConsT : public NavierStokes2DVarSetT {
public: // function
  
  /**
   * Constructor
   * @see NavierStokesPhysicalModel
   */
  HOST_DEVICE NavierStokes2DConsT(NSTerm::DeviceConfigOptions<NOTYPE>* dcoNS, EulerTerm::DeviceConfigOptions<NOTYPE>* dco) : 
    NavierStokes2DVarSetT(dcoNS, dco) {}
 
  /**
   * Constructor
   * @see NavierStokesPhysicalModel
   */
  HOST_DEVICE NavierStokes2DConsT() : NavierStokes2DVarSetT() {}
  
  /**
   * Default destructor
   */
  HOST_DEVICE ~NavierStokes2DConsT() {}
  
  /// Computes the convective flux projected on a normal
  HOST_DEVICE void getFlux(CFreal* state, CFreal* gradients, CFreal* normal, CFreal* flux) 
  {   
    // velocity
    const CFreal u = state[1]/state[0];
    const CFreal v = state[2]/state[0];
    
    // dynamic viscosity
    const CFreal mu = getDynViscosity(state);
    
    // thermal conductivity
    const CFreal lambda = getThermConductivity(state, mu);
    
    // velocity gradients
    const CFreal gradUX = gradients[1*2];
    const CFreal gradUY = gradients[1*2+1];
    const CFreal gradVX = gradients[2*2];
    const CFreal gradVY = gradients[2*2+1];

    // two thirds div u
    const CFreal twoThirdDivV = 2.0/3.0*(gradUX + gradVY);
    
    // stress tensor
    const CFreal coeffTauMu = m_dcoNS->coeffTau*mu;
    const CFreal tauXX = coeffTauMu*(2.*gradUX - twoThirdDivV);
    const CFreal tauXY = coeffTauMu*(gradUY + gradVX);
    const CFreal tauYY = coeffTauMu*(2.*gradVY - twoThirdDivV);
  
    // heat transfer constribution
    const CFreal qFlux = -m_dcoNS->coeffQ*lambda*(gradients[3*2]*normal[XX] + gradients[3*2+1]*normal[YY]);

    // flux
    flux[0] = 0.0;
    flux[1] = tauXX*normal[XX] + tauXY*normal[YY];
    flux[2] = tauXY*normal[XX] + tauYY*normal[YY];
    flux[3] = (tauXX*u + tauXY*v)*normal[XX] + (tauXY*u + tauYY*v)*normal[YY] - qFlux;
  }
  
/////////////////////////////////////////////////////////////////////////////////
  
  HOST_DEVICE void setGradientVars(const CFreal* state, CFreal* values)
{  
  const CFreal R = m_dco->R;
  const CFreal ovCv = (m_dco->gamma - 1.)/R;
  
  const CFreal ovRho = 1./state[0]; 
  values[1] = state[1]*ovRho;
  values[2] = state[2]*ovRho;
    
  const CFreal V2 = values[1]*values[1] + values[2]*values[2];
  values[3] = (state[3]*ovRho - 0.5*V2)*ovCv;
  values[0] = R*state[0]*values[3];
}
  
  //////////////////////////////////////////////////////////////////////////////

HOST_DEVICE CFreal getDynViscosity(const CFreal* state)
{
  // here it is assumed that state is in Cons variables
  const CFreal R = m_dco->R;
  const CFreal cv = R/(m_dco->gamma - 1.);
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal T = (state[3] - 0.5*state[0]*(u*u + v*v))/(state[0]*cv);

  return 0.000001458*pow(T,1.5)/(110.4+T);
  //        /(getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

HOST_DEVICE CFreal getDensity(const CFreal* state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

HOST_DEVICE CFreal getThermConductivity(const CFreal* state, const CFreal& dynViscosity)
{
  //if (Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional()) {
  //  return dynViscosity;
  //}
  return dynViscosity*m_dcoNS->cpOverPrandtl;
}

  
}; // end of class NavierStokes2DConsT

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokes2DConsT_hh

