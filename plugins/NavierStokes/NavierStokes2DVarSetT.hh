#ifndef COOLFluiD_Physics_NavierStokes_NavierStokes2DVarSetT_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokes2DVarSetT_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NSTerm.hh"
#include "NavierStokes/EulerTerm.hh"


using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for the 2D NS physical model ported to GPU.
 *
 * @author Ray Vandenhoeck
 */
class NavierStokes2DVarSetT {
  
public: // classes

  enum {DIM=2, NBEQS=4, DATASIZE = 14};
  typedef EulerTerm PTERM;
  typedef NSTerm DTERM;
  
  /**
   * Constructor
   * @see NavierStokesPhysicalModel
   */
  HOST_DEVICE NavierStokes2DVarSetT(NSTerm::DeviceConfigOptions<NOTYPE>* dcoNS, EulerTerm::DeviceConfigOptions<NOTYPE>* dco) :
    m_dcoNS(dcoNS),
    m_dco(dco) {} 
  
  /**
   * Constructor
   * @see NavierStokesPhysicalModel
   */
  HOST_DEVICE NavierStokes2DVarSetT() {}
  
  /**
   * Set the model data
   */
  HOST_DEVICE void setModelData(NSTerm::DeviceConfigOptions<NOTYPE>* dcoNS, EulerTerm::DeviceConfigOptions<NOTYPE>* dco) 
  {
      m_dcoNS = dcoNS;
      m_dco = dco;
  }
  
  /**
   * Default destructor
   */
  HOST_DEVICE virtual ~NavierStokes2DVarSetT() {}
  
protected:
  
  /// configurable options
  NSTerm::DeviceConfigOptions<NOTYPE>* m_dcoNS;
  
  EulerTerm::DeviceConfigOptions<NOTYPE>* m_dco;
  
}; // end of class NavierStokes2DVarSetT

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokes2DVarSetT_hh

