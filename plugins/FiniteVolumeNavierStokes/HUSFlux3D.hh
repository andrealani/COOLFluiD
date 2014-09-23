#ifndef COOLFluiD_Numerics_FiniteVolume_HUSFlux3D_hh
#define COOLFluiD_Numerics_FiniteVolume_HUSFlux3D_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolumeNavierStokes/HUSFlux2D.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents the HUS flux splitting corresponding
 * to the Euler physical model 3D (in conservative variables)
 *
 * @author Janos Molnar
 * @author Andrea Lani
 * @author Justin Elszasz
 *
 */
template <class UPDATEVAR>   
class HUSFlux3D : public HUSFlux2D<UPDATEVAR> {
public: // classes
  
  /**
   * Constructor
   */
  HUSFlux3D(const std::string& name);

  /**
   * Default destructor
   */
  ~HUSFlux3D();

  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Compute the flux in the current face
   */
  virtual void compute(RealVector& result);
  
private: // helper functions
  
  /**
   * Compute the flux assuming full euler equations + species equations
   * in chemical non equilibrium
   */
  virtual void CNEQ(RealVector& result);
  
  /**
   * Compute the flux assuming thermo-chemical non equilibrium
   */
  virtual void TCNEQ(RealVector& result);
  
private:
  
  /// velocity vectors for left and right states
  RealVector _velL;
  RealVector _velR;
  
  /// vectors for the velocities at the intermediate states
  RealVector _velLS;
  RealVector _velRS;
  
  /// primary tangential direction vector for L and R states
  RealVector _t1L;
  RealVector _t1R;
  
  /// secondary tangential direction vector for L and R states
  RealVector _t2L;
  RealVector _t2R;
  
  /// matrix (and inverses) storing all relative direction information
  RealMatrix _coordsL;
  RealMatrix _coordsR;
  RealMatrix _coordsLinv;
  RealMatrix _coordsRinv;
  
  /// vector storing velocities in relative face coords
  RealVector _velRelL;
  RealVector _velRelR;
  
  /// matrix inverter size 3
  MathTools::MatrixInverterT<3> m_inverter3;
  
}; // end of class HUSFlux3D

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#include "HUSFlux3D.ci"

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_HUSFlux3D_hh
