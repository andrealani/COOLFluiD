#ifndef COOLFluiD_Numerics_FiniteVolume_RoeTCNEQFlux_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeTCNEQFlux_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux for thermo-chemical non-equilibrium applications
 *
 * @author Andrea Lani
 *
 */
template <class BASE , class UPDATEVAR>
class RoeTCNEQFlux : public BASE {
public:

  /**
   * Constructor
   */
  RoeTCNEQFlux(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeTCNEQFlux();

  /**
   * Set up private data
   */
  virtual void setup();
    
private:
  
  /**
   * Perform ad-hoc linearization
   */
  void linearize();
  
private:
  
  /// acquaintance of the concrete variable set
  Common::SafePtr<UPDATEVAR> _upVar;
  
  /// physical model data
  RealVector _lData;

  /// physical model data
  RealVector _rData;
  
  /// species molar masses
  RealVector _mmasses;
   
  /// f_i coefficients
  RealVector _fcoeff;
  
  /// total variation of partial densities over corresponding molar masses 
  RealVector _dRhoiOvMM;
        
}; // end of class RoeTCNEQFlux

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "RoeTCNEQFlux.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeTCNEQFlux_hh
