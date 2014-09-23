#ifndef COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivt_hh
#define COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivt_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/DistanceBasedExtrapolatorGMove.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics {
    
    namespace NavierStokes {
      class EulerTerm;
    }
  }
  
  namespace Framework {
    template <class BASE> class MultiScalarTerm;
    class PhysicalChemicalLibrary;
  }
  
  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object to be used 
 * in combination with BCs that moves the ghost states (nodes) like
 * @see NoSlipWallIsothermalNSPvt
 *
 * @author Andrea Lani
 *
 */
class DistanceBasedExtrapolatorGMoveRhoivt : public DistanceBasedExtrapolatorGMove {
public:

  /**
   * Constructor
   */
  DistanceBasedExtrapolatorGMoveRhoivt(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DistanceBasedExtrapolatorGMoveRhoivt();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);


  /**
   * Set up private data needed by the computation
   */
  virtual void setup();

protected:
    
  /**
   * Apply the nodal extrapolation
   */
  virtual void extrapolate();
  
  /**
   * Transform a given state in another one before accumulating the values
   */
  virtual void transform(const RealVector& in, RealVector& out);
  
  /**
   * Transform back a given state in the original one after accumulating the values
   */
  virtual void transformBack(RealVector& nstate);
  
  /**
   * Apply the boundary condition
   */
  virtual void applyBC();
  
  /// Set the species-related variables in a given state
  void setSpeciesVariables(RealVector& nstate);
  
protected:
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// acquaintance of the physical model term
  Common::SafePtr<Framework::MultiScalarTerm<Physics::NavierStokes::EulerTerm> > _model;
  
  /// temporary state vector
  RealVector _tmpState;

  /// temporary mass fractions
  RealVector _ys;
  
  /// species molar masses array
  RealVector _mmasses;
  
  /// nodal pressure
  CFreal _nodalPressure;
  
  /// temporary pressure value
  CFreal _pressure;
  
}; // end of class DistanceBasedExtrapolatorGMoveRhoivt

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_DistanceBasedExtrapolatorGMoveRhoivt_hh
