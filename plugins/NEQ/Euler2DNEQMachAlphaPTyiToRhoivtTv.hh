#ifndef COOLFluiD_Physics_NEQ_Euler2DNEQMachAlphaPTyiToRhoivtTv_hh
#define COOLFluiD_Physics_NEQ_Euler2DNEQMachAlphaPTyiToRhoivtTv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VarSetTransformer.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a transformer of variables from MachAlphaPTyi
 * to conservative variables
 *
 * @author Andrea Lani
 *
 */
class Euler2DNEQMachAlphaPTyiToRhoivtTv : public Framework::VarSetTransformer {
public:
  
  typedef Framework::MultiScalarTerm<NavierStokes::EulerTerm> NEQTerm;
  
  /**
   * Default constructor without arguments
   */
  Euler2DNEQMachAlphaPTyiToRhoivtTv(Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Euler2DNEQMachAlphaPTyiToRhoivtTv();

  /**
   * Set up
   */
  void setup(const CFuint maxNbTransStates);
  
  /**
   * Transform a state into another one
   */
  void transform(const Framework::State& state, Framework::State& result);
  
  /**
   * Transform a state into another one from reference precomputed
   * values (physical data)associated to the given state
   */
  void transformFromRef(const RealVector& data, Framework::State& result);
  
private:
  
  /// acquaintance of the model
  Common::SafePtr<NEQTerm> _model;
    
  /// array with all different vibrational dimensional temperatures
  RealVector _tvDim;
    
  /// species molar masses
  RealVector _mmasses;
  
  /// species mass fractions
  RealVector _ys;
  
}; // end of class Euler2DNEQMachAlphaPTyiToRhoivtTv

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NEQ_Euler2DNEQMachAlphaPTyiToRhoivtTv_hh
