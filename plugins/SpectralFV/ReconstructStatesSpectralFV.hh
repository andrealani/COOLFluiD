#ifndef COOLFluiD_Numerics_SpectralFV_ReconstructStatesSpectralFV_hh
#define COOLFluiD_Numerics_SpectralFV_ReconstructStatesSpectralFV_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that reconstructs the states in a given set of points
 *
 * @author Kris Van den Abeele
 */
class ReconstructStatesSpectralFV : public SpectralFVMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFVMethodData,ReconstructStatesSpectralFV > PROVIDER;

public:  // methods

  /// Constructor
  ReconstructStatesSpectralFV(const std::string& name);

  /// Destructor
  ~ReconstructStatesSpectralFV();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ReconstructStatesSpectralFV";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// Set up private data and data
  virtual void setup() {}

  /// reconstruct one state
  void reconstructState(const std::vector< Framework::State* >& cellStates,
                        Framework::State& recState,
                        const std::vector< CFreal >& recCoefs,
                        const CFuint nbrCellStates);

  /// reconstruct the states in the given points
  void reconstructStates(const std::vector< Framework::State* >& cellStates,
                         std::vector< Framework::State* >& recStates,
                         const std::vector< std::vector< CFreal > >& recCoefs,
                         const CFuint nbrCellStates);

  /// reconstruct the extra variables in the given points
  void reconstructExtraVars(const std::vector< RealVector* >& cellExtraVars,
                            std::vector< RealVector* >& recExtraVars,
                            const std::vector< std::vector< CFreal > >& recCoefs,
                            const CFuint nbrCellExtraVars);

  /// reconstruct the gradients in the given points
  void reconstructGradients(const std::vector< std::vector< RealVector >* >& cellGradients,
                            std::vector< std::vector< RealVector* > >& recGradients,
                            const std::vector< std::vector< CFreal > >& recCoefs,
                            const CFuint nbrCellGrads);

  /// reconstruct one physical variable in the given points
  void reconstructPhysVar(const CFuint iVar,
                          const std::vector< Framework::State* >& cellStates,
                          std::vector< Framework::State* >& recStates,
                          const std::vector< std::vector< CFreal > >& recCoefs,
                          const CFuint nbrCellStates);

}; // class ReconstructStatesSpectralFV

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_ReconstructStatesSpectralFV_hh

