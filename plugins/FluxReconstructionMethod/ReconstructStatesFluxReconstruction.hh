#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_ReconstructStatesFluxReconstruction_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_ReconstructStatesFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command that reconstructs the states in a given set of points
 *
 * @author Kris Van den Abeele
 */
class ReconstructStatesFluxReconstruction : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      FluxReconstructionSolverData,ReconstructStatesFluxReconstruction > PROVIDER;

public:  // methods

  /// Constructor
  ReconstructStatesFluxReconstruction(const std::string& name);

  /// Destructor
  ~ReconstructStatesFluxReconstruction();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ReconstructStatesFluxReconstruction";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return "ReconstructStatesFluxReconstruction";}
  
  /// Set up private data and data
  virtual void setup();

  /// reconstruct one state
  void reconstructState(const std::vector< Framework::State* >& cellStates,
                        Framework::State& recState,
                        const RealVector& recCoefs,
                        const std::vector< CFuint >& cellStateIdxsForRec);

  /// reconstruct the states in the given points
  void reconstructStates(const std::vector< Framework::State* >& cellStates,
                         std::vector< Framework::State* >& recStates,
                         const RealMatrix& recCoefs,
                         const std::vector< CFuint >& recStateIdxs,
                         const std::vector< CFuint >& recStateMatrixIdx,
                         const std::vector< std::vector< CFuint > >& cellStateIdxs);

  /// reconstruct the extra variables in the given points
  void reconstructExtraVars(const std::vector< RealVector* >& cellExtraVars,
                            std::vector< RealVector* >& recExtraVars,
                            const RealMatrix& recCoefs,
                            const std::vector< CFuint >& recVarIdxs,
                            const std::vector< CFuint >& recVarMatrixIdx,
                            const std::vector< std::vector< CFuint > >& cellVarIdxs);

  /// reconstruct the gradients in the given points
  void reconstructGradients(const std::vector< std::vector< RealVector >* >& cellGradients,
                            std::vector< std::vector< RealVector* > >& recGradients,
                            const RealMatrix& recCoefs,
                            const std::vector< CFuint >& recGradIdxs,
                            const std::vector< CFuint >& recGradMatrixIdxs,
                            const std::vector< std::vector< CFuint > >& cellGradIdxs);

  /// reconstruct one physical variable in the given points
  void reconstructPhysVar(const CFuint iVar,
                          const std::vector< Framework::State* >& cellStates,
                          std::vector< Framework::State* >& recStates,
                          const RealMatrix& recCoefs,
                          const std::vector< CFuint >& recStateIdxs,
                          const std::vector< CFuint >& recStateMatrixIdxs,
                          const std::vector< std::vector< CFuint > >& cellStateIdxs);

  /// reconstruct the gradients in the given points
  void reconstructPhysVarGrad(const CFuint iVar,
                              const std::vector< std::vector< RealVector >* >& cellGradients,
                              std::vector< std::vector< RealVector* > >& recGradients,
                              const RealMatrix& recCoefs,
                              const std::vector< CFuint >& recGradIdxs,
                              const std::vector< CFuint >& recGradMatrixIdxs,
                              const std::vector< std::vector< CFuint > >& cellGradIdxs);

  /// compute solution polynomial gradients in the given points
  void computePolyGradients(const std::vector< Framework::State* >& cellStates,
                            std::vector< std::vector< RealVector* > >& recGradients,
                            const std::vector< std::vector< std::vector< CFreal > > >& derivCoefs,
                            const std::vector< CFuint >& recGradIdxs,
                            const std::vector< std::vector< std::vector< CFuint > > >& cellStateIdxs,
                            const std::vector< RealMatrix >& invJacobMatr);

  /// compute solution polynomial gradients in the given points
  void computePolyGradient(const CFuint iVar,
                           const std::vector< Framework::State* >& cellStates,
                           std::vector< std::vector< RealVector* > >& recGradients,
                           const std::vector< std::vector< std::vector< CFreal > > >& derivCoefs,
                           const std::vector< CFuint >& recGradIdxs,
                           const std::vector< std::vector< std::vector< CFuint > > >& cellStateIdxs,
                           const std::vector< RealMatrix >& invJacobMatr);

private: // data

/// all state gradients with respect to the mapped coordinates
std::vector< RealVector > m_allGradMappedCoord;

/// one gradient with respect to the mapped coordinates
RealVector m_gradMappedCoord;

}; // class ReconstructStatesFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_ReconstructStatesFluxReconstruction_hh

