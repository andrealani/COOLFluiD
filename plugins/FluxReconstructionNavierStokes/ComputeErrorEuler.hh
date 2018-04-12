#ifndef COOLFluiD_FluxReconstructionMethod_ComputeErrorEuler_hh
#define COOLFluiD_FluxReconstructionMethod_ComputeErrorEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"
#include "FluxReconstructionMethod/CellToFaceGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;//3D
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is used to calculate the error of the solution of FR for Euler
 *
 * @author Ray Vandenhoeck
 */
class ComputeErrorEuler : public FluxReconstructionSolverCom {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ComputeErrorEuler(const std::string& name);

  /// Destructor
  ~ComputeErrorEuler();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeErrorEuler";
  }
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

  /// Set up private data and data
  void setup();
  

  /**
   * Computes the error
   */
  void execute();


protected: // data
  
  /// socket for state's
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;//3D
  
  /// variable for physical data
  RealVector m_solPhysData;
  
  /// coefficients for the computation of the cell averaged solution
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;
  
  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;
  
  /// variable for cell
  Framework::GeometricEntity* m_cell;
  
  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;
  
  /// gaus integrator
  TensorProductGaussIntegrator m_tpIntegrator;
  
  /// mapped coords of the gauss quadrature points
  std::vector< RealVector > m_quadPntCoords;
  
  /// coefficients for reconstructing the state in the quadrature points
  std::vector< std::vector< CFreal > > m_quadCoefs;


}; // class ComputeErrorEuler

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_ComputeErrorEuler_hh

