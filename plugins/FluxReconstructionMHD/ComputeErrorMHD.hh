#ifndef COOLFluiD_FluxReconstructionMethod_ComputeErrorMHD_hh
#define COOLFluiD_FluxReconstructionMethod_ComputeErrorMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/TensorProductGaussIntegrator.hh"
#include "FluxReconstructionMethod/CellToFaceGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {
    class FluxReconstructionElementData;
  }

  namespace Physics {
    namespace MHD {
      class MHD2DProjectionVarSet;
    }
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class is used to calculate the error of the solution of FR for MHD
 *
 * @author Rayan Dhib
 */
class ComputeErrorMHD : public FluxReconstructionSolverCom {

public:  // types

  /// Enum for reference solution types
  enum RefSolutionType {
    MANUFACTURED = 0,
    ALFVEN_WAVE = 1,
    MHD_VORTEX = 2
  };

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ComputeErrorMHD(const std::string& name);

  /// Destructor
  ~ComputeErrorMHD();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeErrorMHD";
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
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> m_MHDVarSet;
  
  /// variable for physical data
  RealVector m_solPhysData;
  
  /// coefficients for the computation of the cell averaged solution
  Common::SafePtr< RealVector > m_cellAvgSolCoefs;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

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

  /// local coordinates of the quadrature solution points
  Common::SafePtr< std::vector< RealVector > > m_quadSolPntsLocalCoords;

  /// cell averaging coefficients for quadrature points
  Common::SafePtr< RealVector > m_quadCellAvgSolCoefs;

  /// polynomial values at quadrature points for extrapolation
  std::vector< std::vector< CFreal > > m_solPolyValsAtQuadPnts;

  /// number of quadrature solution points
  CFuint m_nbrQuadSolPnts;

  /// number of equations
  CFuint m_nbrEqs;

  /// flux reconstruction element data for quadrature rule
  FluxReconstructionElementData* m_frLocalDataQuad;

  /// coefficients for reconstructing the state in the quadrature points
  std::vector< std::vector< CFreal > > m_quadCoefs;

  /// total L2 error for all processors
  CFreal m_total_error_rho;

  /// total Li error for all processors
  CFreal m_total_error_rho_Li;

  /// total avg error for all processors
  CFreal m_total_avg_e_rho;

  /// total L1 error for all processors
  CFreal m_total_error_rho_L1;

  /// total magnetic field error for all processors
  CFreal m_total_error_B;

  /// total nb of states for all processors
  CFuint m_total_nbStates;

  /// L2 error 
  CFreal m_error_rho;

  /// Li error
  CFreal m_error_rho_Li;

  /// L1 error (replaces avg error for MHD vortex case)
  CFreal m_error_rho_L1;

  /// Magnetic field error for MHD vortex case
  CFreal m_error_B;

  /// Fixed quadrature order for error integration
  CFuint m_quadOrder;

  /// Reference solution type
  RefSolutionType m_refSolutionType;

  /// Reference solution type string (config option)
  std::string m_refSolutionTypeStr;

  /// Show rate of the error information
  CFuint m_showRate;

  /// ID of the variable to monitor for error computation
  CFuint m_monitoredVar;


protected: // helper functions

  /**
   * Compute analytical solution based on the selected reference solution type
   * @param x,y,z physical coordinates
   * @param analyticalSol output vector with analytical solution values
   */
  void computeAnalyticalSolution(const CFreal x, const CFreal y, const CFreal z, 
                                 RealVector& analyticalSol);

  /**
   * Compute manufactured solution
   * @param x,y,z physical coordinates
   * @param analyticalSol output vector with analytical solution values
   */
  void computeManufacturedSolution(const CFreal x, const CFreal y, const CFreal z,
                                   RealVector& analyticalSol);

  /**
   * Compute Alfven wave solution (placeholder)
   * @param x,y,z physical coordinates
   * @param analyticalSol output vector with analytical solution values
   */
  void computeAlfvenWaveSolution(const CFreal x, const CFreal y, const CFreal z,
                                 RealVector& analyticalSol);

  /**
   * Compute MHD vortex solution
   * @param x,y,z physical coordinates
   * @param analyticalSol output vector with analytical solution values
   */
  void computeMHDVortexSolution(const CFreal x, const CFreal y, const CFreal z,
                                RealVector& analyticalSol);

protected: // data
}; // class ComputeErrorMHD

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_ComputeErrorMHD_hh

