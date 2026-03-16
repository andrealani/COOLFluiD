#ifndef COOLFluiD_FluxReconstructionMethod_OutputGradients_hh
#define COOLFluiD_FluxReconstructionMethod_OutputGradients_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This command copies the FR gradient socket (vector<RealVector> per state)
 * into a flat CFreal socket that can be picked up by DataHandleOutput
 * for visualization in CGNS, Tecplot, or ParaView writers.
 *
 * Configure as:
 *   FluxReconstruction.ComputeErrorCom = OutputGradients
 *   FluxReconstruction.OutputGradients.OutputVarIDs = 0 1 2  (optional, default: all)
 *   FluxReconstruction.OutputGradients.OutputMagnitude = true (default: true)
 *
 * OutputMagnitude = true:  one scalar |grad(var)| per variable (for AMR indicator selection)
 * OutputMagnitude = false: dim components per variable (dx, dy, dz)
 *
 * Then connect the output socket to the writer:
 *   CGNS.Data.DataHandleOutput.SocketNames = gradOutput
 *   CGNS.Data.DataHandleOutput.VariableNames = grad_rho grad_rhoU grad_rhoV ...  (magnitude mode)
 *   CGNS.Data.DataHandleOutput.VariableNames = drho_dx drho_dy drho_dz ...       (component mode)
 *
 * @author Rayan Dhib
 */
class OutputGradients : public FluxReconstructionSolverCom {

public:

  static void defineConfigOptions(Config::OptionList& options);

  OutputGradients(const std::string& name);

  ~OutputGradients();

  static std::string getClassName() { return "OutputGradients"; }

  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

  std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  void setup();

  void unsetup();

  void execute();

private:

  /// sink: gradient data (vector<RealVector> per state, sized nEqs, each RealVector dim-sized)
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;

  /// sink: states (needed for state count)
  Framework::DataSocketSink< Framework::State*, Framework::GLOBAL > socket_states;

  /// source: flat CFreal buffer for DataHandleOutput
  Framework::DataSocketSource< CFreal > socket_gradOutput;

  /// user-selected variable indices (empty = all)
  std::vector< CFuint > m_outputVarIDs;

  /// resolved list of variable indices to output
  std::vector< CFuint > m_varIDs;

  /// number of spatial dimensions
  CFuint m_dim;

  /// number of output values per state
  CFuint m_stride;

  /// output magnitude instead of components
  bool m_outputMagnitude;

}; // class OutputGradients

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_OutputGradients_hh
