#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_ComputeErrorMMSTNEQ_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_ComputeErrorMMSTNEQ_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionNEQ/Freno2021ManufacturedSolution.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Error of the FR solution against the manufactured solution of Freno et al. (2021), Eq. (19) (see
 * Freno2021ManufacturedSolution), for the two-temperature update variables RhoivtTv or LogRhoivLogTTv on quads.
 * Use it as ComputeErrorCom: it runs at every non-linear residual evaluation, on the state the residual is computed
 * from. Over the parallel updatable cells:
 * - quadrature with n x n Gauss-Legendre points per cell (n = QuadPointsPerDir, default P+3), physical coordinates
 *   from the geometric map of the cell, weights times the Jacobian determinant of the map;
 * - solution at a quadrature point: the update-variable polynomial of the cell evaluated there, the logarithmic
 *   variables (ln rho_s, ln T, ln Tv) then exponentiated;
 * - fields rho_s, u, v, T, Tv; e = (solution - manufactured)/bar value (v divided by ub);
 * - L1 = integral |e|, L2 = sqrt(integral e^2), Linf = max |e| over the quadrature points, with the coordinates of
 *   the point where the max is reached; the integral of 1 (domain area) is written as a check;
 * - residual R at the solution points (the rhs holds R times the Jacobian determinant): per equation
 *   max |R| and the integral of |R| with the Gauss-Legendre weights of the solution points.
 * Rank 0 appends one line per evaluation to OutputFile. Needs a dimensional 2D run (reference values 1, reference
 * length 1) and Gauss-Legendre solution points.
 *
 * @author Rayan Dhib
 */
class ComputeErrorMMSTNEQ : public FluxReconstructionSolverCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  ComputeErrorMMSTNEQ(const std::string& name);

  /// Destructor
  ~ComputeErrorMMSTNEQ();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

  /// Set up private data and data
  void setup();

  /// Compute the errors and write them
  void execute();

  /**
   * Gauss-Legendre nodes and weights on [-1, 1]
   * @param n number of points
   * @param x nodes, ascending
   * @param w weights
   */
  static void gaussLegendre(const CFuint n, std::vector<CFreal>& x, std::vector<CFreal>& w);

protected: // data

  /// socket for the residual
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// option: species names in the order of the chemistry library
  std::vector<std::string> m_speciesNames;

  /// option: bar partial densities [kg/m3], same order as the names
  std::vector<CFreal> m_rhoBar;

  /// option: velocity scale ub [m/s]
  CFreal m_uBar;

  /// option: temperature scale Tb [K]
  CFreal m_TBar;

  /// option: vibrational temperature scale Tvb [K]
  CFreal m_TvBar;

  /// option: amplitude e of the perturbations
  CFreal m_eps;

  /// option: name of the output file
  std::string m_outputFile;

  /// option: number of quadrature points per direction (0: P+3)
  CFuint m_nbQuadPerDir;

  /// manufactured solution
  Freno2021ManufacturedSolution m_mms;

  /// number of equations
  CFuint m_nbrEqs;

  /// number of species
  CFuint m_nbSpecies;

  /// true if the update variables are logarithmic (ln rho_s, u, v, ln T, ln Tv)
  bool m_logVariables;

  /// local coordinates of the quadrature points
  std::vector< RealVector > m_quadLocalCoords;

  /// quadrature weights on [-1, 1]^2
  std::vector< CFreal > m_quadWeights;

  /// values of the solution basis polynomials at the quadrature points
  std::vector< std::vector< CFreal > > m_solPolyValsAtQuad;

  /// Gauss-Legendre weights of the solution points on [-1, 1]^2, in the solution point order
  std::vector< CFreal > m_solWeights;

  /// number of evaluations so far
  CFuint m_nbCalls;

}; // class ComputeErrorMMSTNEQ

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluxReconstructionMethod_ComputeErrorMMSTNEQ_hh
